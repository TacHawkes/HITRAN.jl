# Type piracy free as of 0.1.3 to add the false return if key is 'missing'
# This definition of haskey only lives in the HITRAN module and does not monkey-patch
# Base.haskey anymore
haskey(collection, key) = Base.haskey(collection, key)
haskey(line::SQLite.Row, key::Missing) = false

Sij_T(S_ij_ref, T, T_ref, Q_T, Q_T_ref, e_lower, ν_ij)::Float64 = S_ij_ref * Q_T_ref / Q_T * exp(-c_c2 * e_lower / T) * (1 - exp(-c_c2 * ν_ij / T)) / exp(-c_c2 * e_lower / T_ref) * (1 - exp(-c_c2 * ν_ij / T_ref))

γ_Doppler(T, ν_ij, M)::Float64 = ν_ij / c_c_SI * √(2c_NA_SI * c_kB_SI * T * log(2) / (M * 1e-3))

γ_collision_0(γ0_ref, T, T_ref, p, p_ref, n_diluent)::Float64 = γ0_ref * p / p_ref * (T_ref / T)^n_diluent

function get_components(tables::AbstractVector{String})
    # use all components in the tables and get their natural abundance
    subqueries = []
    for table in tables
        push!(subqueries, " 
            SELECT      DISTINCT global_iso_id AS gid,
                        iso.molecule_id AS molecule_id,
                        iso.local_id AS local_id,
                        iso.abundance AS abundance 
            FROM        " * table * " 
            LEFT JOIN   isotopologues iso 
            ON          (global_id=gid) ")
    end
    sql = join(subqueries, " UNION ALL ")
    result = query_local_db(sql)        
    components = Dict{Tuple{Integer,Integer},AbstractFloat}()
    for row in result
        components[(row.molecule_id, row.local_id)] = row.abundance
    end            

    return components, components
end

function get_components(comps::AbstractVector{Tuple{T,T}}) where T <: Integer
    # load the natural abundances
    sql = " SELECT  molecule_id, local_id, abundance 
            FROM    isotopologues
            WHERE 	(molecule_id, local_id)
            IN 		(VALUES " * join(["(" * join("?"^length(t), ',') * ")" for t in comps], ',') * ")"
    result = query_local_db(sql, [i for t in comps for i in t])
    # replace components by new dict
    components = Dict{Tuple{Integer,Integer},AbstractFloat}()
    for row in result
        components[(row.molecule_id, row.local_id)] = row.abundance
    end        
    
    return components, components
end

function get_components(comps::Dict{Tuple{T,T},V}) where T <: Integer where V <: AbstractFloat
    # load the natural abundances
    sql = " SELECT  molecule_id, local_id, abundance 
            FROM    isotopologues
            WHERE 	(molecule_id, local_id)
            IN 		(VALUES " * join(["(" * join("?"^length(t), ',') * ")" for t in keys(comps)], ',') * ")"
    result = query_local_db(sql, [i for t in keys(comps) for i in t])
    # replace components by new dict
    natural_abundances = Dict{Tuple{Integer,Integer},AbstractFloat}()
    for row in result
        natural_abundances[(row.molecule_id, row.local_id)] = row.abundance
    end

    return comps, natural_abundances
end

function get_diluents(diluent::Dict{Symbol, T}, components) where T <: Number
    if length(components) > 1
        throw(ErrorException(
            "For gas mixtures a diluent entry for every component has to be supplied"
            ))
    end

    if (0 <= sum(values(diluent)) <= 1) == false
        throw(ErrorException(
            "Sum of diluent fractions must not exceed 1 or lie below zero"
            ))    
   end

    # convert to different dict format
    return Dict(
        first(keys(components)) => diluent
    )
end

function get_diluents(diluent::Dict{Tuple{T, T}, Dict{Symbol, V}}, components) where T <: Integer where V <: Number
    if length(diluent) > length(components)
        throw(ErrorException("More diluents than components specified."))
    end

    diluents = Dict()

    # water abundance
    has_water = false
    water_abundance = 0.0    
    if 1 in getindex.(keys(components), 1)
        has_water = true

        for (k,v) in components
            if k[1] == 1
                water_abundance += v
            end
        end
    end    

    for (mi_tuple, abundance) in components
        if mi_tuple ∈ keys(diluent)
            # insert provided values
            diluents[mi_tuple] = diluent[mi_tuple]
        else
            # default values and handle H2O            
            if has_water
                diluents[mi_tuple] = Dict(
                    :self => abundance - water_abundance,
                    :air => 1 - abundance - water_abundance,
                    :H2O => water_abundance
                )
            else                
                diluents[mi_tuple] = Dict(
                    :self => abundance,
                    :air => 1 - abundance
                )
            end
        end

        if (0 <= sum(values(diluents[mi_tuple])) <= 1) == false
            throw(ErrorException(
            "Sum of diluent fractions must not exceed 1 or lie below zero"
            ))
        end
    end

    return diluents
end

function get_diluent_names(diluent)
    # get all diluent names
    names = []
    for (comp_mi, dil) in diluent
        for k in keys(dil)
            push!(names, k)
        end
    end

    return names
end

function parse_kwargs(tables;kwargs...)    
    intensity_threshold = get(kwargs, :intensity_threshold, -Inf)
    pressure = convert(Float64, get(kwargs, :pressure, c_p_ref))
    temperature = convert(Float64, get(kwargs, :temperature, c_T_ref))

    ν_range = get(kwargs, :ν_range) do
        # find limits from the used tables
        subqueries = []
        for table in tables
            push!(subqueries, " SELECT nu FROM " * table)
        end
        subquery = join(subqueries, " UNION ALL ")

        sql = "SELECT MIN(nu) AS min_nu, MAX(nu) AS max_nu FROM (" * subquery * ")"        
        result = first(query_local_db(sql))
        ν_range = (result.min_nu, result.max_nu)
    end
    ν_min, ν_max = ν_range    
    ν_step = get(kwargs, :ν_step, 0.01)
    ν_wing = get(kwargs, :ν_wing, 0.)
    ν_wing_hw = get(kwargs, :ν_wing_hw, 50.)    

    components, natural_abundances = get_components(get(kwargs, :components, tables))        
    diluent = get(kwargs, :diluent, nothing)
    if isnothing(diluent)
        diluent = get_diluents(Dict((0,0) => Dict(:self => 1.0)), components)
    else
        diluent = get_diluents(diluent, components)
    end    

    return intensity_threshold, pressure, temperature, ν_range, 
        ν_min, ν_max, ν_step, ν_wing, ν_wing_hw, diluent, components, natural_abundances
end

"""
    α(tables::AbstractVector{String} [, profile=:hartmann_tran; kwargs...])

Computes the absorption coefficient using line-by-line data stored in the
database tables specified in `tables`. The lineshape can be optionally specified
using the `profile` argument and one of the Symbol keys `:hartmann_tran`, `:voigt`, `:sdvoigt`,
`:lorentz`, `:gauss`. If no keyword arguments are specified, they will be automatically
chosen from the tables provided.

# Keyword arguments
- `components`: the components of the gas mixture for the calculation. Can be either a vector of tuples with `(molecule_id, local_iso_id)`
                or a `Dict` with the `(molecule_id, local_iso_id)` tuple as key and the abundance as value. If the vector of tuples is supplied
                the natural abundance will be used, so this makes no sense for gas mixtures other than isotopologues of the same molecule.                
- `intensity_threshold`: the minimum line strength in ``cm^{-1}/(\\text{molecule} \\cdot cm^{-2})``
- `pressure`: the environmental pressure in atmospheres (default: $c_p_ref atm)
- `temperature`: the environmental temperature in Kelvin (default: $c_T_ref K)
- `ν_range`: a tuple of the form (ν_min, ν_max) where ν_min/ν_max is the minimum/maximum wavenumber for the absorption_spectrum respectively in ``cm^{-1}``
- `ν_step`: the wavenumber step in ``cm^{-1}`` (default: 0.01 ``cm^{-1}``)
- `ν_wing`: absolute calculation width of a line in ``cm^{-1}`` (default: 0 ``cm^{-1}``)
- `ν_wing_hw`: relative calculation width of a line in multiples of a half-width (default: 50)
- `diluent`: a `Dict` of the diluting substances, specified as `Symbol`, e.g. `:air` or `:H2O` for the key and the relative concentration as key (default: `Dict(:self => 1.0)`). See the example for details.
"""
function α(tables::AbstractVector{String}, profile=:hartmann_tran;kwargs...)        
    # parse and prepare input arguments
    intensity_threshold, pressure, temperature, ν_range, ν_min, ν_max,
    ν_step, ν_wing, ν_wing_hw, diluent, components, natural_abundances = parse_kwargs(tables;kwargs...)    

    # prepare the common fields and parameters    
    ν = ν_range[1]:ν_step:ν_range[2]
    len = length(ν)
    data = zeros(Float64, len)    
    data_cache = zeros(Float64, len)    
    molar_masses = Dict(
        keys(components) .=> molar_mass(keys(components))
    )    

    # molecule concentration in molecules / cm^3
    # allows conversion from cm^2/molecule -> cm^-1
    volume_conc = pressure * c_atm_factor / c_kB_SI / temperature * 1e-6
    
    # iterate over all tables and calculate line profiles line-by-line
    for table in tables
        # get line data (note: as :standard is enforced when fetching HITRAN
        # source data, it is assumed that at least all these fields exist and
        # no further check is performed)   
        result = query_local_db(
            "SELECT     * 
            FROM        " * table * "
            WHERE 	    (molec_id, local_iso_id)
            IN 		    (VALUES " * join(["(" * join("?"^length(t), ',') * ")" for t in keys(components)], ',') * ")
            AND         nu >= ?
            AND         nu <= ?", vcat([i for t in keys(components) for i in t], [ν_range...]))

        # call profile specific generic preparation function
        if haskey(profile_preparation_map, profile)
            profile_kwargs = profile_preparation_map[profile](; temperature, diluent, column_names=result.names)            
        else        
            profile_kwargs = Dict()
        end        
        
        q_cache = Dict{Integer,AbstractFloat}()
        qref_cache = Dict{Integer,AbstractFloat}()        
        for line in result            
            # partition sum
            q_t = get!(q_cache, line.global_iso_id) do 
                tips(line.global_iso_id, temperature) 
            end
            q_t_ref = get!(qref_cache, line.global_iso_id) do 
                tips(line.global_iso_id, c_T_ref) 
            end

            # environment adjusted line intensity            
            S_ij = Sij_T(line.sw, temperature, c_T_ref, q_t, q_t_ref, line.elower, line.nu)            

            # skip weak lines
            if S_ij < intensity_threshold
                continue
            end         
            
            factor = volume_conc * S_ij * components[(line.molec_id, line.local_iso_id)] / 
                natural_abundances[(line.molec_id, line.local_iso_id)]
            
            # call lineshape function
            lineshape_map[profile](line,
                diluent[(line.molec_id, line.local_iso_id)],
                temperature,
                pressure,
                molar_masses[(line.molec_id, line.local_iso_id)],
                ν,
                ν_wing,
                ν_wing_hw,
                factor,
                data,
                data_cache;
                profile_kwargs...)            
        end        
    end

    return ν, data
end

function hartmann_tran_reference_temperature(T)::Float64
    for (T_range, T_ref) in c_HT_T_ref
        if T >= T_range[1] && T < T_range[2]
            return T_ref
        end
    end
end

function hartmann_tran_profile!(
    out::AbstractVector{T},
    ν::AbstractVector{T},
    ν_0::T,
    ν_VC::V,
    γ_D::T,
    γ_0::T,
    γ_2::T,
    Δ_0::T,
    Δ_2::T,
    η::V
) where T <: AbstractFloat where V <: Complex    
    C_0 = γ_0 + im * Δ_0
    C_2 = γ_2 + im * Δ_2
    C_0t = (1 - η) * (C_0 - 3C_2 / 2.) + ν_VC
    C_2t = (1 - η) * C_2
    ν_a0 = c_c_SI * γ_D / (√(log(2)) * ν_0)
        
    if (abs(C_2t) ≈ 0.)                
        Threads.@threads for i = 1:length(ν)
            Z_m = (im * (ν_0 - ν[i]) + C_0t) / (ν_0 * ν_a0 / c_c_SI)
            wofz_Z_m = wofz(im * Z_m)      
            A_ν = √(π) * c_c_SI / ν_0 / ν_a0 * wofz_Z_m       
            B_ν = √(π) * c_c_SI * ν_a0 / ν_0 * ((1 - Z_m^2) * wofz_Z_m + Z_m / √(π))
            @inbounds out[i] = 1 / π * real(A_ν / (1 - (ν_VC - η * (C_0 - 3C_2 / 2) * A_ν + (η * C_2 / ν_a0^2) * B_ν)))
        end        
    else
        Y = (ν_0 * ν_a0 / 2 / c_c_SI / C_2t)^2
        Threads.@threads for i = 1:length(ν)
            X = (im * (ν_0 - ν[i]) + C_0t) / C_2t        
            Z_p = √(X + Y) + √(Y) 
            Z_m = √(X + Y) - √(Y)

            wofz_Z_m = wofz(im * Z_m)      
            wofz_Z_p = wofz(im * Z_p)      

            A_ν = √(π) * c_c_SI / ν_0 / ν_a0 * (wofz_Z_m - wofz_Z_p)
            B_ν = ν_a0^2 / C_2t^2 * (-1 + √(π) / 2 / √(Y) * (1 - Z_m^2) * wofz_Z_m - √(π) / 2 / √(Y) * (1 - Z_p^2) * wofz_Z_p)        
            @inbounds out[i] = 1 / π * real(A_ν / (1 - (ν_VC - η * (C_0 - 3C_2 / 2) * A_ν + (η * C_2 / ν_a0^2) * B_ν)))
        end
    end
end

function get_line_par(line::SQLite.Row, key::Symbol, ::Type{T}=Float64) where T
    q = SQLite.getquery(line)    
    ind = q.lookup[key]
    if q.types[ind] <: Union{Missing, T}
        SQLite.getvalue(q, ind, T)
    else
        zero(T)
    end  
end

function get_line_par(line::SQLite.Row, key::Symbol, fallback_key::Symbol, ::Type{T}=Float64) where T 
    q = SQLite.getquery(line)    
    ind = q.lookup[key]
    if q.types[ind] <: Union{Missing, T}
        SQLite.getvalue(q, ind, T)
    elseif q.types[q.lookup[fallback_key]] <: Union{Missing, T}
        SQLite.getvalue(q, q.lookup[fallback_key], T)
    else
        zero(T)
    end    
end
get_line_par(line::SQLite.Row, key::Symbol, ::Missing, ::Type{T}=Float64) where T = get_line_par(line, key, T)
get_line_par(::SQLite.Row, ::Missing)::Float64 = c_default_zero
get_line_par(::SQLite.Row, ::Missing, ::Missing)::Float64 = c_default_zero

function hartmann_tran_lineshape(
    line            :: SQLite.Row,
    diluent         :: Dict{Symbol,T},
    temperature     :: T,
    pressure        :: T,
    mass            :: T,        
    ν               :: AbstractRange,
    ν_wing          :: T,
    ν_wing_hw       :: T,
    factor          :: T,
    data            :: AbstractVector{T},
    out_cache       :: AbstractVector{T};
    kwargs...
) where T <: AbstractFloat
    T_ref_HT = get(kwargs, :T_ref_HT, c_T_ref)
    fields = kwargs[:fields]
    voigt_fields = kwargs[:voigt_fields]
    ν_0 = get_line_par(line, :nu)
    γ_D = γ_Doppler(temperature, ν_0, mass)    
    γ_0 = γ_2 = Δ_0 = Δ_2 = zero(Float64)
    ν_VC = zero(ComplexF64)
    η = zero(ComplexF64)

    # loop over all diluents and build combined line parameters
    for (diluent_name, diluent_abundance) in diluent    
        # get Hartmann-Tran or Voigt parameters if available        

        # γ_0 contribution
        γ_0_dil = get_line_par(line, fields[diluent_name][:γ_0], voigt_fields[diluent_name][:γ_0])
        n_dil = get_line_par(line, fields[diluent_name][:n], voigt_fields[diluent_name][:n])        
        if (diluent_name == :self || n_dil == c_default_zero)
            n_dil = get_line_par(line, :n_air)
        end        
        T_ref = (haskey(line, fields[diluent_name][:γ_0]) && haskey(line, fields[diluent_name][:n])) ? T_ref_HT : c_T_ref        
        γ_0t = γ_collision_0(γ_0_dil, temperature, T_ref, pressure, 
            c_p_ref, n_dil)
        γ_0 += diluent_abundance * γ_0t        

        # Δ_0 contribution        
        Δ_0_dil = get_line_par(line, fields[diluent_name][:Δ_0], voigt_fields[diluent_name][:Δ_0])
        Δ_0p_dil = get_line_par(line, fields[diluent_name][:Δ_0p], voigt_fields[diluent_name][:Δ_0p])
        T_ref = (haskey(line, fields[diluent_name][:Δ_0]) && haskey(line, fields[diluent_name][:Δ_0p])) ? T_ref_HT : c_T_ref        
        Δ_0t = Δ_0_dil + Δ_0p_dil * (temperature - T_ref) * pressure / c_p_ref
        Δ_0 += diluent_abundance * Δ_0t

        # γ_2 contribution                
        γ_2_dil = get_line_par(line, fields[diluent_name][:γ_2])
        γ_2t = γ_2_dil * pressure / c_p_ref
        γ_2 += diluent_abundance * γ_2t
        
        #  Δ_2 contribution                
        Δ_2_dil = get_line_par(line, fields[diluent_name][:Δ_2])
        Δ_2t = Δ_2_dil * pressure / c_p_ref
        Δ_2 += diluent_abundance * Δ_2t

        # η contribution        
        η_dil = get_line_par(line, fields[diluent_name][:η])
        η += diluent_abundance * η_dil * (γ_2t - im * Δ_2t)

        # ν_VC contribution        
        ν_VC_dil = get_line_par(line, fields[diluent_name][:ν_VC])
        κ_dil = get_line_par(line, fields[diluent_name][:κ])
        ν_VC += diluent_abundance * ν_VC_dil * (T_ref / temperature)^κ_dil * pressure
        ν_VC -= η_dil * diluent_abundance * (γ_0t - im * Δ_0t)
    end
    
    if η != 0.
        η /= γ_2 - im * Δ_2        
    end

    ν_VC += η * (γ_0 - im * Δ_0)    
    # use absolute or hw wing specification?
    ν_wing_val = max(ν_wing, ν_wing_hw * γ_0, ν_wing_hw * γ_D)

    # find a suitable bisection for integrating the data into the vector
    ind_lo = searchsortedfirst(ν, ν_0 - ν_wing_val)
    ind_hi = searchsortedlast(ν, ν_0 + ν_wing_val)    
    
    hartmann_tran_profile!(out_cache, @view(ν[ind_lo:ind_hi]), ν_0, ν_VC, γ_D, γ_0, γ_2, Δ_0, Δ_2, η)
    for i = 1:(ind_hi - ind_lo + 1)
        data[ind_lo + i - 1] += factor * out_cache[i]
    end    
end

get_symbol(cn::AbstractVector{Symbol}, s1::Symbol, s2::Symbol)::Union{Symbol,Missing} = s1 in cn ? s1 : s2 in cn ? s2 : missing
get_symbol(cn::AbstractVector{Symbol}, s1::Symbol)::Union{Symbol,Missing} = s1 in cn ? s1 : missing

function prepare_hartmann_tran_kwargs(;kwargs...)
    T = get(kwargs, :temperature, c_T_ref)
    diluent = get(kwargs, :diluent) do 
        Dict()
    end
    T_ht = hartmann_tran_reference_temperature(T)
    column_names = get(kwargs, :column_names) do
        Dict()
    end

    # prepare field names for diluents
    fields = Dict{Symbol,Dict{Symbol,Union{Symbol,Missing}}}()

    # get all diluent names    
    for diluent_name in get_diluent_names(diluent)
        fields[diluent_name] = Dict(
            :γ_0 => get_symbol(column_names, Symbol(@sprintf("gamma_HT_0_%s_%d", diluent_name, T_ht)), Symbol(@sprintf("gamma_%s", diluent_name))),
            :n => get_symbol(column_names, Symbol(@sprintf("n_HT_0_%s_%d", diluent_name, T_ht)), Symbol(@sprintf("n_%s", diluent_name))),            
            :Δ_0 => get_symbol(column_names, Symbol(@sprintf("delta_HT_0_%s_%d", diluent_name, T_ht)), Symbol(@sprintf("delta_%s", diluent_name))),            
            :Δ_0p => get_symbol(column_names, Symbol(@sprintf("deltap_HT_0_%s_%d", diluent_name, T_ht)), Symbol(@sprintf("deltap_%s", diluent_name))),            
            :γ_2 => get_symbol(column_names, Symbol(@sprintf("gamma_HT_2_%s_%d", diluent_name, T_ht))),            
            :Δ_2 => get_symbol(column_names, Symbol(@sprintf("delta_HT_2_%s_%d", diluent_name, T_ht))),
            :η => get_symbol(column_names, Symbol(@sprintf("eta_HT_%s", diluent_name))),            
            :ν_VC => get_symbol(column_names, Symbol(@sprintf("nu_HT_%s", diluent_name))),            
            :κ => get_symbol(column_names, Symbol(@sprintf("kappa_HT_%s", diluent_name)))            
        )
    end
    voigt_args = prepare_voigt_kwargs(;kwargs...)
    return Dict(
        :T_ref_HT => T_ht,
        :fields => fields,     
        :voigt_fields => voigt_args[:fields]
    )
end

voigt_profile!(
    out::AbstractVector{T},
    ν::AbstractVector{T},
    ν_0::T,    
    γ_D::T,
    γ_0::T
) where T <: AbstractFloat = hartmann_tran_profile!(out, ν, ν_0, 0. + im * 0., γ_D, γ_0, 0., 0., 0., 0. + im * 0.)

function voigt_lineshape(
    line            :: SQLite.Row,
    diluent         :: Dict{Symbol,T},
    temperature     :: T,
    pressure        :: T,
    mass            :: T,        
    ν               :: AbstractRange,
    ν_wing          :: T,
    ν_wing_hw       :: T,
    factor          :: T,
    data            :: AbstractVector{T},
    out_cache       :: AbstractVector{T};
    kwargs...
) where T <: AbstractFloat 
    fields = kwargs[:fields]    
    # initialize lineshape specific parameters    
    ν_0 = get_line_par(line, :nu)
    γ_D = γ_Doppler(temperature, ν_0, mass)
    γ_0 = Δ_0 = Float64(0.0)    

    # loop over all diluents and build combined line parameters    
    for (diluent_name, diluent_abundance) in diluent    
        # get Voigt parameters if available        

        # γ_0 contribution
        γ_0_dil = get_line_par(line, fields[diluent_name][:γ_0])        
        n_dil::Float64 = get_line_par(line, fields[diluent_name][:n])        
        if (diluent_name == :self || n_dil == c_default_zero)
            n_dil = get_line_par(line, :n_air)
        end                
        γ_0t = γ_collision_0(γ_0_dil, temperature, c_T_ref, pressure, 
            c_p_ref, n_dil)
        γ_0 += diluent_abundance * γ_0t

        # Δ_0 contribution        
        Δ_0_dil = get_line_par(line, fields[diluent_name][:Δ_0])
        Δ_0p_dil = get_line_par(line, fields[diluent_name][:Δ_0p])                
        Δ_0t = Δ_0_dil + Δ_0p_dil * (temperature - c_T_ref) * pressure / c_p_ref
        Δ_0 += diluent_abundance * Δ_0t       
    end
        
    # use absolute or hw wing specification?
    ν_wing_val::Float64 = max(ν_wing, ν_wing_hw * γ_0, ν_wing_hw * γ_D)

    # find a suitable bisection for integrating the data into the vector
    ind_lo = searchsortedfirst(ν, ν_0 - ν_wing_val)
    ind_hi = searchsortedlast(ν, ν_0 + ν_wing_val)  

    voigt_profile!(out_cache, @view(ν[ind_lo:ind_hi]), ν_0 + Δ_0, γ_D, γ_0)
    for i = 1:(ind_hi - ind_lo + 1)
        data[ind_lo + i - 1] += factor * out_cache[i]
    end
end

function prepare_voigt_kwargs(;kwargs...)        
    diluent = get(kwargs, :diluent) do 
        Dict()
    end    
    column_names = get(kwargs, :column_names) do
        Dict()
    end     

    # prepare field names for diluents
    fields = Dict{Symbol,Dict{Symbol,Union{Symbol,Missing}}}()
    for diluent_name in get_diluent_names(diluent)
        fields[diluent_name] = Dict(            
            :γ_0 => get_symbol(column_names, Symbol(@sprintf("gamma_%s", diluent_name))),
            :n => get_symbol(column_names, Symbol(@sprintf("n_%s", diluent_name))),
            :Δ_0 => get_symbol(column_names, Symbol(@sprintf("delta_%s", diluent_name))), 
            :Δ_0p => get_symbol(column_names, Symbol(@sprintf("deltap_%s", diluent_name)))
        )
    end

    return Dict(        
        :fields => fields        
    )
end

speed_dependent_voigt_profile!(
    out::AbstractVector{T},
    ν::AbstractVector{T},
    ν_0::T,    
    γ_D::T,
    γ_0::T,
    γ_2::T,
    Δ_0::T       
 ) where T <: AbstractFloat = hartmann_tran_profile!(out, ν, ν_0, 0.0 + 0.0*im, γ_D, γ_0, γ_2, Δ_0, 0., 0.0 + 0.0*im)

 function speed_dependent_voigt_lineshape(
    line            :: SQLite.Row,
    diluent         :: Dict{Symbol,T},
    temperature     :: T,
    pressure        :: T,
    mass            :: T,        
    ν               :: AbstractRange,
    ν_wing          :: T,
    ν_wing_hw       :: T,
    factor          :: T,
    data            :: AbstractVector{T},
    out_cache       :: AbstractVector{T};
    kwargs...
) where T <: AbstractFloat 
    fields = kwargs[:fields]    
    # initialize lineshape specific parameters    
    ν_0 = get_line_par(line, :nu)
    γ_D = γ_Doppler(temperature, ν_0, mass)
    γ_0 = Δ_0 = γ_2 = Float64(0.0)    

    # loop over all diluents and build combined line parameters    
    for (diluent_name, diluent_abundance) in diluent    
        # get Voigt parameters if available        

        # γ_0 contribution
        γ_0_dil = get_line_par(line, fields[diluent_name][:γ_0])        
        n_dil::Float64 = get_line_par(line, fields[diluent_name][:n])        
        if (diluent_name == :self || n_dil == c_default_zero)
            n_dil = get_line_par(line, :n_air)
        end                
        γ_0t = γ_collision_0(γ_0_dil, temperature, c_T_ref, pressure, 
            c_p_ref, n_dil)
        γ_0 += diluent_abundance * γ_0t

        # Δ_0 contribution        
        Δ_0_dil = get_line_par(line, fields[diluent_name][:Δ_0])
        Δ_0p_dil = get_line_par(line, fields[diluent_name][:Δ_0p])                
        Δ_0t = Δ_0_dil + Δ_0p_dil * (temperature - c_T_ref) * pressure / c_p_ref
        Δ_0 += diluent_abundance * Δ_0t    
        
        # γ_2 contribution
        γ_2_dil = get_line_par(line, fields[diluent_name][:γ_2])
        γ_2t = γ_2_dil * pressure / c_p_ref
        γ_2 += diluent_abundance * γ_2t
    end
        
    # use absolute or hw wing specification?
    ν_wing_val::Float64 = max(ν_wing, ν_wing_hw * γ_0, ν_wing_hw * γ_D)

    # find a suitable bisection for integrating the data into the vector
    ind_lo = searchsortedfirst(ν, ν_0 - ν_wing_val)
    ind_hi = searchsortedlast(ν, ν_0 + ν_wing_val)  

    speed_dependent_voigt_profile!(out_cache, @view(ν[ind_lo:ind_hi]), ν_0, γ_D, γ_0, γ_2, Δ_0)
    for i = 1:(ind_hi - ind_lo + 1)
        data[ind_lo + i - 1] += factor * out_cache[i]
    end
end

 function prepare_speed_dependent_voigt_kwargs(;kwargs...)        
    diluent = get(kwargs, :diluent) do 
        Dict()
    end    
    column_names = get(kwargs, :column_names) do
        Dict()
    end     

    # prepare field names for diluents
    fields = Dict{Symbol,Dict{Symbol,Union{Symbol,Missing}}}()
    for diluent_name in get_diluent_names(diluent)
        fields[diluent_name] = Dict(            
            :γ_0 => get_symbol(column_names, Symbol(@sprintf("gamma_%s", diluent_name))),
            :n => get_symbol(column_names, Symbol(@sprintf("n_%s", diluent_name))),
            :Δ_0 => get_symbol(column_names, Symbol(@sprintf("delta_%s", diluent_name))), 
            :Δ_0p => get_symbol(column_names, Symbol(@sprintf("deltap_%s", diluent_name))),
            :γ_2 => get_symbol(column_names, Symbol(@sprintf("sd_%s", diluent_name)))
        )
    end

    return Dict(        
        :fields => fields        
    )
end

function lorentz_profile!(
    out::AbstractVector{T},
    ν::AbstractVector{T},
    ν_0::T,
    γ_0::T
) where T <: AbstractFloat     
    for i = 1:length(ν)
        out[i] = γ_0 / (π * (γ_0^2 + (ν[i] - ν_0)^2))
    end
end

function lorentz_lineshape(
    line            :: SQLite.Row,
    diluent         :: Dict{Symbol,T},
    temperature     :: T,
    pressure        :: T,
    mass            :: T,        
    ν               :: AbstractRange,
    ν_wing          :: T,
    ν_wing_hw       :: T,
    factor          :: T,
    data            :: AbstractVector{T},
    out_cache       :: AbstractVector{T};
    kwargs...
) where T <: AbstractFloat 
    fields = kwargs[:fields]    
    # initialize lineshape specific parameters    
    ν_0 = get_line_par(line, :nu)    
    γ_0 = Δ_0 = Float64(0.0)        

    # loop over all diluents and build combined line parameters
    for (diluent_name, diluent_abundance) in diluent                        
        # γ_0 contribution
        γ_0_dil = get_line_par(line, fields[diluent_name][:γ_0])        
        n_dil::Float64 = get_line_par(line, fields[diluent_name][:n])        
        if (diluent_name == :self || n_dil == c_default_zero)
            n_dil = get_line_par(line, :n_air)
        end                
        γ_0t = γ_collision_0(γ_0_dil, temperature, c_T_ref, pressure, 
            c_p_ref, n_dil)
        γ_0 += diluent_abundance * γ_0t

        # Δ_0 contribution        
        Δ_0_dil = get_line_par(line, fields[diluent_name][:Δ_0])
        Δ_0p_dil = get_line_par(line, fields[diluent_name][:Δ_0p])                
        Δ_0t = Δ_0_dil + Δ_0p_dil * (temperature - c_T_ref) * pressure / c_p_ref
        Δ_0 += diluent_abundance * Δ_0t
    end

    # use absolute or hw wing specification?
    ν_wing_val::Float64 = max(ν_wing, ν_wing_hw * γ_0)

    # find a suitable bisection for integrating the data into the vector
    ind_lo = searchsortedfirst(ν, ν_0 - ν_wing_val)
    ind_hi = searchsortedlast(ν, ν_0 + ν_wing_val)  

    lorentz_profile!(out_cache, @view(ν[ind_lo:ind_hi]), ν_0 + Δ_0, γ_0)
    for i = 1:(ind_hi - ind_lo + 1)
        data[ind_lo + i - 1] += factor * out_cache[i]
    end
end

function gauss_profile!(
    out::AbstractVector{T},
    ν::AbstractVector{T},
    ν_0::T,
    γ_D::T    
) where T <: AbstractFloat
    for i = 1:length(ν)
        out[i] = √(log(2) / π) / γ_D * exp(-log(2) * ((ν[i] - ν_0) / γ_D)^2)  
    end
end

function gauss_lineshape(
    line            :: SQLite.Row,
    diluent         :: Dict{Symbol,T},
    temperature     :: T,
    pressure        :: T,
    mass            :: T,        
    ν               :: AbstractRange,
    ν_wing          :: T,
    ν_wing_hw       :: T,
    factor          :: T,
    data            :: AbstractVector{T},
    out_cache       :: AbstractVector{T};
    kwargs...
) where T <: AbstractFloat           
    # initialize lineshape specific parameters    
    ν_0 = get_line_par(line, :nu)
    γ_D = γ_Doppler(temperature, ν_0, mass)
    Δ_0 = get_line_par(line, :delta_air) * pressure / c_p_ref    

    # use absolute or hw wing specification?
    ν_wing_val = max(ν_wing, ν_wing_hw * γ_D)

    ind_lo = searchsortedfirst(ν, ν_0 - ν_wing_val)
    ind_hi = searchsortedlast(ν, ν_0 + ν_wing_val)

    gauss_profile!(out_cache, @view(ν[ind_lo:ind_hi]), ν_0 + Δ_0, γ_D)
    for i = 1:(ind_hi - ind_lo + 1)
        data[ind_lo + i - 1] += factor * out_cache[i]
    end
end

const profile_map = Dict(
    :gauss => gauss_profile!,
    :lorentz => lorentz_profile!,
    :voigt => voigt_profile!,    
    :sdvoigt => speed_dependent_voigt_profile!,    
    :hartmann_tran => hartmann_tran_profile!,        
)

const lineshape_map = Dict(    
    :hartmann_tran => hartmann_tran_lineshape,    
    :voigt => voigt_lineshape,
    :sdvoigt => speed_dependent_voigt_lineshape,
    :lorentz => lorentz_lineshape,
    :gauss => gauss_lineshape   
)

const profile_preparation_map = Dict(
    :hartmann_tran => prepare_hartmann_tran_kwargs,
    :voigt => prepare_voigt_kwargs,
    :sdvoigt => prepare_speed_dependent_voigt_kwargs,
    :lorentz => prepare_voigt_kwargs,    
)

# Fadeeva function
wofz(z) = erfcx(-im * z)