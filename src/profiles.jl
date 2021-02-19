
function get_line_parameter!(target, key, line::SQLite.Row, parameter::Symbol, fallback::Symbol=:none, default=c_default_zero, parameter_type::Type=AbstractFloat)
    if haskey(line, parameter) && typeof(line[parameter]) <: parameter_type
        setfield!(target, key, line[parameter])
    elseif haskey(line, fallback) && typeof(line[fallback]) <: parameter_type
        setfield!(target, key, line[fallback])
    else
        setfield!(target, key, default)                
    end
end

Sij_T(S_ij_ref, T, T_ref, Q_T, Q_T_ref, e_lower, ν_ij) = S_ij_ref * Q_T_ref / Q_T * exp(-c_c2 * e_lower / T) * (1 - exp(-c_c2 * ν_ij / T)) / exp(-c_c2 * e_lower / T_ref) * (1 - exp(-c_c2 * ν_ij / T_ref))

γ_Doppler(T, ν_ij, M) = ν_ij / c_c_SI * √(2c_NA_SI * c_kB_SI * T * log(2) / (M * 1e-3))

γ_collision_0(γ0_ref, T, T_ref, p, p_ref, n_diluent) = γ0_ref * p / p_ref * (T_ref / T)^n_diluent

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

function parse_kwargs(tables;kwargs...)    
    db = get(kwargs, :db, current_db())
    
    intensity_threshold = get(kwargs, :intensity_threshold, nothing)
    pressure = get(kwargs, :pressure, c_p_ref)    
    temperature = get(kwargs, :temperature, c_T_ref)    

    ν_range = get(kwargs, :ν_range, nothing)
    if ν_range === nothing
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
    ν_min = ν_range[1]
    ν_max = ν_range[2]
    ν_step = get(kwargs, :ν_step, 0.01)
    ν_wing = get(kwargs, :ν_wing, 0.)
    ν_wing_hw = get(kwargs, :ν_wing_hw, 50.)

    diluent = get(kwargs, :diluent, Dict(:self => 1.))
    !all(0 .<= values(diluent) .<= 1) && error("Diluent must not exceed 0 or 1 for individual fractions")    

    components, natural_abundances = get_components(get(kwargs, :components, tables))        

    return db, intensity_threshold, pressure, temperature, ν_range, 
        ν_min, ν_max, ν_step, ν_wing, ν_wing_hw, diluent, components, natural_abundances
end

"""
    α(tables::AbstractVector{String} [, profile=:hartmann_tran; kwargs...])

Computes the absorption coefficient using line-by-line data stored in the
database tables specified in `tables`. The lineshape can be optionally specified
using the `profile` argument and one of the Symbol keys `:hartmann_tran`, `:voigt`,
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
- `diluent`: a `Dict` of the diluting substances, specified as `Symbol`, e.g. `:air` or `:H2O` for the key and the relative concentration as key (default: `Dict(:self => 1.0)`)
"""
function α(tables::AbstractVector{String}, profile=:hartmann_tran;kwargs...)
    # parse and prepare input arguments
    db, intensity_threshold, pressure, temperature, ν_range, ν_min, ν_max,
    ν_step, ν_wing, ν_wing_hw, diluent, components, natural_abundances = parse_kwargs(tables;kwargs...)    

    # prepare the common fields and parameters    
    ν = ν_range[1]:ν_step:ν_range[2]
    len = length(ν)
    data = zeros(len)    
    data_cache = zeros(len)    
    molar_masses = Dict(
        keys(components) .=> molar_mass(keys(components))
    )

    # molecule concentration in molecules / cm^3
    # allows conversion from cm^2/molecule -> cm^-1
    volume_conc = pressure * c_atm_factor / c_kB_SI / temperature * 1e-6

    # call profile specific generic preparation function
    if haskey(profile_preparation_map, profile)
        profile_kwargs = profile_preparation_map[profile](; temperature, diluent)
    else        
        profile_kwargs = Dict()
    end
    
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
        parameters = result.names
        for line in result            
            # partition sum
            q_t = tips(line.global_iso_id, temperature)
            q_t_ref = tips(line.global_iso_id, c_T_ref)

            # environment adjusted line intensity
            S_ij = Sij_T(line.sw, temperature, c_T_ref, q_t, q_t_ref, line.elower, line.nu)

            # skip weak lines
            if intensity_threshold !== nothing && S_ij < intensity_threshold
                continue
            end

            γ_D = γ_Doppler(temperature, line.nu, molar_masses[(line.molec_id, line.local_iso_id)])
            
            # call lineshape function
            lineshape_map[profile](line,
                diluent,
                temperature,
                pressure,
                γ_D,
                ν,
                ν_wing,
                ν_wing_hw,
                volume_conc * S_ij * components[(line.molec_id, line.local_iso_id)] / 
                natural_abundances[(line.molec_id, line.local_iso_id)],
                data,
                data_cache;
                profile_kwargs...)            
        end
    end

    return ν, data
end

function hartmann_tran_reference_temperature(T)
    for (T_range, T_ref) in c_HT_T_ref
        if T > T_range[1] && T < T_range[2]
            return T_ref
        end
    end
end

mutable struct HartmannTranLineParameters{T <: AbstractFloat,V <: Complex}
    ν_0::T
    ν_VC::V
    γ_D::T
    γ_0::T
    γ_2::T
    Δ_0::T
    Δ_2::T
    η::V
end

mutable struct HartmannTranDiluentParameters{T <: AbstractFloat}        
    T_ref::T
    γ_D::T
    γ_0::T
    γ_0t::T
    n::T
    γ_2::T
    γ_2t::T
    Δ_0::T
    Δ_0p::T
    Δ_0t::T
    Δ_2::T
    Δ_2t::T
    ν_VC::T
    κ::T
    η::T
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
        for i = 1:length(ν)
            Z_m = (im * (ν_0 - ν[i]) + C_0t) / (ν_0 * ν_a0 / c_c_SI)
            wofz_Z_m = wofz(im * Z_m)      
            A_ν = √(π) * c_c_SI / ν_0 / ν_a0 * wofz_Z_m       
            B_ν = √(π) * c_c_SI * ν_a0 / ν_0 * ((1 - Z_m^2) * wofz_Z_m + Z_m / √(π))
            out[i] = 1 / π * real(A_ν / (1 - (ν_VC - η * (C_0 - 3C_2 / 2) * A_ν + (η * C_2 / ν_a0^2) * B_ν)))
        end        
    else
        Y = (ν_0 * ν_a0 / 2 / c_c_SI / C_2t)^2
        for i = 1:length(ν)
            X = (im * (ν_0 - ν[i]) + C_0t) / C_2t        
            Z_p = √(X + Y) + √(Y) 
            Z_m = √(X + Y) - √(Y)

            wofz_Z_m = wofz(im * Z_m)      
            wofz_Z_p = wofz(im * Z_p)      

            A_ν = √(π) * c_c_SI / ν_0 / ν_a0 * (wofz_Z_m - wofz_Z_p)
            B_ν = ν_a0^2 / C_2t^2 * (-1 + √(π) / 2 / √(Y) * (1 - Z_m^2) * wofz_Z_m - √(π) / 2 / √(Y) * (1 - Z_p^2) * wofz_Z_p)        
            out[i] = 1 / π * real(A_ν / (1 - (ν_VC - η * (C_0 - 3C_2 / 2) * A_ν + (η * C_2 / ν_a0^2) * B_ν)))
        end
    end
end

hartmann_tran_profile!(
    out::AbstractVector{T},
    ν::AbstractVector{T},
    line::HartmannTranLineParameters) where T <: AbstractFloat = hartmann_tran_profile!(out, ν, line.ν_0, line.ν_VC, line.γ_D, line.γ_0, line.γ_2, line.Δ_0, line.Δ_2, line.η)

function hartmann_tran_lineshape(
    line,
    diluent,
    temperature,
    pressure,
    γ_D,
    ν,
    ν_wing,
    ν_wing_hw,
    factor,
    data,
    out_cache;
    kwargs...
)
    T_ref_HT = get(kwargs, :T_ref_HT, c_T_ref)
    fields = kwargs[:fields]
    line_parameters::HartmannTranLineParameters = kwargs[:line_parameters]
    diluent_parameters::HartmannTranDiluentParameters = kwargs[:diluent_parameters]
    # initialize lineshape specific parameters    
    line_parameters.ν_0 = line.nu
    line_parameters.γ_D = γ_D
    line_parameters.γ_0 = c_default_zero
    line_parameters.γ_2 = c_default_zero
    line_parameters.Δ_0 = c_default_zero
    line_parameters.Δ_2 = c_default_zero
    line_parameters.η = c_default_zero
    line_parameters.ν_VC = c_default_zero

    # loop over all diluents and build combined line parameters
    for (diluent_name, diluent_abundance) in diluent
        # get Hartmann-Tran or Voigt parameters if available

        diluent_parameters.T_ref = T_ref_HT
        # γ_0 contribution              
        get_line_parameter!(diluent_parameters, :γ_0, line, fields[diluent_name][:γ_0], fields[diluent_name][:γ_0_voigt])                      
        get_line_parameter!(diluent_parameters, :n, line, fields[diluent_name][:n], fields[diluent_name][:n_voigt], line.n_air)
        if (diluent_name == "self" && diluent_parameters.n == c_default_zero)
            diluent_parameters.n = line.n_air
        end        
        diluent_parameters.T_ref = (haskey(line, fields[diluent_name][:γ_0]) && haskey(line, fields[diluent_name][:n])) ? T_ref_HT : c_T_ref        
        diluent_parameters.γ_0t = γ_collision_0(diluent_parameters.γ_0, temperature, diluent_parameters.T_ref, pressure, 
            c_p_ref, diluent_parameters.n)
        line_parameters.γ_0 += diluent_abundance * diluent_parameters.γ_0t

        # Δ_0 contribution        
        get_line_parameter!(diluent_parameters, :Δ_0, line, fields[diluent_name][:Δ_0], fields[diluent_name][:Δ_0_voigt])        
        get_line_parameter!(diluent_parameters, :Δ_0p, line, fields[diluent_name][:Δ_0p], fields[diluent_name][:Δ_0p_voigt])
        diluent_parameters.T_ref = (haskey(line, fields[diluent_name][:Δ_0]) && haskey(line, fields[diluent_name][:Δ_0p])) ? T_ref_HT : c_T_ref        
        diluent_parameters.Δ_0t = (diluent_parameters.Δ_0 + (diluent_parameters.Δ_0p) * (temperature - diluent_parameters.T_ref) * pressure / c_p_ref)
        line_parameters.Δ_0 += diluent_abundance * diluent_parameters.Δ_0t

        # γ_2 contribution        
        get_line_parameter!(diluent_parameters, :γ_2, line, fields[diluent_name][:γ_2], :none)
        diluent_parameters.γ_2t = diluent_parameters.γ_2 * pressure / c_p_ref
        line_parameters.γ_2 += diluent_abundance * diluent_parameters.γ_2t
        
        #  Δ_2 contribution
        get_line_parameter!(diluent_parameters, :Δ_2, line, fields[diluent_name][:Δ_2], :none)
        diluent_parameters.Δ_2t = diluent_parameters.Δ_2 * pressure / c_p_ref
        line_parameters.Δ_2 += diluent_abundance * diluent_parameters.Δ_2t

        # η contribution
        get_line_parameter!(diluent_parameters, :η, line, fields[diluent_name][:η], :none)
        line_parameters.η += diluent_abundance * diluent_parameters.η * (diluent_parameters.γ_2t - im * diluent_parameters.Δ_2t)

        # ν_VC contribution
        get_line_parameter!(diluent_parameters, :ν_VC, line, fields[diluent_name][:ν_VC], :none)
        get_line_parameter!(diluent_parameters, :κ, line, fields[diluent_name][:κ], :none)
        line_parameters.ν_VC += diluent_abundance * diluent_parameters.ν_VC * (diluent_parameters.T_ref / temperature)^diluent_parameters.κ * pressure
        line_parameters.ν_VC -= diluent_parameters.η * diluent_abundance * (diluent_parameters.γ_0t - im * diluent_parameters.Δ_0t)
    end
    
    if line_parameters.η != 0.
        line_parameters.η /= line_parameters.γ_2 - im * line_parameters.Δ_2        
    end

    line_parameters.ν_VC += line_parameters.η * (line_parameters.γ_0 - im * line_parameters.Δ_0)    
    # use absolute or hw wing specification?
    ν_wing_val = max(ν_wing, ν_wing_hw * line_parameters.γ_0, ν_wing_hw * γ_D)

    # find a suitable bisection for integrating the data into the vector
    ind_lo = searchsortedfirst(ν, line.nu - ν_wing_val)
    ind_hi = searchsortedlast(ν, line.nu + ν_wing_val)    
    
    hartmann_tran_profile!(out_cache, @view(ν[ind_lo:ind_hi]), line_parameters)
    for i = 1:(ind_hi - ind_lo + 1)
        data[ind_lo + i - 1] += factor * out_cache[i]
    end    
end

function prepare_hartmann_tran_kwargs(;kwargs...)
    T = get(kwargs, :temperature, c_T_ref)
    diluent = get(kwargs, :diluent, Dict())
    T_ht = hartmann_tran_reference_temperature(T)

    # prepare field names for diluents
    fields = Dict{Symbol,Dict{Symbol,Symbol}}()
    for diluent_name in keys(diluent)
        fields[diluent_name] = Dict(
            :γ_0 => Symbol(@sprintf("gamma_HT_0_%s_%d", diluent_name, T_ht)),
            :γ_0_voigt => Symbol(@sprintf("gamma_%s", diluent_name)),
            :n => Symbol(@sprintf("n_HT_0_%s_%d", diluent_name, T_ht)),
            :n_voigt => Symbol(@sprintf("n_%s", diluent_name)),
            :Δ_0 => Symbol(@sprintf("delta_HT_0_%s_%d", diluent_name, T_ht)),
            :Δ_0_voigt => Symbol(@sprintf("delta_%s", diluent_name)),
            :Δ_0p => Symbol(@sprintf("deltap_HT_0_%s_%d", diluent_name, T_ht)),
            :Δ_0p_voigt => Symbol(@sprintf("deltap_%s", diluent_name)),
            :γ_2 => Symbol(@sprintf("gamma_HT_2_%s_%d", diluent_name, T_ht)),
            :Δ_2 => Symbol(@sprintf("delta_HT_2_%s_%d", diluent_name, T_ht)),
            :η => Symbol(@sprintf("eta_HT_%s", diluent_name)),
            :ν_VC => Symbol(@sprintf("nu_HT_%s", diluent_name)),
            :κ => Symbol(@sprintf("kappa_HT_%s", diluent_name))
        )
    end

    return Dict(
        :T_ref_HT => T_ht,
        :fields => fields,
        :line_parameters => HartmannTranLineParameters(0., 0.0 * im, 0., 0., 0., 0., 0., 0.0 * im),
        :diluent_parameters => HartmannTranDiluentParameters(c_T_ref, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.0, 0., 0.0)
    )
end

voigt_profile!(
    out::AbstractVector{T},
    ν::AbstractVector{T},
    ν_0::T,    
    γ_D::T,
    γ_0::T
) where T <: AbstractFloat = hartmann_tran_profile!(out, ν, ν_0, 0. + im * 0., γ_D, γ_0, 0., 0., 0., 0. + im * 0.)

voigt_profile!(
    out::AbstractVector{T},
    ν::AbstractVector{T},
    line::HartmannTranLineParameters) where T <: AbstractFloat = voigt_profile!(out, ν, line.ν_0 + line.Δ_0, line.γ_D, line.γ_0)

function voigt_lineshape(
    line,
    diluent,
    temperature,
    pressure,
    γ_D,
    ν,
    ν_wing,
    ν_wing_hw,
    factor,
    data,
    out_cache;
    kwargs...
)    
    fields = kwargs[:fields]
    line_parameters::HartmannTranLineParameters = kwargs[:line_parameters]
    diluent_parameters::HartmannTranDiluentParameters = kwargs[:diluent_parameters]
    # initialize lineshape specific parameters    
    line_parameters.ν_0 = line.nu
    line_parameters.γ_D = γ_D
    line_parameters.γ_0 = c_default_zero
    line_parameters.Δ_0 = c_default_zero

    # loop over all diluents and build combined line parameters
    for (diluent_name, diluent_abundance) in diluent                        
        # γ_0 contribution        
        get_line_parameter!(diluent_parameters, :γ_0, line, fields[diluent_name][:γ_0])
        get_line_parameter!(diluent_parameters, :n, line, fields[diluent_name][:n], :none, line.n_air)
        if (diluent_name == "self" && diluent_parameters.n == c_default_zero)
            diluent_parameters.n = line.n_air
        end                        
        line_parameters.γ_0 += diluent_abundance * γ_collision_0(diluent_parameters.γ_0, temperature, 
                            c_T_ref, pressure, c_p_ref, diluent_parameters.n)

        # Δ_0 contribution
        get_line_parameter!(diluent_parameters, :Δ_0, line, fields[diluent_name][:Δ_0])
        get_line_parameter!(diluent_parameters, :Δ_0p, line, fields[diluent_name][:Δ_0p]) 
        line_parameters.Δ_0 += diluent_abundance * (diluent_parameters.Δ_0 + (diluent_parameters.Δ_0p) * (temperature - c_T_ref) * pressure / c_p_ref)                        
    end
    
    # use absolute or hw wing specification?
    ν_wing_val = max(ν_wing, ν_wing_hw * line_parameters.γ_0, ν_wing_hw * γ_D)

    ind_lo = searchsortedfirst(ν, line.nu - ν_wing_val)
    ind_hi = searchsortedlast(ν, line.nu + ν_wing_val)

    voigt_profile!(out_cache, @view(ν[ind_lo:ind_hi]), line_parameters)
    for i = 1:(ind_hi - ind_lo + 1)
        data[ind_lo + i - 1] += factor * out_cache[i]
    end
end

function prepare_voigt_kwargs(;kwargs...)    
    diluent = get(kwargs, :diluent, Dict())    

    # prepare field names for diluents
    fields = Dict{Symbol,Dict{Symbol,Symbol}}()
    for diluent_name in keys(diluent)
        fields[diluent_name] = Dict(            
            :γ_0 => Symbol(@sprintf("gamma_%s", diluent_name)),            
            :n => Symbol(@sprintf("n_%s", diluent_name)),  
            :Δ_0 => Symbol(@sprintf("delta_%s", diluent_name)),            
            :Δ_0p => Symbol(@sprintf("deltap_%s", diluent_name))            
        )
    end

    return Dict(        
        :fields => fields,
        :line_parameters => HartmannTranLineParameters(0., 0.0 * im, 0., 0., 0., 0., 0., 0.0 * im),
        :diluent_parameters => HartmannTranDiluentParameters(c_T_ref, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.0, 0., 0.0)
    )
end

speed_dependent_voigt_profile(
    ν::AbstractVector{T},
    ν_0::T,    
    γ_D::T,
    γ_0::T,
    γ_2::T,
    Δ_0::T,
    Δ_2::T    
) where T <: AbstractFloat = hartmann_tran_profile(ν, ν_0, 0., γ_D, γ_0, γ_2, Δ_0, Δ_2, 0.)

speed_dependent_rautian_profile(
    ν::AbstractVector{T},
    ν_0::T,
    ν_VC::T,
    γ_D::T,
    γ_0::T,
    γ_2::T,
    Δ_0::T,
    Δ_2::T
) where T <: AbstractFloat = hartmann_tran_profile(ν, ν_0, ν_VC, γ_D, γ_0, γ_2, Δ_0, Δ_2, 0.)

rautian_profile(
    ν::AbstractVector{T},
    ν_0::T,
    ν_VC::T,
    γ_D::T,
    γ_0::T,
    Δ_0::T,
) where T <: AbstractFloat = hartmann_tran_profile(ν, ν_0, ν_VC, γ_D, γ_0, 0., Δ_0, 0., 0.)

function lorentz_profile!(
    out::AbstractVector{T},
    ν::AbstractVector{T},
    ν_0::T,
    γ_0::T
) where T <: AbstractFloat     
    for i=1:length(ν)
        out[i] = γ_0 / (π * (γ_0^2 + (ν[i] - ν_0)^2))
    end
end

function lorentz_lineshape(
    line,
    diluent,
    temperature,
    pressure,
    γ_D,
    ν,
    ν_wing,
    ν_wing_hw,
    factor,
    data,
    out_cache;
    kwargs...
)    
    fields = kwargs[:fields]
    line_parameters::HartmannTranLineParameters = kwargs[:line_parameters]
    diluent_parameters::HartmannTranDiluentParameters = kwargs[:diluent_parameters]
    # initialize lineshape specific parameters    
    line_parameters.ν_0 = line.nu    
    line_parameters.γ_0 = c_default_zero
    line_parameters.Δ_0 = c_default_zero

    # loop over all diluents and build combined line parameters
    for (diluent_name, diluent_abundance) in diluent                        
        # γ_0 contribution        
        get_line_parameter!(diluent_parameters, :γ_0, line, fields[diluent_name][:γ_0])
        get_line_parameter!(diluent_parameters, :n, line, fields[diluent_name][:n], :none, line.n_air)
        if (diluent_name == "self" && diluent_parameters.n == c_default_zero)
            diluent_parameters.n = line.n_air
        end                        
        line_parameters.γ_0 += diluent_abundance * γ_collision_0(diluent_parameters.γ_0, temperature, 
                        c_T_ref, pressure, c_p_ref, diluent_parameters.n)

        # Δ_0 contribution
        get_line_parameter!(diluent_parameters, :Δ_0, line, fields[diluent_name][:Δ_0])
        get_line_parameter!(diluent_parameters, :Δ_0p, line, fields[diluent_name][:Δ_0p]) 
        line_parameters.Δ_0 += diluent_abundance * (diluent_parameters.Δ_0 + (diluent_parameters.Δ_0p) * (temperature - c_T_ref) * pressure / c_p_ref)                        
    end

    # use absolute or hw wing specification?
    ν_wing_val = max(ν_wing, ν_wing_hw * line_parameters.γ_0)

    ind_lo = searchsortedfirst(ν, line.nu - ν_wing_val)
    ind_hi = searchsortedlast(ν, line.nu + ν_wing_val)

    lorentz_profile!(out_cache, @view(ν[ind_lo:ind_hi]), line_parameters.ν_0 + line_parameters.Δ_0, line_parameters.γ_0)
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
    for i=1:length(ν)
        out[i] = √(log(2) / π) / γ_D * exp(-log(2) * ((ν[i] - ν_0) / γ_D)^2)  
    end
end

function gauss_lineshape(
    line,
    diluent,
    temperature,
    pressure,
    γ_D,
    ν,
    ν_wing,
    ν_wing_hw,
    factor,
    data,
    out_cache;
    kwargs...
)        
    line_parameters::HartmannTranLineParameters = kwargs[:line_parameters]
    # initialize lineshape specific parameters    
    line_parameters.ν_0 = line.nu        
    get_line_parameter!(line_parameters, :Δ_0, line, :delta_air) * pressure / c_p_ref  

    # use absolute or hw wing specification?
    ν_wing_val = max(ν_wing, ν_wing_hw * line_parameters.γ_D)

    ind_lo = searchsortedfirst(ν, line.nu - ν_wing_val)
    ind_hi = searchsortedlast(ν, line.nu + ν_wing_val)

    gauss_profile!(out_cache, @view(ν[ind_lo:ind_hi]), line_parameters.ν_0 + line_parameters.Δ_0, γ_D)
    for i = 1:(ind_hi - ind_lo + 1)
        data[ind_lo + i - 1] += factor * out_cache[i]
    end
end

function prepare_gauss_kwargs(;kwargs...)        
    return Dict(                
        :line_parameters => HartmannTranLineParameters(0., 0.0 * im, 0., 0., 0., 0., 0., 0.0 * im)        
    )
end

profile_map = Dict(
    :gauss => gauss_profile!,
    :lorentz => lorentz_profile!,
    :voigt => voigt_profile!,
    :speed_dependent_voigt => speed_dependent_voigt_profile, # unused for now
    :hartmann_tran => hartmann_tran_profile!,
    :rautian => rautian_profile, # unused for now
    :speed_dependent_rautian => speed_dependent_rautian_profile # unused for now
)

lineshape_map = Dict(    
    :hartmann_tran => hartmann_tran_lineshape,
    :voigt => voigt_lineshape,
    :lorentz => lorentz_lineshape,
    :gauss => gauss_lineshape   
)

profile_preparation_map = Dict(
    :hartmann_tran => prepare_hartmann_tran_kwargs,
    :voigt => prepare_voigt_kwargs,
    :lorentz => prepare_voigt_kwargs,
    :gauss => prepare_gauss_kwargs 
)

# Fadeeva function
wofz(z) = erfcx(-im * z)