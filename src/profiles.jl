Sij_T(S_ij_ref, T, T_ref, Q_T, Q_T_ref, e_lower, ν_ij)::Float64 = S_ij_ref * Q_T_ref / Q_T *
                                                                  exp(-c_c2 * e_lower / T) *
                                                                  (1 -
                                                                   exp(-c_c2 * ν_ij / T)) /
                                                                  exp(-c_c2 * e_lower /
                                                                      T_ref) * (1 -
                                                                   exp(-c_c2 * ν_ij / T_ref))

γ_Doppler(T, ν_ij, M)::Float64 = ν_ij / c_c_SI *
                                 √(2c_NA_SI * c_kB_SI * T * c_log2 / (M * 1e-3))

γ_collision_0(γ0_ref, T, T_ref, p, p_ref, n_diluent)::Float64 = γ0_ref * p / p_ref *
                                                                (T_ref / T)^n_diluent

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
    components = Dict{Tuple{Int, Int}, Float64}()
    for row in result
        components[(row.molecule_id, row.local_id)] = row.abundance
    end

    return components, components
end

function get_components(comps::AbstractVector{Tuple{T, T}}) where {T <: Integer}
    # load the natural abundances
    sql = " SELECT  molecule_id, local_id, abundance
            FROM    isotopologues
            WHERE 	(molecule_id, local_id)
            IN 		(VALUES " *
          join(["(" * join("?"^length(t), ',') * ")" for t in comps], ',') * ")"
    result = query_local_db(sql, [i for t in comps for i in t])
    # replace components by new dict
    components = Dict{Tuple{Int, Int}, Float64}()
    for row in result
        components[(row.molecule_id, row.local_id)] = row.abundance
    end

    return components, components
end

function get_components(comps::Dict{Tuple{T, T}, V}) where {T <: Integer
                                                            } where {V <:
                                                                     AbstractFloat}
    # load the natural abundances
    sql = " SELECT  molecule_id, local_id, abundance
            FROM    isotopologues
            WHERE 	(molecule_id, local_id)
            IN 		(VALUES " *
          join(["(" * join("?"^length(t), ',') * ")" for t in keys(comps)], ',') * ")"
    result = query_local_db(sql, [i for t in keys(comps) for i in t])
    # replace components by new dict
    natural_abundances = Dict{Tuple{Int, Int}, Float64}()
    for row in result
        natural_abundances[(row.molecule_id, row.local_id)] = row.abundance
    end

    return comps, natural_abundances
end

function get_diluents(diluent::Dict{Symbol, T}, components) where {T <: Number}
    if length(components) > 1
        throw(ErrorException("For gas mixtures a diluent entry for every component has to be supplied"))
    end

    if (0.0 ≲ sum(values(diluent)) ≲ 1.0) == false
        throw(ErrorException("Sum of diluent fractions must not exceed 1 or lie below zero"))
    end

    # convert to different dict format
    return Dict{Tuple{Int, Int}, Dict{Symbol, Float64}}(first(keys(components)) => diluent)
end

function get_diluents(diluent::Dict{Tuple{Int, Int}, Dict{Symbol, V}},
                      components) where {V <: Number}
    if length(diluent) > length(components)
        throw(ErrorException("More diluents than components specified."))
    end

    diluents = Dict{Tuple{Int, Int}, Dict{Symbol, Float64}}()

    # water abundance
    has_water = false
    water_abundance = 0.0
    if 1 in getindex.(keys(components), 1)
        has_water = true

        for (k, v) in components
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
                diluents[mi_tuple] = Dict(:self => abundance * (1 - water_abundance),
                                          :air => (1 - abundance) * (1 - water_abundance),
                                          :H2O => water_abundance)
            else
                diluents[mi_tuple] = Dict(:self => abundance,
                                          :air => 1 - abundance)
            end
        end

        if (0.0 ≲ sum(values(diluents[mi_tuple])) ≲ 1.0) == false
            throw(ErrorException("Sum of diluent fractions must not exceed 1 or lie below zero"))
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

function parse_kwargs(tables; kwargs...)
    intensity_threshold = get(kwargs, :intensity_threshold, -Inf)
    pressure = convert(Float64, get(kwargs, :pressure, c_p_ref))
    temperature = convert(Float64, get(kwargs, :temperature, c_T_ref))

    ν = get(kwargs, :ν, nothing)
    ν_step = convert(Float64, get(kwargs, :ν_step, 0.01))
    if ν === nothing
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
        ν_range = (convert(Float64, ν_range[1]), convert(Float64, ν_range[2]))
        ν = ν_range[1]:ν_step:ν_range[2]
    else
        ν_range = (minimum(ν), maximum(ν))
    end

    ν_min, ν_max = ν_range
    ν_wing = convert(Float64, get(kwargs, :ν_wing, 0.0))
    ν_wing_hw = convert(Float64, get(kwargs, :ν_wing_hw, 50.0))

    components, natural_abundances = get_components(get(kwargs, :components, tables))
    diluent = get(kwargs, :diluent, nothing)
    if isnothing(diluent)
        diluent = get_diluents(Dict((0, 0) => Dict(:self => 1.0)), components)
    else
        diluent = get_diluents(diluent, components)
    end

    return intensity_threshold, pressure, temperature, ν, ν_range,
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
- `ν`: a vector or range of wavenumbers for which the absorption should be calculated. Useful for reusing the output ν vector for multiple calculations. If ν is supplied ν_range will be ignored.[]
- `ν_range`: a tuple of the form (ν_min, ν_max) where ν_min/ν_max is the minimum/maximum wavenumber for the absorption_spectrum respectively in ``cm^{-1}``
- `ν_step`: the wavenumber step in ``cm^{-1}`` (default: 0.01 ``cm^{-1}``)
- `ν_wing`: absolute calculation width of a line in ``cm^{-1}`` (default: 0 ``cm^{-1}``)
- `ν_wing_hw`: relative calculation width of a line in multiples of a half-width (default: 50)
- `diluent`: a `Dict` of the diluting substances, specified as `Symbol`, e.g. `:air` or `:H2O` for the key and the relative concentration as key (default: `Dict(:self => 1.0)`). See the example for details.
"""
function α(tables::AbstractVector{String}, profile = :hartmann_tran; kwargs...)
    # parse and prepare input arguments
    intensity_threshold, pressure, temperature, ν, ν_range, ν_min, ν_max,
    ν_step, ν_wing, ν_wing_hw, diluent, components, natural_abundances = parse_kwargs(tables;
                                                                                      kwargs...)

    molar_masses = molar_mass(keys(components))

    α(tables,
      profile,
      components,
      diluent,
      intensity_threshold,
      pressure,
      temperature,
      ν,
      ν_range,
      ν_wing,
      ν_wing_hw,
      natural_abundances,
      molar_masses)
end

function α(tables,
           profile,
           components,
           diluent,
           intensity_threshold,
           pressure,
           temperature,
           ν,
           ν_range,
           ν_wing,
           ν_wing_hw,
           natural_abundances,
           molar_masses)
    # allocate output
    data = zeros(eltype(ν), length(ν))

    # lineshape function
    lineshape::Function = lineshape_map[profile]

    # molecule concentration in molecules / cm^3
    # allows conversion from cm^2/molecule -> cm^-1
    volume_conc = pressure * c_atm_factor / c_kB_SI / temperature * 1e-6

    # iterate over all tables and calculate line profiles line-by-line
    for table in tables
        # get line data (note: as :standard is enforced when fetching HITRAN
        # source data, it is assumed that at least all these fields exist and
        # no further check is performed)
        result = query_local_db("SELECT     *
                                FROM        " * table * "
                                WHERE 	    (molec_id, local_iso_id)
                                IN 		    (VALUES " *
                                join(["(" * join("?"^length(t), ',') * ")"
                                      for t in keys(components)], ',') * ")
AND         nu >= ?
AND         nu <= ?", vcat([i for t in keys(components) for i in t], [ν_range...]))

        # call profile specific generic preparation function
        if haskey(profile_preparation_map, profile)
            profile_args = profile_preparation_map[profile](; temperature, diluent,
                                                            query = result)
        else
            profile_args = ()
        end

        # field indices
        ind_sw = lookup_symbol(result, :sw)
        ind_nu = lookup_symbol(result, :nu)
        ind_elower = lookup_symbol(result, :elower)

        ind_gid = lookup_symbol(result, :global_iso_id)
        ind_mid = lookup_symbol(result, :molec_id)
        ind_lid = lookup_symbol(result, :local_iso_id)

        last_id = 0
        q_t = q_t_ref = 0.0
        for line::SQLite.Row in result
            gid = getvalue(result, ind_gid, Int)
            mid = getvalue(result, ind_mid, Int)
            lid = getvalue(result, ind_lid, Int)
            MI = (mid, lid)
            # partition sum
            if last_id != gid
                q_t = tips(gid, temperature)
                q_t_ref = tips(gid, c_T_ref)
            end

            # environment adjusted line intensity
            S_ij = Sij_T(getvalue(result, ind_sw, Float64),
                         temperature,
                         c_T_ref,
                         q_t,
                         q_t_ref,
                         getvalue(result, ind_elower, Float64),
                         getvalue(result, ind_nu, Float64))

            # skip weak lines
            if S_ij < intensity_threshold
                continue
            end

            factor = volume_conc * S_ij * components[MI] /
                     natural_abundances[MI]

            # call lineshape function
            lineshape(result,
                      diluent[MI],
                      temperature,
                      pressure,
                      molar_masses[MI],
                      ν,
                      ν_wing,
                      ν_wing_hw,
                      factor,
                      data,
                      profile_args)
        end
    end

    return ν, data
end

α(tables::String, profile = :hartmann_tran; kwargs...) = α([tables], profile; kwargs...)

function hartmann_tran_reference_temperature(T)::Float64
    for (T_range, T_ref) in c_HT_T_ref
        if T >= T_range[1] && T < T_range[2]
            return T_ref
        end
    end
end

function hartmann_tran_profile!(out,
                                ν,
                                ν_0,
                                ν_VC,
                                γ_D,
                                γ_0,
                                γ_2,
                                Δ_0,
                                Δ_2,
                                η,
                                factor = 1.0)
    C_0 = γ_0 + im * Δ_0
    C_2 = γ_2 + im * Δ_2
    C_0t = (1.0 - η) * (C_0 - 3.0 * C_2 / 2.0) + ν_VC
    C_2t = (1.0 - η) * C_2
    ν_a0 = c_c_SI * γ_D / (√(c_log2) * ν_0)

    c_ν = √π * c_c_SI / ν_0

    if (abs(C_2t) ≈ 0.0)
        @batch per=core for i in 1:length(ν)
            Z_m = (im * (ν_0 - ν[i]) + C_0t) / (ν_0 * ν_a0 / c_c_SI)
            wofz_Z_m = wofz(im * Z_m)
            A_ν = c_ν / ν_a0 * wofz_Z_m
            B_ν = c_ν * ν_a0 * ((1.0 - Z_m^2.0) * wofz_Z_m + Z_m / √(π))
            @inbounds out[i] += factor / π * real(A_ν /
                                     (1.0 - (ν_VC - η * (C_0 - 3.0 * C_2 / 2.0) * A_ν +
                                       (η * C_2 / ν_a0^2.0) * B_ν)))
        end
    else
        Y = (ν_0 * ν_a0 / 2.0 / c_c_SI / C_2t)^2
        @batch per=core for i in 1:length(ν)
            X = (im * (ν_0 - ν[i]) + C_0t) / C_2t
            Z_p = √(X + Y) + √(Y)
            Z_m = √(X + Y) - √(Y)

            wofz_Z_m = wofz(im * Z_m)
            wofz_Z_p = wofz(im * Z_p)

            A_ν = √(π) * c_c_SI / ν_0 / ν_a0 * (wofz_Z_m - wofz_Z_p)
            B_ν = ν_a0^2 / C_2t^2 * (-1.0 + √(π) / 2.0 / √(Y) * (1.0 - Z_m^2.0) * wofz_Z_m -
                   √(π) / 2.0 / √(Y) * (1.0 - Z_p^2.0) * wofz_Z_p)
            @inbounds out[i] += factor / π * real(A_ν /
                                     (1.0 - (ν_VC - η * (C_0 - 3.0 * C_2 / 2.0) * A_ν +
                                       (η * C_2 / ν_a0^2.0) * B_ν)))
        end
    end
end

# specialized getvalue functions which resemble the SQLite.getvalue functions but optimized for the type
# of access we have in this module
function getvalue(q::SQLite.Query, col::Int, ::Type{T}) where {T}
    handle = SQLite._get_stmt_handle(q.stmt)
    t = SQLite.C.sqlite3_column_type(handle, col)
    if t == SQLite.C.SQLITE_NULL
        return zero(T)
    else
        return SQLite.sqlitevalue(T, handle, col)
    end
end

function getvalue(q::SQLite.Query, col::Int, col_fb::Int, ::Type{T}) where {T}
    handle = SQLite._get_stmt_handle(q.stmt)
    t = SQLite.C.sqlite3_column_type(handle, col)
    t_fb = SQLite.C.sqlite3_column_type(handle, col_fb)
    if t == SQLite.C.SQLITE_NULL && t_fb == SQLite.C.SQLITE_NULL
        return zero(T)
    elseif t != SQLite.C.SQLITE_NULL
        return SQLite.sqlitevalue(T, handle, col)
    else
        return SQLite.sqlitevalue(T, handle, col_fb)
    end
end

function get_line_parameter(q::SQLite.Query, index::Int,
                            ::Type{T} = Float64) where {T <: AbstractFloat}
    getvalue(q, index, T)
end
function get_line_parameter(q::SQLite.Query, index::Int, fallback_index::Int,
                            ::Type{T} = Float64) where {T <: AbstractFloat}
    getvalue(q, index, fallback_index, T)
end
function get_line_parameter(q::SQLite.Query, ::Missing, index::Int,
                            ::Type{T} = Float64) where {T <: AbstractFloat}
    get_line_parameter(q, index, T)
end
function get_line_parameter(q::SQLite.Query, index::Int, ::Missing,
                            ::Type{T} = Float64) where {T <: AbstractFloat}
    get_line_parameter(q, index, T)
end
function get_line_parameter(::SQLite.Query, ::Missing, ::Missing,
                            ::Type{T} = Float64) where {T <: AbstractFloat}
    zero(T)
end
function get_line_parameter(::SQLite.Query, ::Missing,
                            ::Type{T} = Float64) where {T <: AbstractFloat}
    zero(T)
end

function hartmann_tran_lineshape(q,
                                 diluent,
                                 temperature,
                                 pressure,
                                 mass,
                                 ν,
                                 ν_wing,
                                 ν_wing_hw,
                                 factor,
                                 data,
                                 ht_args)
    ν_0 = get_line_parameter(q, q.lookup[:nu])
    ind_air::Int = q.lookup[:n_air]
    γ_D = γ_Doppler(temperature, ν_0, mass)
    γ_0 = γ_2 = Δ_0 = Δ_2 = zero(Float64)
    ν_VC = η = zero(ComplexF64)

    # loop over all diluents and build combined line parameters
    γ_0_dil = n_dil = T_ref = γ_0t = Δ_0_dil = Δ_0p_dil = Δ_0t = γ_2_dil = γ_2t = zero(Float64)
    Δ_2_dil = Δ_2t = η_dil = ν_VC_dil = κ_dil = zero(Float64)
    for (diluent_name, diluent_abundance) in diluent
        # get Hartmann-Tran or Voigt parameters if available
        ht_dil = ht_args.fields[diluent_name]
        vg_dil = ht_args.voigt_args.fields[diluent_name]
        # γ_0 contribution
        γ_0_dil = get_line_parameter(q, ht_dil.γ_0, vg_dil.γ_0)
        n_dil = get_line_parameter(q, ht_dil.n, vg_dil.n)
        if (diluent_name == :self || n_dil == c_default_zero)
            n_dil = get_line_parameter(q, ind_air)
        end
        T_ref = (ht_dil.γ_0 !== missing && ht_dil.n !== missing) ? ht_args.T_ht : c_T_ref
        γ_0t = γ_collision_0(γ_0_dil, temperature, T_ref, pressure,
                             c_p_ref, n_dil)
        γ_0 += diluent_abundance * γ_0t

        # Δ_0 contribution
        Δ_0_dil = get_line_parameter(q, ht_dil.Δ_0, vg_dil.Δ_0)
        Δ_0p_dil = get_line_parameter(q, ht_dil.Δ_0p, vg_dil.Δ_0p)
        T_ref = (ht_dil.Δ_0 !== missing && ht_dil.Δ_0p !== missing) ? ht_args.T_ht : c_T_ref
        Δ_0t = Δ_0_dil + Δ_0p_dil * (temperature - T_ref) * pressure / c_p_ref
        Δ_0 += diluent_abundance * Δ_0t

        # γ_2 contribution
        γ_2_dil = get_line_parameter(q, ht_dil.γ_2)
        γ_2t = γ_2_dil * pressure / c_p_ref
        γ_2 += diluent_abundance * γ_2t

        #  Δ_2 contribution
        Δ_2_dil = get_line_parameter(q, ht_dil.Δ_2)
        Δ_2t = Δ_2_dil * pressure / c_p_ref
        Δ_2 += diluent_abundance * Δ_2t

        # η contribution
        η_dil = get_line_parameter(q, ht_dil.η)
        η += diluent_abundance * η_dil * (γ_2t - im * Δ_2t)

        # ν_VC contribution
        ν_VC_dil = get_line_parameter(q, ht_dil.ν_VC)
        κ_dil = get_line_parameter(q, ht_dil.κ)
        ν_VC += diluent_abundance * ν_VC_dil * (T_ref / temperature)^κ_dil * pressure
        ν_VC -= η_dil * diluent_abundance * (γ_0t - im * Δ_0t)
    end
    if η != 0.0
        η /= γ_2 - im * Δ_2
    end

    ν_VC += η * (γ_0 - im * Δ_0)
    # use absolute or hw wing specification?
    ν_wing_val = max(ν_wing, ν_wing_hw * γ_0, ν_wing_hw * γ_D)

    # find a suitable bisection for integrating the data into the vector
    ind_lo = searchsortedfirst(ν, ν_0 - ν_wing_val)
    ind_hi = searchsortedlast(ν, ν_0 + ν_wing_val)

    hartmann_tran_profile!(@view(data[ind_lo:ind_hi]), @view(ν[ind_lo:ind_hi]), ν_0, ν_VC,
                           γ_D, γ_0, γ_2, Δ_0, Δ_2, η, factor)
end

function lookup_symbol(q::SQLite.Query, s1::Symbol)
    s1 in keys(q.lookup) && !isa(Any, q.types[q.lookup[s1]]) ? q.lookup[s1] : missing
end

abstract type ProfileArguments end

struct VoigtDiluentIndices
    γ_0::Union{Missing, Int}
    n::Union{Missing, Int}
    Δ_0::Union{Missing, Int}
    Δ_0p::Union{Missing, Int}
end

struct VoigtProfileArguments <: ProfileArguments
    fields::Dict{Symbol, VoigtDiluentIndices}
end

struct HartmannTranDiluentIndices
    γ_0::Union{Missing, Int}
    n::Union{Missing, Int}
    Δ_0::Union{Missing, Int}
    Δ_0p::Union{Missing, Int}
    γ_2::Union{Missing, Int}
    Δ_2::Union{Missing, Int}
    η::Union{Missing, Int}
    ν_VC::Union{Missing, Int}
    κ::Union{Missing, Int}
end

struct HartmannTranProfileArguments <: ProfileArguments
    T_ht::Float64
    fields::Dict{Symbol, HartmannTranDiluentIndices}
    voigt_args::VoigtProfileArguments
end

function prepare_hartmann_tran_kwargs(; kwargs...)
    T = get(kwargs, :temperature, c_T_ref)
    diluent = get(kwargs, :diluent) do
        Dict()
    end
    T_ht = hartmann_tran_reference_temperature(T)
    q = get(kwargs, :query, nothing)

    # prepare field names for diluents
    fields = Dict{Symbol, HartmannTranDiluentIndices}()
    #Dict{Symbol,Union{Int,Missing}}
    # get all diluent names
    for diluent_name in get_diluent_names(diluent)
        # lookup fields

        fields[diluent_name] = HartmannTranDiluentIndices(lookup_symbol(q,
                                                                        Symbol("gamma_HT_0_$(diluent_name)_$(T_ht)")),
                                                          lookup_symbol(q,
                                                                        Symbol("n_HT_0_$(diluent_name)_$(T_ht)")),
                                                          lookup_symbol(q,
                                                                        Symbol("delta_HT_0_$(diluent_name)_$(T_ht)")),
                                                          lookup_symbol(q,
                                                                        Symbol("deltap_HT_0_$(diluent_name)_$(T_ht)")),
                                                          lookup_symbol(q,
                                                                        Symbol("gamma_HT_2_$(diluent_name)_$(T_ht)")),
                                                          lookup_symbol(q,
                                                                        Symbol("delta_HT_2_$(diluent_name)_$(T_ht)")),
                                                          lookup_symbol(q,
                                                                        Symbol("eta_HT_$(diluent_name)")),
                                                          lookup_symbol(q,
                                                                        Symbol("nu_HT_$(diluent_name)")),
                                                          lookup_symbol(q,
                                                                        Symbol("kappa_HT_$(diluent_name)")))
    end
    voigt_args = prepare_voigt_kwargs(; kwargs...)
    return HartmannTranProfileArguments(T_ht,
                                        fields,
                                        voigt_args)
end

function voigt_profile!(out,
                        ν,
                        ν_0,
                        γ_D,
                        γ_0,
                        factor = 1.0)
    hartmann_tran_profile!(out, ν, ν_0, 0.0 + im * 0.0, γ_D, γ_0, 0.0, 0.0, 0.0,
                           0.0 + im * 0.0, factor)
end

function voigt_lineshape(q,
                         diluent,
                         temperature,
                         pressure,
                         mass,
                         ν,
                         ν_wing,
                         ν_wing_hw,
                         factor,
                         data,
                         voigt_args)
    ν_0::Float64 = get_line_parameter(q, q.lookup[:nu])
    ind_air::Int = q.lookup[:n_air]
    # initialize lineshape specific parameters
    γ_D = γ_Doppler(temperature, ν_0, mass)
    γ_0 = Δ_0 = Float64(0.0)

    # loop over all diluents and build combined line parameters
    for (diluent_name, diluent_abundance) in diluent
        # get Voigt parameters if available
        vg_dil = voigt_args.fields[diluent_name]

        # γ_0 contribution
        γ_0_dil = get_line_parameter(q, vg_dil.γ_0)
        n_dil = get_line_parameter(q, vg_dil.n)
        if (diluent_name == :self || n_dil == c_default_zero)
            n_dil = get_line_parameter(q, ind_air)
        end
        γ_0t = γ_collision_0(γ_0_dil, temperature, c_T_ref, pressure,
                             c_p_ref, n_dil)
        γ_0 += diluent_abundance * γ_0t

        # Δ_0 contribution
        Δ_0_dil = get_line_parameter(q, vg_dil.Δ_0)
        Δ_0p_dil = get_line_parameter(q, vg_dil.Δ_0p)
        Δ_0t = Δ_0_dil + Δ_0p_dil * (temperature - c_T_ref) * pressure / c_p_ref
        Δ_0 += diluent_abundance * Δ_0t
    end

    # use absolute or hw wing specification?
    ν_wing_val = max(ν_wing, ν_wing_hw * γ_0, ν_wing_hw * γ_D)

    # find a suitable bisection for integrating the data into the vector
    ind_lo = searchsortedfirst(ν, ν_0 - ν_wing_val)
    ind_hi = searchsortedlast(ν, ν_0 + ν_wing_val)

    voigt_profile!(@view(data[ind_lo:ind_hi]), @view(ν[ind_lo:ind_hi]), ν_0 + Δ_0, γ_D, γ_0,
                   factor)
end

function prepare_voigt_kwargs(; kwargs...)
    diluent = get(kwargs, :diluent) do
        Dict()
    end
    q = get(kwargs, :query, nothing)

    # prepare field names for diluents
    fields = Dict{Symbol, VoigtDiluentIndices}()
    for diluent_name in get_diluent_names(diluent)
        fields[diluent_name] = VoigtDiluentIndices(lookup_symbol(q,
                                                                 Symbol("gamma_$(diluent_name)")),
                                                   lookup_symbol(q,
                                                                 Symbol("n_$(diluent_name)")),
                                                   lookup_symbol(q,
                                                                 Symbol("delta_$(diluent_name)")),
                                                   lookup_symbol(q,
                                                                 Symbol("deltap_$(diluent_name)")))
    end

    return VoigtProfileArguments(fields)
end

function speed_dependent_voigt_profile!(out,
                                        ν,
                                        ν_0,
                                        γ_D,
                                        γ_0,
                                        γ_2,
                                        Δ_0,
                                        factor = 1.0)
    hartmann_tran_profile!(out, ν, ν_0, 0.0 + 0.0 * im, γ_D, γ_0, γ_2, Δ_0, 0.0,
                           0.0 + 0.0 * im, factor)
end

function speed_dependent_voigt_lineshape(q,
                                         diluent,
                                         temperature,
                                         pressure,
                                         mass,
                                         ν,
                                         ν_wing,
                                         ν_wing_hw,
                                         factor,
                                         data,
                                         sd_args)
    ν_0::Float64 = get_line_parameter(q, q.lookup[:nu])
    ind_air::Int = q.lookup[:n_air]
    # initialize lineshape specific parameters
    γ_D = γ_Doppler(temperature, ν_0, mass)
    γ_0 = Δ_0 = γ_2 = Float64(0.0)

    # loop over all diluents and build combined line parameters
    for (diluent_name, diluent_abundance) in diluent
        # get Voigt parameters if available
        sd_dil = sd_args.fields[diluent_name]

        # γ_0 contribution
        γ_0_dil = get_line_parameter(q, sd_dil.γ_0)
        n_dil = get_line_parameter(q, sd_dil.n)
        if (diluent_name == :self || n_dil == c_default_zero)
            n_dil = get_line_parameter(q, ind_air)
        end
        γ_0t = γ_collision_0(γ_0_dil, temperature, c_T_ref, pressure,
                             c_p_ref, n_dil)
        γ_0 += diluent_abundance * γ_0t

        # Δ_0 contribution
        Δ_0_dil = get_line_parameter(q, sd_dil.Δ_0)
        Δ_0p_dil = get_line_parameter(q, sd_dil.Δ_0p)
        Δ_0t = Δ_0_dil + Δ_0p_dil * (temperature - c_T_ref) * pressure / c_p_ref
        Δ_0 += diluent_abundance * Δ_0t

        # γ_2 contribution
        γ_2_dil = get_line_parameter(q, sd_dil.γ_2)
        γ_2t = γ_2_dil * pressure / c_p_ref
        γ_2 += diluent_abundance * γ_2t
    end

    # use absolute or hw wing specification?
    ν_wing_val = max(ν_wing, ν_wing_hw * γ_0, ν_wing_hw * γ_D)

    # find a suitable bisection for integrating the data into the vector
    ind_lo = searchsortedfirst(ν, ν_0 - ν_wing_val)
    ind_hi = searchsortedlast(ν, ν_0 + ν_wing_val)

    speed_dependent_voigt_profile!(@view(data[ind_lo:ind_hi]), @view(ν[ind_lo:ind_hi]), ν_0,
                                   γ_D, γ_0, γ_2, Δ_0, factor)
end

struct SpeedDependentVoigtDiluentIndices
    γ_0::Union{Missing, Int}
    n::Union{Missing, Int}
    Δ_0::Union{Missing, Int}
    Δ_0p::Union{Missing, Int}
    γ_2::Union{Missing, Int}
end

struct SpeedDependentVoigtProfileArguments <: ProfileArguments
    fields::Dict{Symbol, SpeedDependentVoigtDiluentIndices}
end

function prepare_speed_dependent_voigt_kwargs(; kwargs...)
    diluent = get(kwargs, :diluent) do
        Dict()
    end
    q = get(kwargs, :query, nothing)

    # prepare field names for diluents
    fields = Dict{Symbol, SpeedDependentVoigtDiluentIndices}()
    for diluent_name in get_diluent_names(diluent)
        fields[diluent_name] = SpeedDependentVoigtDiluentIndices(lookup_symbol(q,
                                                                               Symbol("gamma_$(diluent_name)")),
                                                                 lookup_symbol(q,
                                                                               Symbol("n_$(diluent_name)")),
                                                                 lookup_symbol(q,
                                                                               Symbol("delta_$(diluent_name)")),
                                                                 lookup_symbol(q,
                                                                               Symbol("deltap_$(diluent_name)")),
                                                                 lookup_symbol(q,
                                                                               Symbol("sd_$(diluent_name)")))
    end

    return SpeedDependentVoigtProfileArguments(fields)
end

function lorentz_profile!(out,
                          ν,
                          ν_0,
                          γ_0,
                          factor = 1.0)
    @inbounds @simd for i in 1:length(ν)
        out[i] += factor * γ_0 / (π * (γ_0^2 + (ν[i] - ν_0)^2))
    end
end

function lorentz_lineshape(q,
                           diluent,
                           temperature,
                           pressure,
                           mass,
                           ν,
                           ν_wing,
                           ν_wing_hw,
                           factor,
                           data,
                           lorentz_args)
    ν_0::Float64 = get_line_parameter(q, q.lookup[:nu])
    ind_air::Int = q.lookup[:n_air]
    # initialize lineshape specific parameters
    γ_0 = Δ_0 = Float64(0.0)

    # loop over all diluents and build combined line parameters
    for (diluent_name, diluent_abundance) in diluent
        lor_dil = lorentz_args.fields[diluent_name]

        # γ_0 contribution
        γ_0_dil = get_line_parameter(q, lor_dil.γ_0)
        n_dil = get_line_parameter(q, lor_dil.n)
        if (diluent_name == :self || n_dil == c_default_zero)
            n_dil = get_line_parameter(q, ind_air)
        end
        γ_0t = γ_collision_0(γ_0_dil, temperature, c_T_ref, pressure,
                             c_p_ref, n_dil)
        γ_0 += diluent_abundance * γ_0t

        # Δ_0 contribution
        Δ_0_dil = get_line_parameter(q, lor_dil.Δ_0)
        Δ_0p_dil = get_line_parameter(q, lor_dil.Δ_0p)
        Δ_0t = Δ_0_dil + Δ_0p_dil * (temperature - c_T_ref) * pressure / c_p_ref
        Δ_0 += diluent_abundance * Δ_0t
    end

    # use absolute or hw wing specification?
    ν_wing_val::Float64 = max(ν_wing, ν_wing_hw * γ_0)

    # find a suitable bisection for integrating the data into the vector
    ind_lo = searchsortedfirst(ν, ν_0 - ν_wing_val)
    ind_hi = searchsortedlast(ν, ν_0 + ν_wing_val)

    lorentz_profile!(@view(data[ind_lo:ind_hi]), @view(ν[ind_lo:ind_hi]), ν_0 + Δ_0, γ_0,
                     factor)
end

function gauss_profile!(out,
                        ν,
                        ν_0,
                        γ_D,
                        factor = 1.0)
    @inbounds @simd for i in 1:length(ν)
        out[i] += factor * √(c_log2 / π) / γ_D * exp(-c_log2 * ((ν[i] - ν_0) / γ_D)^2)
    end
end

function gauss_lineshape(q,
                         diluent,
                         temperature,
                         pressure,
                         mass,
                         ν,
                         ν_wing,
                         ν_wing_hw,
                         factor,
                         data,
                         ::Any)
    # initialize lineshape specific parameters
    ν_0::Float64 = get_line_parameter(q, q.lookup[:nu])
    γ_D = γ_Doppler(temperature, ν_0, mass)
    Δ_0 = get_line_parameter(q, q.lookup[:delta_air]) * pressure / c_p_ref

    # use absolute or hw wing specification?
    ν_wing_val = max(ν_wing, ν_wing_hw * γ_D)

    ind_lo = searchsortedfirst(ν, ν_0 - ν_wing_val)
    ind_hi = searchsortedlast(ν, ν_0 + ν_wing_val)

    gauss_profile!(@view(data[ind_lo:ind_hi]), @view(ν[ind_lo:ind_hi]), ν_0 + Δ_0, γ_D,
                   factor)
end

const profile_map = (gauss = gauss_profile!,
                     lorentz = lorentz_profile!,
                     voigt = voigt_profile!,
                     sdvoigt = speed_dependent_voigt_profile!,
                     hartmann_tran = hartmann_tran_profile!)

const lineshape_map = (hartmann_tran = hartmann_tran_lineshape,
                       voigt = voigt_lineshape,
                       sdvoigt = speed_dependent_voigt_lineshape,
                       lorentz = lorentz_lineshape,
                       gauss = gauss_lineshape)

const profile_preparation_map = (hartmann_tran = prepare_hartmann_tran_kwargs,
                                 voigt = prepare_voigt_kwargs,
                                 sdvoigt = prepare_speed_dependent_voigt_kwargs,
                                 lorentz = prepare_voigt_kwargs)

# Faddeeva function
wofz(z) = erfcx(-im * z)
