module HITRAN

using Printf
using CSV
using Downloads
using DSP
using FileIO
using JLD2
using SHA
using SpecialFunctions
using SQLite

const module_path = dirname(pathof(HITRAN))

include("constants.jl")
include("util.jl")
include("meta/components.jl")
include("meta/parameters.jl")
include("meta/environments.jl")
include("meta/tips.jl")
include("database.jl")
include("profiles.jl")
include("spectra.jl")
include("instruments.jl")

export  current_db,
        open_database,
        fetch!,
        iso_id,
        download_HITRAN,
        query_local_db,
        tips,
        α,
        default_environments,
        moist_air,
        absorption_spectrum,
        transmittance_spectrum,
        optical_depth,
        wavelength_to_wavenumber,
        frequency_to_wavenumber,
        apply_instrument_function        

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing    
    Base.precompile(Tuple{typeof(iso_id),Vector{String}})
    Base.precompile(Tuple{Core.kwftype(typeof(α)),NamedTuple{(:components, :temperature, :ν_range, :ν_step), Tuple{Dict{Tuple{Int64, Int64}, Float64}, Float64, Tuple{Int64, Int64}, Float64}},typeof(α),Vector{String}})
    Base.precompile(Tuple{Core.kwftype(typeof(prepare_hartmann_tran_kwargs)),NamedTuple{(:temperature, :diluent, :query), Tuple{Float64, Dict{Tuple{Int64, Int64}, Dict{Symbol, Float64}}, SQLite.Query}},typeof(prepare_hartmann_tran_kwargs)})
    Base.precompile(Tuple{typeof(instrument_lorentzian),StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}},Float64})
    Base.precompile(Tuple{typeof(hartmann_tran_lineshape),SQLite.Query,Dict{Symbol, Float64},Float64,Float64,Float64,StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}},Float64,Float64,Float64,Vector{Float64},HartmannTranProfileArguments})
    Base.precompile(Tuple{typeof(merge_groups),Symbol,Vararg{Symbol, N} where N})
    Base.precompile(Tuple{typeof(fetch!),SQLite.DB,String,Vector{Int64},Int64,Int64,Vector{String}})
    Base.precompile(Tuple{typeof(moist_air),Int64,Float64,Float64})
    Base.precompile(Tuple{typeof(α),Vector{String},Symbol,Dict{Tuple{Int64, Int64}, Float64},Dict{Tuple{Int64, Int64}, Dict{Symbol, Float64}},Float64,Float64,Float64,StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}},Tuple{Float64, Float64},Float64,Float64,Dict{Tuple{Int64, Int64}, Float64},Dict{Tuple{Int64, Int64}, Float64}})
    Base.precompile(Tuple{typeof(fetch!),String,Vector{Int64},Int64,Int64,Vector{Symbol}})
    Base.precompile(Tuple{typeof(apply_instrument_function),StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}},Vector{Float64},Symbol,Float64,Float64})
end        

_precompile_()

end
