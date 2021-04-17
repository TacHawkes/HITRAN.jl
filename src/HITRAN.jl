__precompile__()

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

end
