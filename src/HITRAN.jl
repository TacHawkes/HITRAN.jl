__precompile__()

module HITRAN

using Printf, CSV, Downloads, SQLite, FileIO, JLD2,
     Interpolations, SpecialFunctions

const module_path = dirname(pathof(HITRAN))

include("constants.jl")
include("meta/components.jl")
include("meta/parameters.jl")
include("meta/environments.jl")
include("meta/tips.jl")
include("database.jl")
include("profiles.jl")
include("spectra.jl")

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
        optical_depth
end
