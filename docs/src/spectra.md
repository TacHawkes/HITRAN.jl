# Functions for absorption spectrum calculcation

The most important function of the module is the [`α`](@ref) function.
Before using i, you have to initialise a database and use the [`fetch!`](@ref) function
to retrieve line-by-line data.

## The absorption coefficient

```@docs
    α(tables::AbstractVector{String}, profile=:hartmann_tran;kwargs...)
``` 

## Convenience functions

The following three functions are provided as a convenience to convert the absorption
coefficient to either an absorption spectrum, transmittance spectrum or an optical depth.

```@docs
    absorption_spectrum(α::AbstractVector, len)
``` 

```@docs
    transmittance_spectrum(α::AbstractVector, len)
``` 

```@docs
    optical_depth(α::AbstractVector, len)
``` 