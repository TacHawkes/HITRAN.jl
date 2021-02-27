# Functions for absorption spectrum calculation

The most important function of the module is the [`α`](@ref) function.
Before using i, you have to initialise a database and use the [`fetch!`](@ref) function
to retrieve line-by-line data.

## The absorption coefficient

```@docs
    α(tables::AbstractVector{String}, profile=:hartmann_tran;kwargs...)
``` 

## Environments

There are two default environments defined at the moment, which ease the calculation of atmospheric spectra.

### Dry air
The first environment is dry air which can be accessed using `default_environments[:dry_air]`. The following composition is used [^Picard2008]. Only HITRAN listed gases are used:

[^Picard2008]: A., Picard, R.S., Davis, M., Gläser and K., Fujii (2008), Revised formula for the density of moist air (CIPM-2007), Metrologia 45, 149–155 (2008).

| Gas | Volume (ppmv) |
| :---   |       ---: |
| ``N_2`` |   780,848 |
| ``O_2`` |   209,390 |
| ``CO_2`` |      400 |
| ``CH_4`` |       1.5 |
| ``H_2`` |        0.5 |
| ``N_2O`` |       0.3 |
| ``CO``  |        0.2 |


### Moist air

The second environment model is moist air which takes the relative humidity into the account.
You can use the function `moist_air` to get a 
composition dictionary with the correct water concentration.

```@docs
    moist_air
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