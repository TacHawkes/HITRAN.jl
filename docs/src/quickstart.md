# Atmospheric transmission at the oxygen A-band

The oxygen A-band is a spectral band of molecular oxygen $O_2$ in the wavelength region of about 759 nanometers to 771 nanometers. Molecular oxygen is quite remarkable because its ground state is a triplet configuration for the electrons and as everything would oxidize quite rapidly in our atmosphere if this would not be the case, we should be happy for this fact. The A-band is basically formed by the electronic transition from the ground state to the second excited singlet state in its lowest vibrational state. The band consists of a lot of individual transitions involving different rotational quantum states. For further reading pick your favorite spectroscopy textbook and have a look at the details (it is a rabbit hole with no end... and that's the fun part ;)).

Let's say we want to know the A-band spectrum at standard HITRAN conditions ($T=296\,\text{K}, p=1.0\,\text{atm}$) and what is the maximum transmission through a 100 m air column. For simplicity we assume that the conditions (and abundance does not change along this air column.

## Getting started

We start by loading the module:

```@repl
using HITRAN
```

Now we could specify a database using the function `open_database` or just stick to the default database. If no database is specified, the package will just use a default database in a file named `HITRAN.sqlite` in your current environment folder.

## Populating the local database

We need line-by-line data, so let's fetch it using:

```@setup o2_demo
using HITRAN, Plots
```

```@repl o2_demo
fetch!("StdAtm", iso_id(["N2", "O2", "CO2", "H2O", "CH4"]), 12900, 13200, [:standard, :ht_self]);
```

Not this command needs some explanation. First of all it is named `fetch!`with exclamation mark because it modifies the underlying database and is therefore considered "mutating". The parameters are the following:

| parameter | description |      example value |
| :---   |    :---:    |       ---: |
| `name`    |   A user-chosen identifier used as table name in the database  |      `"StdAtm"` |
| `global_iso_ids`   | A list of [global isotopologue IDs](https://hitran.org/docs/iso-meta/). The `iso_id` function is a convenience function which allows to lookup iso ids by formula (either molecule or isotopologue) or by molecule/local id| `iso_id(["N2", "O2", "CO2", "H2O", "CH4"])` |
| `ν_min`   | The minimum wavenumber (in $cm^{-1}$) to consider | `12900` |
| `ν_max`   | The maximum wavenumber (in $cm^{-1}$) to consider | `13200` |
| `parameters`   | A list of HITRAN line parameters to get. As nobody wants to remember these, there are shortcuts defined as symbols. The default HITRAN parameter set can be specified with :standard. In this case we additionally want to fetch the Hartmann-Tran paramters (if available). | `[:standard, :ht_self]` |

!!! note "Why are we fetching N2, CO2, H2O and CH4 instead of just O2?"
    For this example it would be sufficient to just fetch O2 spectral lines as there are only some weak H2O lines in the A-band. But it showcases how a table containing multiple components can be setup.

Now that we have a table called "StdAtm" (for Standard Atmosphere) with some of the more relevant atmospheric constituents, we can already calculate a spectrum. One more important detail: The fetch commands generates a hash of the url used to query the HITRAN database. If you make a consecutive `fetch!` call with exactly the same parameters, no new download will be initiated. This is also handy because you do not need to comment/uncomment `fetch!` commands in your code to disable download. As soon as you change the parameters (for example the wavenumber limit) new data will be downloaded.

## Calculating a spectrum

The most simple approach would be to just straight away call this command:

```@repl o2_demo
wavenumbers, absorption_coefficient = α(
    ["StdAtm"]
)
```

As you can see, we get a range of wavenumber values and values for the absorption coefficient right away. But wait, how does the module know the atmospheric composition we are interested in? In short: It does not! If you just specify a table name to the `α` function it will just calculate all spectral lines withing the table assuming their natural isotopologue abundance and does not take the proper gas mixture into account. See the function description of [`α`](@ref) for all details. For this example, we can use a neat shortcut:

```@repl o2_demo
wavenumbers, absorption_coefficient = α(
    ["StdAtm"];
    components=default_environments[:dry_air]
)
```

The `default_environment` variable is a `Dict` containing standard mixtures. For now this is only dry air by using the key `:dry_air`.
See [Environments](@ref) for details for the environments.

!!! warning Diluent parameter
    The diluent parameter behaves a little bit differently for HITRAN.jl 0.1.1 and greater and differs from HAPI.
    The problem with the HAPI diluent specification is that the diluent `:self` will apply to ALL gases in a mixture and is therefore
    wrong by design. To fix this, HITRAN.jl can actually work out the diluent itself given a gas mixture. The `self` portion will be
    set to abundance specified in components and the `air` portion is attributed to the remaining fraction. If H2O is part of the mixture
    the `H2O` diluent will be set accordingly and the influence of water on the collisional broadening will be taken into account (if the necessary data is supplied by the HITRAN database).
    
    The behaviour of the diluent parameter matches the HAPI behaviour for a single component specification. If you want to calculate
    a gas mixture, you have to provide a `Dict` with the molecule/local isotopologue id as key and another `Dict`as value containing the diluent
    as usual. It is possible to provide only diluent information for some components, the module will work out the other diluents automatically.


### The final code

To sum it all up, the correct code (besides the warning above) for this example is:

```@repl o2_demo
wavenumbers, absorption_coefficient = α(["StdAtm"];
        components=(default_environments[:dry_air])                
    )
```

This specifies that our oxygen is dilutet by the surrounding air mass which tells the module how to properly calculate the collision-broadening of the spectral lines.

To get the transmission along our air column we can use the little helper function and specifying the air path length in centimeters:

```@repl o2_demo
transmission = transmittance_spectrum(absorption_coefficient, 100e2);
```

## Full code & plotting

Let's put everything so far together and create a plot:

```@example
using HITRAN, Plots

fetch!("StdAtm", iso_id(["N2", "O2", "CO2", "H2O", "CH4"]), 12900, 13200, [:standard, :ht_self]);
wavenumbers, absorption_coefficient = α(["StdAtm"];
        components=(default_environments[:dry_air]),        
        diluent=Dict(:self => 0.209390, :air => 1 - 0.209390)        
    )
transmission = transmittance_spectrum(absorption_coefficient, 100e2)

plot(
    wavenumbers,
    transmission, 
    xlabel="Wavenumbers [1/cm]", 
    ylabel="Transmission", 
    title="Transmission along a 100 m air column",
    leg=false
)
savefig("air-plot.svg"); nothing # hide
```

![](air-plot.svg)