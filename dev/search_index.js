var documenterSearchIndex = {"docs":
[{"location":"spectra/#Functions-for-absorption-spectrum-calculation","page":"Calculating spectra","title":"Functions for absorption spectrum calculation","text":"","category":"section"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"The most important function of the module is the α function. Before using i, you have to initialise a database and use the fetch! function to retrieve line-by-line data.","category":"page"},{"location":"spectra/#The-absorption-coefficient","page":"Calculating spectra","title":"The absorption coefficient","text":"","category":"section"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"    α(tables::AbstractVector{String}, profile=:hartmann_tran;kwargs...)","category":"page"},{"location":"spectra/#HITRAN.α","page":"Calculating spectra","title":"HITRAN.α","text":"α(tables::AbstractVector{String} [, profile=:hartmann_tran; kwargs...])\n\nComputes the absorption coefficient using line-by-line data stored in the database tables specified in tables. The lineshape can be optionally specified using the profile argument and one of the Symbol keys :hartmann_tran, :voigt, :sdvoigt, :lorentz, :gauss. If no keyword arguments are specified, they will be automatically chosen from the tables provided.\n\nKeyword arguments\n\ncomponents: the components of the gas mixture for the calculation. Can be either a vector of tuples with (molecule_id, local_iso_id)               or a Dict with the (molecule_id, local_iso_id) tuple as key and the abundance as value. If the vector of tuples is supplied               the natural abundance will be used, so this makes no sense for gas mixtures other than isotopologues of the same molecule.                \nintensity_threshold: the minimum line strength in cm^-1(textmolecule cdot cm^-2)\npressure: the environmental pressure in atmospheres (default: 1.0 atm)\ntemperature: the environmental temperature in Kelvin (default: 296.0 K)\nν_range: a tuple of the form (νmin, νmax) where νmin/νmax is the minimum/maximum wavenumber for the absorption_spectrum respectively in cm^-1\nν_step: the wavenumber step in cm^-1 (default: 0.01 cm^-1)\nν_wing: absolute calculation width of a line in cm^-1 (default: 0 cm^-1)\nν_wing_hw: relative calculation width of a line in multiples of a half-width (default: 50)\ndiluent: a Dict of the diluting substances, specified as Symbol, e.g. :air or :H2O for the key and the relative concentration as key (default: Dict(:self => 1.0)). See the example for details.\n\n\n\n\n\n","category":"function"},{"location":"spectra/#Environments","page":"Calculating spectra","title":"Environments","text":"","category":"section"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"There are two default environments defined at the moment, which ease the calculation of atmospheric spectra.","category":"page"},{"location":"spectra/#Dry-air","page":"Calculating spectra","title":"Dry air","text":"","category":"section"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"The first environment is dry air which can be accessed using default_environments[:dry_air]. The following composition is used [Picard2008]. Only HITRAN listed gases are used:","category":"page"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"[Picard2008]: A., Picard, R.S., Davis, M., Gläser and K., Fujii (2008), Revised formula for the density of moist air (CIPM-2007), Metrologia 45, 149–155 (2008).","category":"page"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"Gas Volume (ppmv)\nN_2 780,848\nO_2 209,390\nCO_2 400\nCH_4 1.5\nH_2 0.5\nN_2O 0.3\nCO 0.2","category":"page"},{"location":"spectra/#Moist-air","page":"Calculating spectra","title":"Moist air","text":"","category":"section"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"The second environment model is moist air which takes the relative humidity into the account. You can use the function moist_air to get a  composition dictionary with the correct water concentration.","category":"page"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"    moist_air","category":"page"},{"location":"spectra/#HITRAN.moist_air","page":"Calculating spectra","title":"HITRAN.moist_air","text":"moist_air(humidity [, pressure=c_p_ref, temp=c_T_ref])\n\nReturns a component list for moist air at relative humidity with the corresponding abundances of all components. The pressure (in atm) and the temperature (in K) have to be provided otherwise the HITRAN defaults will be used.\n\ninfo: Valid range\nPlease note that the underlying model for the saturation vapor pressure uses separate models for water and ice. It should provide reasonable values within the range of 200 to 400 K and between 0.6 to 1.1 atm.\n\n\n\n\n\n","category":"function"},{"location":"spectra/#Instrument-functions","page":"Calculating spectra","title":"Instrument functions","text":"","category":"section"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"    apply_instrument_function","category":"page"},{"location":"spectra/#HITRAN.apply_instrument_function","page":"Calculating spectra","title":"HITRAN.apply_instrument_function","text":"apply_instrument_function(ν, α[, instrument_function=:rectangular, instrument_wing=10.0, instrument_resolution=0.1])\n\nApplies an instrument function to the given input spectrum. \n\nArguments\n\nν: The wavenumber vector\nα: The calculated absorption coefficient using α\ninstrument_function (optional): A Symbol describing one of the instrument functions below\ninstrument_wing (optional): The half-width of the range for calculating the instrument function in cm^-1\ninstrument_resolution (optional): The full-width of the instrument resolution in cm^-1\n\nOutput\n\nReturns a new vector with the spectrum influenced by the instrument function\n\nInstrument functions\n\nThe following instrument functions I(x Δ) are supported. Here xis the coordinate for evaluating the function, whose range is given by instrument_wing. Δ is the resolution parameter instrument_resolution. Use the stated symbol as value for the argument instrument_function.\n\nSymbol Equation Description\n:rectangular begincases frac1Δ  lvert x rvert leq fracΔ2  0  lvert x rvert  fracΔ2 endcases A rectangular instrument function\n:triangular begincases frac1Δ (1 - fraclvert x rvertΔ)  lvert x rvert leq Δ  0  lvert x rvert  Δ endcases A triangular instrument function\n:gaussian frac2Δ sqrtfracmathrmln2pi mathrmexp left (- mathrmln2 left ( frac2xΔright)^2 right ) A Gaussian instrument function (e.g. a broadband source)\n:lorentzian fracΔ2pi frac1x^2+left(fracΔ2right)^2 A Lorentzian instrument function (e.g. a single frequency laser)\n:cosine begincases fracpi4Δ cos left ( frac pi lvert x rvert2Δ right )  lvert x rvert leq Δ  0  lvert x rvert  Δ endcases A cosine instrument function\n:diffraction frac1Δ mathrmsinc^2 left(  fracpi xΔ right) A diffraction (sinc-type) instrument function\n:michelson frac2Δ mathrmsinc left(  frac2 pi xΔ right) A Michelson interferometer-type instrument function (e.g. FTIR)\n\n\n\n\n\n","category":"function"},{"location":"spectra/#Overview-of-instrument-functions","page":"Calculating spectra","title":"Overview of instrument functions","text":"","category":"section"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"The following graph shows all supported instrument functions for a resolution of 01 cm^-1.","category":"page"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"using HITRAN, Plots\n\nx = -0.2:0.001:0.2\nplot()\nfor (s, fn) in HITRAN.instrument_functions\n    plot!(x, fn(x, 0.1), label=String(s), lw=2)\nend\n\nsavefig(\"plot.svg\")","category":"page"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"(Image: )","category":"page"},{"location":"spectra/#Convenience-functions","page":"Calculating spectra","title":"Convenience functions","text":"","category":"section"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"The following three functions are provided as a convenience to convert the absorption coefficient to either an absorption spectrum, transmittance spectrum or an optical depth.","category":"page"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"    absorption_spectrum","category":"page"},{"location":"spectra/#HITRAN.absorption_spectrum","page":"Calculating spectra","title":"HITRAN.absorption_spectrum","text":"absorption_spectrum(α, len)\n\nComputes the absorption spectrum for the given length len in centimeters. The vector α should be calculated using the corresponding α function.\n\n\n\n\n\n","category":"function"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"    transmittance_spectrum","category":"page"},{"location":"spectra/#HITRAN.transmittance_spectrum","page":"Calculating spectra","title":"HITRAN.transmittance_spectrum","text":"transmittance_spectrum(α, len)\n\nComputes the transmittance spectrum for the given length len in centimeters. The vector α should be calculated using the corresponding α function.\n\n\n\n\n\n","category":"function"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"    optical_depth","category":"page"},{"location":"spectra/#HITRAN.optical_depth","page":"Calculating spectra","title":"HITRAN.optical_depth","text":"optical_depth(α, len)\n\nComputes the optical depth for the given length len in centimeters. The vector α should be calculated using the corresponding α function.\n\n\n\n\n\n","category":"function"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"    frequency_to_wavenumber","category":"page"},{"location":"spectra/#HITRAN.frequency_to_wavenumber","page":"Calculating spectra","title":"HITRAN.frequency_to_wavenumber","text":"frequency_to_wavenumber(x)\n\nConverts a given frequency in Hz to wavenumber in cm^-1\n\n\n\n\n\n","category":"function"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"    wavelength_to_wavenumber","category":"page"},{"location":"spectra/#HITRAN.wavelength_to_wavenumber","page":"Calculating spectra","title":"HITRAN.wavelength_to_wavenumber","text":"wavelength_to_wavenumber(x)\n\nConverts a given wavelength in m to wavenumber in cm^-1\n\n\n\n\n\n","category":"function"},{"location":"quickstart/#Atmospheric-transmittance-at-the-oxygen-A-band","page":"Quickstart/Example","title":"Atmospheric transmittance at the oxygen A-band","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"The oxygen A-band is a spectral band of molecular oxygen O_2 in the wavelength region of about 759 nanometers to 771 nanometers. Molecular oxygen is quite remarkable because its ground state is a triplet configuration for the electrons and as everything would oxidize quite rapidly in our atmosphere if this would not be the case, we should be happy for this fact. The A-band is basically formed by the electronic transition from the ground state to the second excited singlet state in its lowest vibrational state. The band consists of a lot of individual transitions involving different rotational quantum states. For further reading pick your favorite spectroscopy textbook and have a look at the details (it is a rabbit hole with no end... and that's the fun part ;)).","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"Let's say we want to know the A-band spectrum at standard HITRAN conditions (T=296textK p=10textatm) and what is the maximum transmittance through a 100 m air column. For simplicity we assume that the conditions (and abundance does not change along this air column.","category":"page"},{"location":"quickstart/#Getting-started","page":"Quickstart/Example","title":"Getting started","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"We start by loading the module:","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"using HITRAN","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"Now we could specify a database using the function open_database or just stick to the default database. If no database is specified, the package will just use a default database in a file named HITRAN.sqlite in your current environment folder.","category":"page"},{"location":"quickstart/#Populating-the-local-database","page":"Quickstart/Example","title":"Populating the local database","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"We need line-by-line data, so let's fetch it using:","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"using HITRAN, Plots","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"fetch!(\"StdAtm\", iso_id([\"N2\", \"O2\", \"CO2\", \"H2O\", \"CH4\"]), 12900, 13200, [:standard, :ht_self]);","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"Now this command needs some explanation. First of all it is named fetch!with exclamation mark because it modifies the underlying database and is therefore considered \"mutating\". The parameters are the following:","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"parameter description example value\nname A user-chosen identifier used as table name in the database \"StdAtm\"\nglobal_iso_ids A list of global isotopologue IDs. The iso_id function is a convenience function which allows to lookup iso ids by formula (either molecule or isotopologue) or by molecule/local id iso_id([\"N2\", \"O2\", \"CO2\", \"H2O\", \"CH4\"])\nν_min The minimum wavenumber (in cm^-1) to consider 12900\nν_max The maximum wavenumber (in cm^-1) to consider 13200\nparameters A list of HITRAN line parameters to get. As nobody wants to remember these, there are shortcuts defined as symbols. The default HITRAN parameter set can be specified with :standard. In this case we additionally want to fetch the Hartmann-Tran paramters (if available). [:standard, :ht_self]","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"note: Why are we fetching N2, CO2, H2O and CH4 instead of just O2?\nFor this example it would be sufficient to just fetch O2 spectral lines as there are only some weak H2O lines in the A-band. But it showcases how a table containing multiple components can be setup.","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"Now that we have a table called \"StdAtm\" (for Standard Atmosphere) with some of the more relevant atmospheric constituents, we can already calculate a spectrum. One more important detail: The fetch commands generates a hash of the url used to query the HITRAN database. If you make a consecutive fetch! call with exactly the same parameters, no new download will be initiated. This is also handy because you do not need to comment/uncomment fetch! commands in your code to disable download. As soon as you change the parameters (for example the wavenumber limit) new data will be downloaded.","category":"page"},{"location":"quickstart/#Calculating-a-spectrum","page":"Quickstart/Example","title":"Calculating a spectrum","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"The most simple approach would be to just straight away call this command:","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"wavenumbers, absorption_coefficient = α(\n    [\"StdAtm\"]\n)","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"As you can see, we get a range of wavenumber values and values for the absorption coefficient right away. But wait, how does the module know the atmospheric composition we are interested in? In short: It does not! If you just specify a table name to the α function it will just calculate all spectral lines withing the table assuming their natural isotopologue abundance and does not take the proper gas mixture into account. See the function description of α for all details. For this example, we can use a neat shortcut:","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"wavenumbers, absorption_coefficient = α(\n    [\"StdAtm\"];\n    components=default_environments[:dry_air]\n)","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"The default_environment variable is a Dict containing standard mixtures. For now this is only dry air by using the key :dry_air. See Environments for details for the environments.","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"warning: Diluent parameter\nThe diluent parameter behaves a little bit differently for HITRAN.jl 0.1.1 and greater and differs from HAPI. The problem with the HAPI diluent specification is that the diluent :self will apply to ALL gases in a mixture and is therefore wrong by design. To fix this, HITRAN.jl can actually work out the diluent itself given a gas mixture. The self portion will be set to abundance specified in components and the air portion is attributed to the remaining fraction. If H2O is part of the mixture the H2O diluent will be set accordingly and the influence of water on the collisional broadening will be taken into account (if the necessary data is supplied by the HITRAN database).The behaviour of the diluent parameter matches the HAPI behaviour for a single component specification. If you want to calculate a gas mixture, you have to provide a Dict with the molecule/local isotopologue id as key and another Dictas value containing the diluent as usual. It is possible to provide only diluent information for some components, the module will work out the other diluents automatically.","category":"page"},{"location":"quickstart/#The-final-code","page":"Quickstart/Example","title":"The final code","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"To sum it all up, the correct code (besides the warning above) for this example is:","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"wavenumbers, absorption_coefficient = α([\"StdAtm\"];\n        components=(default_environments[:dry_air])                \n    )","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"This specifies that our oxygen is dilutet by the surrounding air mass which tells the module how to properly calculate the collision-broadening of the spectral lines.","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"To get the transmittance along our air column we can use the little helper function and specifying the air path length in centimeters:","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"transmittance = transmittance_spectrum(absorption_coefficient, 100e2);","category":"page"},{"location":"quickstart/#Full-code-and-plotting","page":"Quickstart/Example","title":"Full code & plotting","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"Let's put everything so far together and create a plot:","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"using HITRAN, Plots\n\nfetch!(\"StdAtm\", iso_id([\"N2\", \"O2\", \"CO2\", \"H2O\", \"CH4\"]), 12900, 13200, [:standard, :ht_self]);\nwavenumbers, absorption_coefficient = α([\"StdAtm\"];\n        components=default_environments[:dry_air]        \n    )\ntransmittance = transmittance_spectrum(absorption_coefficient, 100e2)\n\nplot(\n    wavenumbers,\n    transmittance,\n    xlabel=\"Wavenumbers [1/cm]\", \n    ylabel=\"Transmission\", \n    title=\"Transmission along a 100 m air column\",\n    leg=false\n)","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"using HITRAN, Plots\n\nfetch!(\"StdAtm\", iso_id([\"N2\", \"O2\", \"CO2\", \"H2O\", \"CH4\"]), 12900, 13200, [:standard, :ht_self]);\nwavenumbers, absorption_coefficient = α([\"StdAtm\"];\n        components=default_environments[:dry_air]        \n    )\ntransmittance = transmittance_spectrum(absorption_coefficient, 100e2)\n\nplot(\n    wavenumbers,\n    transmittance, \n    xlabel=\"Wavenumbers [1/cm]\", \n    ylabel=\"Transmission\", \n    title=\"Transmission along a 100 m air column\",\n    leg=false\n)\nsavefig(\"air-plot.svg\"); nothing # hide","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"(Image: )","category":"page"},{"location":"#Introduction-and-Motivation","page":"Introduction & Motivation","title":"Introduction and Motivation","text":"","category":"section"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"HITRAN.jl is a Julia package for calculating spectral features using the HITRAN database[Gordon2017]. The package follows a very similar approach to the HITRAN Application Programming Interface (HAPI) [Kochanov2016]. In fact the basic workflow and methods work almost the same which should make switching for HAPI users easy. The package fetches spectral line data from HITRANOnline and stores it in a local SQLite database for structured access. This data can be used to calculate spectral lineshapes while taking into account the gas mixture and environmetal conditions (pressure, temperature).","category":"page"},{"location":"#Motivation","page":"Introduction & Motivation","title":"Motivation","text":"","category":"section"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"This package is an \"exercise\" package for me in order to learn Julia while I am a main time scientific Python user. My idea was to create something useful for the Julia community. As my research history, present and future is related to measuring spectra and comparing them to HITRAN spectral simulations I thought Julia could need its own/native spectral simulation package. There is a package called Jultran.jl but it seems abandonded and incomplete.","category":"page"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"Please note that this package is not part of my employment but a mere hobby project.","category":"page"},{"location":"#Features","page":"Introduction & Motivation","title":"Features","text":"","category":"section"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"Native Julia implementation\nDownload and manage HITRAN spectral line data in a local SQLite database\nCalculate temperature dependent line strengths using precalculated total internal partition sums from HITRAN[Gamache2017]\nCalculate absorption coefficients using the flexible Hartmann-Tran [Ngo2013] line shape model (also used for the Voigt lineshape). Simple lineshapes like Lorentz and Gauss/Doppler are also included.\nSimple and convenient syntax\nReasonable memory footprint and performance (roughly 1 order of magnitude faster than HAPI for the shown O2 A-band example)\nOptimized for performance and large databases","category":"page"},{"location":"#Installation","page":"Introduction & Motivation","title":"Installation","text":"","category":"section"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"Install the package using the package manager:","category":"page"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"] add HITRAN","category":"page"},{"location":"#Citation","page":"Introduction & Motivation","title":"Citation","text":"","category":"section"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"If you find this package useful and use it for your research, please cite it as:","category":"page"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"Kliebisch, Oliver. (2021, March 31). HITRAN.jl - A julia package for calculation absorption spectral using the HITRAN database. Zenodo. http://doi.org/10.5281/zenodo.4653317","category":"page"},{"location":"#Notable-differences-to-HAPI-and-missing-features","page":"Introduction & Motivation","title":"Notable differences to HAPI and missing features","text":"","category":"section"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"Total internal partition sums use the precalculated values converted to a JLD2 binary format instead of the original Python/Fortran wrapper of [Gamache2017]\nCODATA2018[Tiesinga2019] values are used for physical constants throughout the package and SI is used instead of cgs wherever possible\nThe Faddeeva function is calculated using Julia's SpecialFunctions.jl package which is in turn using openlibm. Surprisingly (or not !?) this is about 30% faster than the algorithms described in [Humlíček1982] and [Tran2013][Tran2014]\nThere are no convenience functions at the moment to filter the local database. However, direct SQL access to the database is possible which allows to create views to be used for spectral calculations. Therefore you can use the full power of SQL(ite) to sub-select/filter HITRAN data.\nNo line-mixing\nNo choice for other TIPS implementations\nSQL schema / queries are not optimized for high performance database queries.","category":"page"},{"location":"#References","page":"Introduction & Motivation","title":"References","text":"","category":"section"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"[Gordon2017]: I. E. Gordon, L. S. Rothman, C. Hill, R. V. Kochanov, Y. Tan, et al., The HITRAN2016 molecular spectroscopic database, J. Quant. Spectrosc. Radiat. Transfer 203, 3-69 (2017)","category":"page"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"[Kochanov2016]: R.V. Kochanov, I.E. Gordon, L.S. Rothman, P. Wcislo, C. Hill, J.S. Wilzewski, HITRAN Application Programming Interface (HAPI): A comprehensive approach to working with spectroscopic data, J. Quant. Spectrosc. Radiat. Transfer 177, 15-30 (2016). Link to article","category":"page"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"[Gamache2017]: Robert R. Gamache, Christopher Roller, Eldon Lopes, Iouli E. Gordon, Laurence S. Rothman, Oleg L. Polyansky, Nikolai F. Zobov, Aleksandra A. Kyuberis, Jonathan Tennyson, Sergei N. Yurchenko, Attila G. Császár, Tibor Furtenbacher, Xinchuan Huang, David W. Schwenke, Timothy J. Lee, Brian J. Drouin, Sergei A. Tashkun, Valery I. Perevalov, Roman V. Kochanov, Total internal partition sums for 166 isotopologues of 51 molecules important in planetary atmospheres: Application to HITRAN2016 and beyond, J. Quant. Spectrosc. Radiat. Transfer 203, 70-87 (2017). Link to article","category":"page"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"[Ngo2013]: N. H. Ngo, D. Lisak, H. Tran, J.-M. Hartmann, An isolated line-shape model to go beyond the Voigt profile in spectroscopic databases and radiative transfer codes, J. Quant. Spectrosc. Radiat. Transfer 129, 89-100 (2013). Link to article","category":"page"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"[Humlíček1982]: J. Humlíček, Optimized computation of the voigt and complex probability functions, J. Quant. Spectrosc. Radiat. Transfer 27, 437-444 (1982). Link to article","category":"page"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"[Tran2013]: H. Tran, N. H. Ngo, J.-M. Hartmann, Efficient computation of some speed-dependent isolated line profiles, J. Quant. Spectrosc. Radiat. Transfer 129, 199-203 (2013). Link to article","category":"page"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"[Tran2014]: H. Tran, N. H. Ngo, J.-M. Hartmann, Erratum to \"Efficient computation of some speed-dependent isolated line profiles\", J. Quant. Spectrosc. Radiat. Transfer 134, 104 (2014). Link to article","category":"page"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"[Tiesinga2019]: Eite Tiesinga, Peter J. Mohr, David B. Newell, and Barry N. Taylor: CODATA Recommended Values of the Fundamental Physical Constants: 2018 (2019)","category":"page"},{"location":"database/#Creating-and-maintaining-the-local-HITRAN-database-cache","page":"Local HITRAN database","title":"Creating and maintaining the local HITRAN database cache","text":"","category":"section"},{"location":"database/","page":"Local HITRAN database","title":"Local HITRAN database","text":"The local HITRAN database cache is stored using the SQLite format. You can directly use and manipulate this database in order to filter for specific properties. In future versions some easy-to-use Julia filter functions will be added. For now you can create a view on a table (or multiple tables) and store them as a view which can be used by the α function as a source table.","category":"page"},{"location":"database/#Basic-database-usage","page":"Local HITRAN database","title":"Basic database usage","text":"","category":"section"},{"location":"database/","page":"Local HITRAN database","title":"Local HITRAN database","text":"    open_database(file_path::String)","category":"page"},{"location":"database/#HITRAN.open_database-Tuple{String}","page":"Local HITRAN database","title":"HITRAN.open_database","text":"open_database(file_path::String)\n\nOpens the SQLite database at the given file path and sets it as the current database.  Only use this if you want to use multiple different database.\n\n\n\n\n\n","category":"method"},{"location":"database/","page":"Local HITRAN database","title":"Local HITRAN database","text":"    fetch!","category":"page"},{"location":"database/#HITRAN.fetch!","page":"Local HITRAN database","title":"HITRAN.fetch!","text":"fetch!([db,] name, global_ids, ν_min, ν_max, parameters)\n\nFetches new data from HITRANonline and stores it in the current database in the table given by name. If the table with the given parameters already exists, no data download will be initiated.\n\nArguments\n\ndb: The database to use for storage (optional)\nname: The table name which can subsequently used as source table for the α function\nglobal_ids: The global isotopologue ids to consider. You can also provide a tuple (or an array of tuples) with molecule_id, local_id as identifiers\nν_min: The minimum wavenumber in cm^-1 to consider\nν_max: The minimum wavenumber in cm^-1 to consider\nparameters: A list of parameters to fetch. You can use parameter groups using Symbols as shortcuts, e.g. :standard for all HITRAN standard parameters. \n\n\n\n\n\n","category":"function"},{"location":"database/#Parameter-groups","page":"Local HITRAN database","title":"Parameter groups","text":"","category":"section"},{"location":"database/","page":"Local HITRAN database","title":"Local HITRAN database","text":"There is a number of parameter groups you can use to define parameter lists for fetching line-by-line data. For now it is recommended to use :standard and combine it with :ht_air. A list of all parameter groups will be listed in the future.","category":"page"},{"location":"database/#Isotopologue-IDs","page":"Local HITRAN database","title":"Isotopologue IDs","text":"","category":"section"},{"location":"database/","page":"Local HITRAN database","title":"Local HITRAN database","text":"    iso_id","category":"page"},{"location":"database/#HITRAN.iso_id","page":"Local HITRAN database","title":"HITRAN.iso_id","text":"iso_id([db::SQLite.DB,] M::T, I::T) where T <: Union{Integer, AbstractVector{Integer}}\n\nReturns the global isotopologue IDs for the given molecule ids and local ids provided either as single value or as array for both parameters.\n\nArguments\n\nM: molecule id or array of ids\nI: local isotopologue id or array of ids\n\n\n\n\n\niso_id([db::SQLite.DB,] formulas::T) where T <: Union{String, AbstractVector{String}}\n\nReturns the global isotopologue IDs for the given molecule or isotopologue formulas.\n\n\n\n\n\n","category":"function"},{"location":"database/","page":"Local HITRAN database","title":"Local HITRAN database","text":"","category":"page"}]
}
