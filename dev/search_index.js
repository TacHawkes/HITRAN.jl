var documenterSearchIndex = {"docs":
[{"location":"spectra/#Functions-for-absorption-spectrum-calculcation","page":"Calculating spectra","title":"Functions for absorption spectrum calculcation","text":"","category":"section"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"    α(tables::AbstractVector{String}, profile=:hartmann_tran;kwargs...)","category":"page"},{"location":"spectra/#HITRAN.α","page":"Calculating spectra","title":"HITRAN.α","text":"α(tables::AbstractVector{String} [, profile=:hartmann_tran; kwargs...])\n\nComputes the absorption coefficient using line-by-line data stored in the database tables specified in tables. The lineshape can be optionally specified using the profile argument and one of the Symbol keys :hartmann_tran, :voigt, :lorentz, :gauss. If no keyword arguments are specified, they will be automatically chosen from the tables provided.\n\nKeyword arguments\n\ncomponents: the components of the gas mixture for the calculation. Can be either a vector of tuples with (molecule_id, local_iso_id)               or a Dict with the (molecule_id, local_iso_id) tuple as key and the abundance as value. If the vector of tuples is supplied               the natural abundance will be used, so this makes no sense for gas mixtures other than isotopologues of the same molecule.                \nintensity_threshold: the minimum line strength in cm^-1(textmolecule cdot cm^-2)\npressure: the environmental pressure in atmospheres (default: 1.0 atm)\ntemperature: the environmental temperature in Kelvin (default: 296.0 K)\nν_range: a tuple of the form (νmin, νmax) where νmin/νmax is the minimum/maximum wavenumber for the absorption_spectrum respectively in `cm^{-1}\nν_step: the wavenumber step in cm^-1 (default: 0.01 cm^-1)\nν_wing: absolute calculation width of a line in cm^-1 (default: 0 cm^-1)\nν_wing_hw: relative calculation width of a line in multiples of a half-width (default: 50)\ndiluent: a Dict of the diluting substances, specified as Symbol, e.g. :air or :H2O for the key and the relative concentration as key (default: Dict(:self => 1.0))\n\n\n\n\n\n","category":"function"},{"location":"spectra/","page":"Calculating spectra","title":"Calculating spectra","text":"","category":"page"},{"location":"quickstart/#Atmospheric-transmission-at-the-oxygen-A-band","page":"Quickstart/Example","title":"Atmospheric transmission at the oxygen A-band","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"The oxygen A-band is a spectral band of molecular oxygen O_2 in the wavelength region of about 759 nanometers to 771 nanometers. Molecular oxygen is quite remarkable because its ground state is a triplet configuration for the electrons and as everything would oxidize quite rapidly in our atmosphere if this would not be the case, we should be happy for this fact. The A-band is basically formed by the electronic transition from the ground state to the second excited singlet state in its lowest vibrational state. The band consists of a lot of individual transitions involving different rotational quantum states. For further reading pick your favorite spectroscopy textbook and have a look at the details (it is a rabbit hole with no end... and that's the fun part ;)).","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"Let's say we want to know the A-band spectrum at standard HITRAN conditions (T=296textK p=10textatm) and what is the maximum transmission through a 100 m air column. For simplicity we assume that the conditions (and abundance does not change along this air column.","category":"page"},{"location":"quickstart/#Getting-started","page":"Quickstart/Example","title":"Getting started","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"We start by loading the module:","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"using HITRAN","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"Now we could specify a database using the function open_database or just stick to the default database. If no database is specified, the package will just use a default database in a file named HITRAN.sqlite in your current environment folder.","category":"page"},{"location":"quickstart/#Populating-the-local-database","page":"Quickstart/Example","title":"Populating the local database","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"We need line-by-line data, so let's fetch it using:","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"using HITRAN, Plots","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"fetch!(\"StdAtm\", iso_id([\"N2\", \"O2\", \"CO2\", \"H2O\", \"CH4\"]), 12900, 13200, [:standard, :ht_self]);","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"Not this command needs some explanation. First of all it is named fetch!with exclamation mark because it modifies the underlying database and is therefore considered \"mutating\". The parameters are the following:","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"parameter description example value\nname A user-chosen identifier used as table name in the database \"StdAtm\"\nglobal_iso_ids A list of global isotopologue IDs. The iso_id function is a convenience function which allows to lookup iso ids by formula (either molecule or isotopologue) or by molecule/local id iso_id([\"N2\", \"O2\", \"CO2\", \"H2O\", \"CH4\"])\nν_min The minimum wavenumber (in cm^-1) to consider 12900\nν_max The maximum wavenumber (in cm^-1) to consider 13200\nparameters A list of HITRAN line parameters to get. As nobody wants to remember these, there are shortcuts defined as symbols. The default HITRAN parameter set can be specified with :standard. In this case we additionally want to fetch the Hartmann-Tran paramters (if available). [:standard, :ht_self]","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"note: Why are we fetching N2, CO2, H2O and CH4 instead of just O2?\nFor this example it would be sufficient to just fetch O2 spectral lines as there are only some weak H2O lines in the A-band. But it showcases how a table containing multiple components can be setup.","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"Now that we have a table called \"StdAtm\" (for Standard Atmosphere) with some of the more relevant atmospheric constituents, we can already calculate a spectrum. One more important detail: The fetch commands generates a hash of the url used to query the HITRAN database. If you make a consecutive fetch! call with exactly the same parameters, no new download will be initiated. This is also handy because you do not need to comment/uncomment fetch! commands in your code to disable download. As soon as you change the parameters (for example the wavenumber limit) new data will be downloaded.","category":"page"},{"location":"quickstart/#Calculating-a-spectrum","page":"Quickstart/Example","title":"Calculating a spectrum","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"The most simple approach would be to just straight away call this command:","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"wavenumbers, absorption_coefficient = α(\n    [\"StdAtm\"]\n)","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"As you can see, we get a range of wavenumber values and values for the absorption coefficient right away. But wait, how does the module know the atmospheric composition we are interested in? In short: It does not! If you just specify a table name to the α function it will just calculate all spectral lines withing the table assuming their natural isotopologue abundance and does not take the proper gas mixture into account. See the function description of α for all details. For this example, we can use a neat shortcut:","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"wavenumbers, absorption_coefficient = α(\n    [\"StdAtm\"];\n    components=default_environments[:dry_air]\n)","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"The default_environment variable is a Dict containing standard mixtures. For now this is only dry air by using the key :dry_air. Great, what is left? Actually to match the behaviour of the HAPI, we also have to specify how our species is diluted by the environment.","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"warning: Calculate diluted mixtures separately and sum them afterwards\nThe following diluent specification is actually wrong for the H2O lines overlapping the A-band spectral region. Because they are to weak, we get away with this here. If you are in a spectral band where you have to properly adjust for the self-/air-broadening, you should calculate the spectra for each component one by one and add them later. Specifying a :self dilution applies to all specified components.","category":"page"},{"location":"quickstart/#The-final-code","page":"Quickstart/Example","title":"The final code","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"To sum it all up, the correct code (besides the warning above) for this example is:","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"wavenumbers, absorption_coefficient = α([\"StdAtm\"];\n        components=(default_environments[:dry_air]),        \n        diluent=Dict(:self => 0.209390, :air => 1 - 0.209390)        \n    )","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"This specifies that our oxygen is dilutet by the surrounding air mass which tells the module how to properly calculate the collision-broadening of the spectral lines.","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"To get the transmission along our air column we can use the little helper function and specifying the air path length in centimeters:","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"transmission = transmittance_spectrum(absorption_coefficient, 100e2);","category":"page"},{"location":"quickstart/#Full-code-and-plotting","page":"Quickstart/Example","title":"Full code & plotting","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"Let's put everything so far together and create a plot:","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"using HITRAN, Plots\n\nfetch!(\"StdAtm\", iso_id([\"N2\", \"O2\", \"CO2\", \"H2O\", \"CH4\"]), 12900, 13200, [:standard, :ht_self]);\nwavenumbers, absorption_coefficient = α([\"StdAtm\"];\n        components=(default_environments[:dry_air]),        \n        diluent=Dict(:self => 0.209390, :air => 1 - 0.209390)        \n    )\ntransmission = transmittance_spectrum(absorption_coefficient, 100e2)\n\nplot(\n    wavenumbers,\n    transmission, \n    xlabel=\"Wavenumbers [1/cm]\", \n    ylabel=\"Transmission\", \n    title=\"Transmission along a 100 m air column\",\n    leg=false\n)\nsavefig(\"air-plot.svg\"); nothing # hide","category":"page"},{"location":"quickstart/","page":"Quickstart/Example","title":"Quickstart/Example","text":"(Image: )","category":"page"},{"location":"#Introduction-and-Motivation","page":"Introduction & Motivation","title":"Introduction and Motivation","text":"","category":"section"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"HITRAN.jl is a Julia package for calculating spectral features using the HITRAN database[Gordon2017]. The package follows a very similar approach to the HITRAN Application Programming Interface (HAPI) [Kochanov2016]. In fact the basic workflow and methods work almost the same which should make switching for HAPI users easy. The package fetches spectral line data from HITRANOnline and stores it in a local SQLite database for structured access. This data can be used to calculate spectral lineshapes while taking into account the gas mixture and environmetal conditions (pressure, temperature).","category":"page"},{"location":"#Motivation","page":"Introduction & Motivation","title":"Motivation","text":"","category":"section"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"This package is an \"exercise\" package for me in order to learn Julia while I am a main time scientific Python user. My idea was to create something useful for the Julia community. As my research history, present and future is related to measuring spectra and comparing them to HITRAN spectral simulations I thought Julia could need its own/native spectral simulation package. There is a package called Jultran.jl but it seems abandonded and incomplete.","category":"page"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"Please note that this package is not part of my employment but a mere hobby project.","category":"page"},{"location":"#Features","page":"Introduction & Motivation","title":"Features","text":"","category":"section"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"Native Julia implementation\nDownload and manage HITRAN spectral line data in a local SQLite database\nCalculate temperature dependent line strengths using precalculated total internal partition sums from HITRAN[Gamache2017]\nCalculate absorption coefficients using the flexible Hartmann-Tran [Ngo2013] line shape model (also used to imlement speed-dependent Voigt/Rautian-Sobelman lineshapes, as well as regular Voigt and Rautian-Sobelman lineshapes). Simple lineshapes like Lorentz and Gauss/Doppler are also included.\nSimple and convenient syntax\nReasonable memory footprint and performance","category":"page"},{"location":"#Installation","page":"Introduction & Motivation","title":"Installation","text":"","category":"section"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"Install the package using the package manager:","category":"page"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"] add HITRAN","category":"page"},{"location":"#Notable-differences-to-HAPI-and-missing-features","page":"Introduction & Motivation","title":"Notable differences to HAPI and missing features","text":"","category":"section"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"Total internal partition sums use the precalculated values converted to a JLD2 binary format instead of the original Python/Fortran wrapper of [Gamache2017]\nCODATA2018[Tiesinga2019] values are used for physical constants throughout the package and SI is used over cgs wherever possible\nThe Faddeeva function is calculated using Julia's SpecialFunctions.jl package which is in turn using openlibm. This will maybe change in a future version to use the (presumably!?) better performing algorithm described in [Humlíček1982] and [Tran2013][Tran2014]\nThere are no convenience functions at the moment to filter the local database. However, direct SQL access to the database is possible which allows to create views to be used for spectral calculations. Therefore you can use the full power of SQL(ite) to sub-select/filter HITRAN data.\nNo instrument functions yet\nNo line-mixing\nNo choice for other TIPS implementations\nSQL schema / queries are not optimized for high performance database queries. However, the performance limiting factor is the lineshape math and not the database access","category":"page"},{"location":"#References","page":"Introduction & Motivation","title":"References","text":"","category":"section"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"[Gordon2017]: I. E. Gordon, L. S. Rothman, C. Hill, R. V. Kochanov, Y. Tan, et al., The HITRAN2016 molecular spectroscopic database, J. Quant. Spectrosc. Radiat. Transfer 203, 3-69 (2017)","category":"page"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"[Kochanov2016]: R.V. Kochanov, I.E. Gordon, L.S. Rothman, P. Wcislo, C. Hill, J.S. Wilzewski, HITRAN Application Programming Interface (HAPI): A comprehensive approach to working with spectroscopic data, J. Quant. Spectrosc. Radiat. Transfer 177, 15-30 (2016). Link to article","category":"page"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"[Gamache2017]: Robert R. Gamache, Christopher Roller, Eldon Lopes, Iouli E. Gordon, Laurence S. Rothman, Oleg L. Polyansky, Nikolai F. Zobov, Aleksandra A. Kyuberis, Jonathan Tennyson, Sergei N. Yurchenko, Attila G. Császár, Tibor Furtenbacher, Xinchuan Huang, David W. Schwenke, Timothy J. Lee, Brian J. Drouin, Sergei A. Tashkun, Valery I. Perevalov, Roman V. Kochanov, Total internal partition sums for 166 isotopologues of 51 molecules important in planetary atmospheres: Application to HITRAN2016 and beyond, J. Quant. Spectrosc. Radiat. Transfer 203, 70-87 (2017). Link to article","category":"page"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"[Ngo2013]: N. H. Ngo, D. Lisak, H. Tran, J.-M. Hartmann, An isolated line-shape model to go beyond the Voigt profile in spectroscopic databases and radiative transfer codes, J. Quant. Spectrosc. Radiat. Transfer 129, 89-100 (2013). Link to article","category":"page"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"[Humlíček1982]: J. Humlíček, Optimized computation of the voigt and complex probability functions, J. Quant. Spectrosc. Radiat. Transfer 27, 437-444 (1982). Link to article","category":"page"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"[Tran2013]: H. Tran, N. H. Ngo, J.-M. Hartmann, Efficient computation of some speed-dependent isolated line profiles, J. Quant. Spectrosc. Radiat. Transfer 129, 199-203 (2013). Link to article","category":"page"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"[Tran2014]: H. Tran, N. H. Ngo, J.-M. Hartmann, Erratum to \"Efficient computation of some speed-dependent isolated line profiles\", J. Quant. Spectrosc. Radiat. Transfer 134, 104 (2014). Link to article","category":"page"},{"location":"","page":"Introduction & Motivation","title":"Introduction & Motivation","text":"[Tiesinga2019]: Eite Tiesinga, Peter J. Mohr, David B. Newell, and Barry N. Taylor: CODATA Recommended Values of the Fundamental Physical Constants: 2018 (2019)","category":"page"}]
}
