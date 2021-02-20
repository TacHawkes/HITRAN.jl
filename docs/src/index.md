# Introduction and Motivation

HITRAN.jl is a Julia package for calculating spectral features using the HITRAN database[^Gordon2017]. The package follows a very similar approach to the [HITRAN Application Programming Interface (HAPI)](https://github.com/hitranonline/hapi) [^Kochanov2016]. In fact the basic workflow and methods work almost the same which should make switching for HAPI users easy.
The package fetches spectral line data from HITRANOnline and stores it in a local SQLite database for structured access. This data can be used to calculate spectral lineshapes while taking into account the gas mixture and environmetal conditions (pressure, temperature).

## Motivation

This package is an "exercise" package for me in order to learn Julia while I am a main time scientific Python user. My idea was to create something useful for the Julia community. As my research history, present and future is related to measuring spectra and comparing them to HITRAN spectral simulations I thought Julia could need its own/native spectral simulation package. There is a package called [Jultran.jl](https://github.com/jsbj/Jultran.jl) but it seems abandonded and incomplete.

Please note that this package is not part of my employment but a mere hobby project.

## Features

* Native Julia implementation
* Download and manage HITRAN spectral line data in a local SQLite database
* Calculate temperature dependent line strengths using precalculated total internal partition sums from HITRAN[^Gamache2017]
* Calculate absorption coefficients using the flexible Hartmann-Tran [^Ngo2013] line shape model (also used to imlement speed-dependent Voigt/Rautian-Sobelman lineshapes, as well as regular Voigt and Rautian-Sobelman lineshapes). Simple lineshapes like Lorentz and Gauss/Doppler are also included.
* Simple and convenient syntax
* Reasonable memory footprint and performance (roughly 1 order of magnitude faster than HAPI for the shown O2 A-band example)

## Installation

Install the package using the package manager:

```julia
] add HITRAN
```

## Notable differences to HAPI and missing features

* Total internal partition sums use the [precalculated] (https://hitran.org/docs/iso-meta/) values converted to a JLD2 binary format instead of the original Python/Fortran wrapper of [^Gamache2017]
* CODATA2018[^Tiesinga2019] values are used for physical constants throughout the package and SI is used instead of cgs wherever possible
* The Faddeeva function is calculated using Julia's [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl) package which is in turn using openlibm. Surprisingly (or not !?) this is about 30% faster than the algorithms described in [^Humlíček1982] and [^Tran2013][^Tran2014]
* There are no convenience functions at the moment to filter the local database. However, direct SQL access to the database is possible which allows to create views to be used for spectral calculations. Therefore you can use the full power of SQL(ite) to sub-select/filter HITRAN data.
* No instrument functions yet
* No line-mixing
* No choice for other TIPS implementations
* SQL schema / queries are not optimized for high performance database queries. However, the performance limiting factor is the lineshape math and not the database access

## References

[^Gordon2017]: I. E. Gordon, L. S. Rothman, C. Hill, R. V. Kochanov, Y. Tan, et al., *The HITRAN2016 molecular spectroscopic database*, J. Quant. Spectrosc. Radiat. Transfer 203, 3-69 (2017)

[^Kochanov2016]: R.V. Kochanov, I.E. Gordon, L.S. Rothman, P. Wcislo, C. Hill, J.S. Wilzewski, *HITRAN Application Programming Interface (HAPI): A comprehensive approach to working with spectroscopic data*, J. Quant. Spectrosc. Radiat. Transfer 177, 15-30 (2016). [Link to article](http://dx.doi.org/10.1016/j.jqsrt.2016.03.005)

[^Gamache2017]: Robert R. Gamache, Christopher Roller, Eldon Lopes, Iouli E. Gordon, Laurence S. Rothman, Oleg L. Polyansky, Nikolai F. Zobov, Aleksandra A. Kyuberis, Jonathan Tennyson, Sergei N. Yurchenko, Attila G. Császár, Tibor Furtenbacher, Xinchuan Huang, David W. Schwenke, Timothy J. Lee, Brian J. Drouin, Sergei A. Tashkun, Valery I. Perevalov, Roman V. Kochanov, *Total internal partition sums for 166 isotopologues of 51 molecules important in planetary atmospheres: Application to HITRAN2016 and beyond*, J. Quant. Spectrosc. Radiat. Transfer 203, 70-87 (2017). [Link to article] (https://doi.org/10.1016/j.jqsrt.2017.03.045)

[^Ngo2013]: N. H. Ngo, D. Lisak, H. Tran, J.-M. Hartmann, *An isolated line-shape model to go beyond the Voigt profile in spectroscopic databases and radiative transfer codes*, J. Quant. Spectrosc. Radiat. Transfer 129, 89-100 (2013). [Link to article] (http://www.sciencedirect.com/science/article/pii/S0022407313002598)

[^Humlíček1982]: J. Humlíček, *Optimized computation of the voigt and complex probability functions*, J. Quant. Spectrosc. Radiat. Transfer 27, 437-444 (1982). [Link to article] (https://doi.org/10.1016/0022-4073(82)90078-4)

[^Tran2013]: H. Tran, N. H. Ngo, J.-M. Hartmann, *Efficient computation of some speed-dependent isolated line profiles*, J. Quant. Spectrosc. Radiat. Transfer 129, 199-203 (2013). [Link to article] (http://www.sciencedirect.com/science/article/pii/S0022407313002598)

[^Tran2014]: H. Tran, N. H. Ngo, J.-M. Hartmann, *Erratum to "Efficient computation of some speed-dependent isolated line profiles"*, J. Quant. Spectrosc. Radiat. Transfer 134, 104 (2014). [Link to article] (http://www.sciencedirect.com/science/article/pii/S0022407313004445)

[^Tiesinga2019]: Eite Tiesinga, Peter J. Mohr, David B. Newell, and Barry N. Taylor: *CODATA Recommended Values of the Fundamental Physical Constants: 2018* (2019)