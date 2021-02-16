# HITRAN.jl

[![Build Status](https://ci.appveyor.com/api/projects/status/github/TacHawkes/HITRAN.jl?svg=true)](https://ci.appveyor.com/project/TacHawkes/HITRAN-jl)
[![Coverage](https://codecov.io/gh/TacHawkes/HITRAN.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/TacHawkes/HITRAN.jl)

HITRAN.jl is a Julia package for calculating spectral features using the HITRAN database. The package follows a very similar approach to the [HITRAN Application Programming Interface (HAPI)](https://github.com/hitranonline/hapi). In fact the basic workflow and methods work almost the same which should make switching for HAPI users easy.
The package fetches spectral line data from HITRANOnline and stores it in a local SQLite database for structured access. This data can be used to calculate spectral lineshapes while taking into account the gas mixture and environmetal conditions (pressure, temperature).

## Documentation

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://tachawkes.github.io/HITRAN.jl/dev)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://tachawkes.github.io/HITRAN.jl/stable)

## Installation

Install the package using the package manager:

```julia
] add HITRAN
```

## References

* I. E. Gordon, L. S. Rothman, C. Hill, R. V. Kochanov, Y. Tan, et al., *The HITRAN2016 molecular spectroscopic database*, J. Quant. Spectrosc. Radiat. Transfer 203, 3-69 (2017)

* R.V. Kochanov, I.E. Gordon, L.S. Rothman, P. Wcislo, C. Hill, J.S. Wilzewski, *HITRAN Application Programming Interface (HAPI): A comprehensive approach to working with spectroscopic data*, J. Quant. Spectrosc. Radiat. Transfer 177, 15-30 (2016). [Link to article](http://dx.doi.org/10.1016/j.jqsrt.2016.03.005)