# HITRAN constants

# reference temperature [K]
const c_T_ref = 296.0

# reference pressure [atm]
const c_p_ref = 1.0

# conversion factor atm to Pa
const c_atm_factor = 101325.0

#=
Hartmann-Tran reference temperatures for line parameters

Reference:
    Wcisło P, Gordon IE, Tran H, Tan Y, Hu S-M, Campargue A, et al.
    The im- plementation of non-Voigt line profiles in the HITRAN
    database: H2 case study. J Quant Spectrosc Radiat Transf 2016;177:75–91.
    http://dx.doi.org/10.1016/j.jqsrt.2016.01.024.
=#
const c_HT_T_ref = Base.ImmutableDict((0.0, 100.0) => 50.0,
                                      (100.0, 200.0) => 150.0,
                                      (200.0, 400.0) => 296.0,
                                      (400.0, Inf) => 700.0)

## SI constants according to CODATA2018

# speed of light in m/s (SI)
const c_c_SI = 299792458.0

# Avogrado constant (dimensionless) (SI)
const c_NA_SI = 6.02214076e23

# Boltzmann constant in J/K (SI)
const c_kB_SI = 1.380649e-23

# Planck-constant in J*s (SI)
const c_h_SI = 6.62607015e-34

# c2 = hc/k (cgs units)
const c_c2 = 1e2 * c_h_SI * c_c_SI / c_kB_SI

# natural log of 2
const c_log2 = log(2.0)

### profile constants
const c_default_zero = zero(Float64)

#=
Constants for saturation water vapor pressure taken from:

IAPWS R7-97(2012)
Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam
August 2007
=#
const c_n1 = 0.11670521452767e4
const c_n2 = -0.72421316703206e6
const c_n3 = -0.17073846940092e2
const c_n4 = 0.12020824702470e5
const c_n5 = -0.32325550322333e7
const c_n6 = 0.14915108613530e2
const c_n7 = -0.48232657361591e4
const c_n8 = 0.40511340542057e6
const c_n9 = -0.23855557567849
const c_n10 = 0.65017534844798e3

#=

Constants for ice saturation vaper pressure taken from:

Huang, Jianhua. "A Simple Accurate Formula for Calculating
Saturation Vapor Pressure of Water and Ice", Journal of Applied Meteorology and
Climatology 57, 6 (2018): 1265-1272

=#

const c_a = 43.494
const c_b = 6545.8
const c_d1 = 278.0
const c_d2 = 868.0

const c_α = 1.00062
const c_β = 3.14e-8
const c_γ = 5.60e-7
