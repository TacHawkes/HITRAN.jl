# HITRAN constants

# reference temperature [K]
const c_T_ref = 296.

# reference pressure [atm]
const c_p_ref = 1.

# conversion factor atm to Pa
const c_atm_factor = 101325

#= 
Hartmann-Tran reference temperatures for line parameters

Reference: 
    Wcisło P, Gordon IE, Tran H, Tan Y, Hu S-M, Campargue A, et al.
    The im- plementation of non-Voigt line profiles in the HITRAN 
    database: H2 case study. J Quant Spectrosc Radiat Transf 2016;177:75–91.
    http://dx.doi.org/10.1016/j.jqsrt.2016.01.024. 
=#
const c_HT_T_ref = Base.ImmutableDict(
    (0., 100.) => 50.,
    (100., 200.) => 150.,
    (200., 400.) => 296.,
    (400., Inf) => 700.)

## SI constants according to CODATA2018

# speed of light in m/s (SI)
const c_c_SI = 299792458

# Avogrado constant (dimensionless) (SI)
const c_NA_SI = 6.02214076e23

# Boltzmann constant in J/K (SI)
const c_kB_SI = 1.380649e-23

# Planck-constant in J*s (SI)
const c_h_SI = 6.62607015e-34

# c2 = hc/k (cgs units)
const c_c2 = 1e2*c_h_SI * c_c_SI / c_kB_SI

### profile constants
const c_default_zero = 0.0