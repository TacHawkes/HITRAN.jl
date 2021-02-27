const default_environments = Dict(
    # only the main isotopes are used right now
    :dry_air => Dict(
        (22, 1) =>    0.780848,       # N2
        (7, 1)  =>    0.209390,        # O2
        (2, 1)  =>    0.0004,          # CO2
        (6, 1)  =>    0.0000015,        # CH4 
        (45, 1) =>    0.0000005,       # H2
        (4, 1)  =>    0.0000003,         # N2O                
        (5, 1)  =>    0.0000002         # N2O 
    )
)

function p_s_h2o(temp)
    if (temp <= 273.15)
        return p_s_h2o_s(temp)
    else
        return p_s_h2o_l(temp)
    end
end

function p_s_h2o_l(temp)
    ϑ = temp + c_n9 / (temp - c_n10)
    A = ϑ^2 + c_n1*ϑ + c_n2
    B = c_n3*ϑ^2 + c_n4*ϑ + c_n5
    C = c_n6*ϑ^2 + c_n7*ϑ + c_n8
    return 1e6 * (2C / (-B + sqrt(B^2 - 4A*C)))^4
end

function p_s_h2o_s(temp)
    temp = kelvin_to_celsius(temp)
    return exp(c_a - c_b / (temp + c_d1)) / (temp + c_d2)^2
end

function f(p, temp)
    return c_α + c_β * p + c_γ * kelvin_to_celsius(temp)^2
end

kelvin_to_celsius(temp) = temp - 273.15
celsius_to_kelvin(temp) = temp + 273.15

"""
    moist_air(humidity [, pressure=c_p_ref, temp=c_T_ref])

Returns a component list for moist air at relative humidity with the corresponding abundances
of all components. The pressure (in atm) and the temperature (in K) have to be provided otherwise
the HITRAN defaults will be used.

!!! info "Valid range"
    Please note that the underlying model for the saturation vapor pressure uses separate
    models for water and ice. It should provide reasonable values within the range of 200 to 400 K
    and between 0.6 to 1.1 atm.
"""
function moist_air(humidity, pressure=c_p_ref, temp=c_T_ref)
    # copy the dry air environment and add the 
    # right amount of H2O
    env = copy(default_environments[:dry_air])

    # water volume fraction
    x_w = humidity / 100. * f(pressure*c_atm_factor, temp) * p_s_h2o(temp) / pressure / c_atm_factor

    for (k, v) in env
        env[k] *= (1 - x_w)
    end    
    env[(1, 1)] = x_w
    return env
end