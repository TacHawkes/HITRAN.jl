# some trivial helpers

"""
    frequency_to_wavenumber(x)

Converts a given frequency in Hz to wavenumber in ``cm^{-1}``
"""
frequency_to_wavenumber(x::Real) = x / (1e2*c_c_SI)

"""
    wavelength_to_wavenumber(x)

Converts a given wavelength in m to wavenumber in ``cm^{-1}``
"""
wavelength_to_wavenumber(x::Real) = 1e-2 / x

# float tolerant leq method
≲(a, b) = (a <= b) || a ≈ b