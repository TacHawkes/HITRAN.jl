"""
    absorption_spectrum(α::AbstractVector, len)

Computes the absorption spectrum for the given length `len` in centimeters.
The vector `α` should be calculated using the corresponding [`α`](@ref) function.
"""
absorption_spectrum(α::AbstractVector, len) = @. 1 - exp(-α*len)

"""
    transmittance_spectrum(α::AbstractVector, len)

Computes the transmittance spectrum for the given length `len` in centimeters.
The vector `α` should be calculated using the corresponding [`α`](@ref) function.
"""
transmittance_spectrum(α::AbstractVector, len) = @. exp(-α*len)

"""
    optical_depth(α::AbstractVector, len)

Computes the optical depth for the given length `len` in centimeters.
The vector `α` should be calculated using the corresponding [`α`](@ref) function.
"""
optical_depth(α::AbstractVector, len) = @. α*len