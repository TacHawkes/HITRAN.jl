absorption_spectrum(α::AbstractVector, len) = @. 1 - exp(-α*len)
transmittance_spectrum(α::AbstractVector, len) = @. exp(-α*len)
optical_depth(α::AbstractVector, len) = @. α*len