
function instrument_rectangular(x, res)
    y = zeros(eltype(x), length(x))
    y[-res / 2 .<= x .< res / 2] .= 1 / res
    return y
end

function instrument_triangular(x, res)
    y = zeros(eltype(x), length(x))
    y[abs.(x) .<= res] .= 1 / res * ( 1 .- abs.(@view(x[abs.(x) .<= res])) ./ res)
    return y
end

instrument_gaussian(x, res) = @. 2*√(log(2)/π)/res * exp(-4*log(2)*(x/res)^2)

instrument_lorentzian(x, res) = @. res / π / (x^2 + res^2)

function instrument_cosine(x, res) 
    y = zeros(eltype(x), length(x))
    y[abs.(x) .<= res] = @. cos(π*@view(x[abs.(x) .<= res])/(2*res))/res*pi/4
    return y
end

instrument_diffraction(x, res) = @. 1/res * sinc(x/res)^2

instrument_michelson(x, res) = @. 2/res * sinc(2/res*x)

const instrument_functions = Dict{Symbol, Function}(
    :rectangular => instrument_rectangular,
    :triangular => instrument_triangular,
    :gaussian => instrument_gaussian,
    :lorentzian => instrument_lorentzian,
    :cosine => instrument_cosine,
    :diffraction => instrument_diffraction,
    :michelson => instrument_michelson
)

"""
    apply_instrument_function(ν, α[, instrument_function=:rectangular, instrument_wing=10.0, instrument_resolution=0.1])

Applies an instrument function to the given input spectrum. 

# Arguments
- `ν`: The wavenumber vector
- `α`: The calculated absorption coefficient using [`α`](@ref)
- `instrument_function` (optional): A Symbol describing one of the instrument functions below
- `instrument_wing` (optional): The half-width of the range for calculating the instrument function in ``cm^{-1}``
- `instrument_resolution` (optional): The full-width of the instrument resolution in ``cm^{-1}``

# Output

Returns a new vector with the spectrum influenced by the instrument function

# Instrument functions

The following instrument functions are supported. Use the stated symbol as value
for the argument `instrument_function`.

| Name |  Symbol | Description |
| :---   | :--- |      ---: |
| Rectangular | `:rectangular` |  A rectangular instrument function (e.g. a slit) |
| Triangular | `:triangular` |  A triangular instrument function |
| Gaussian | `:gaussian` |  A Gaussian instrument function (e.g. a broadband source) |
| Lorentzian | `:lorentzian` |  A Lorentzian instrument function (e.g. a single frequency laser) |
| Cosine | `:cosine` |  A cosine instrument function |
| Diffraction | `:diffraction` |  A diffraction (sinc-type) instrument function |
| Michelson | `:michelson` |  A Michelson interferometer-type instrument function (e.g. FTIR) |

"""
function apply_instrument_function(
    ν::AbstractVector{T},
    α::AbstractVector{T},
    instrument_function::Symbol=:rectangular,
    instrument_wing::AbstractFloat=10.0,
    instrument_resolution::AbstractFloat=0.1
) where T <: AbstractFloat
    δν = ν[2] - ν[1]
    x = -instrument_wing:δν:instrument_wing
    ifn = instrument_functions[instrument_function](x, instrument_resolution)    
    ifn ./= sum(ifn)
    α_inst = conv(α, ifn)
    ind_lo = length(ifn) ÷ 2 + 1
    ind_hi = length(α_inst) - ind_lo + 2    
    return α_inst[ind_lo:ind_hi]
end