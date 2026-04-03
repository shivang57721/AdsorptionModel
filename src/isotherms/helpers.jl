# Shared isotherm utility functions

Psat_H2O(T) = 611.21 * exp((18.678 - T / 234.5) * T / (T + 273.15 - 16.01))

# Smooth positivity (AD/stiff-friendly)
@inline function _softplus(x, β::T=convert(typeof(x), 500)) where {T}
    z = clamp(β*x, oftype(x, -50), oftype(x, 50))
    return z > zero(z) ? x + inv(β)*log1p(exp(-z)) : inv(β)*log1p(exp(z))
end

@inline safeexp(x) = exp(clamp(x, oftype(x,-100), oftype(x,100)))
