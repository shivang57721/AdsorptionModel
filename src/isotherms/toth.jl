# Toth isotherm models for CO2

function Toth_WADST_isotherm_CO2_wet(gas, q_H2O, params)
    q∞, b₀, ΔH₀, τ₀, α, T_ref      = params.q∞, params.b₀, params.ΔH₀, params.τ₀, params.α, params.T_ref
    b₀_wet, ΔH₀_wet, τ₀_wet, α_wet = params.b₀_wet, params.ΔH₀_wet, params.τ₀_wet, params.α_wet
    q∞_wet, A                      = params.q∞_wet, params.A
    R = oftype(gas.T, 8.314)

    b_dry = b₀     * safeexp(-ΔH₀     / (R*gas.T))
    b_wet = b₀_wet * safeexp(-ΔH₀_wet / (R*gas.T))

    τ_dry = τ₀     + α     * (1 - T_ref / gas.T)
    τ_wet = τ₀_wet + α_wet * (1 - T_ref / gas.T)

    p_CO2p = _softplus(gas.p_CO2)
    q_H2Op = _softplus(q_H2O)

    bp_dry = b_dry * p_CO2p
    bp_wet = b_wet * p_CO2p

    denom_dry = (1 + bp_dry ^ τ_dry) ^ (1 / τ_dry)
    denom_wet = (1 + bp_wet ^ τ_wet) ^ (1 / τ_wet)

    q_dry = q∞     * bp_dry / denom_dry
    q_wet = q∞_wet * bp_wet / denom_wet

    w_wet = safeexp(-A / q_H2Op)
    q_star = (1 - w_wet) * q_dry + w_wet * q_wet

    return q_star
end

function Toth_isotherm_CO2_modified_H2O_Stampi(gas, q_H2O, params)
    q∞    = params.q∞
    b₀    = params.b₀
    ΔH₀   = params.ΔH₀
    τ₀    = params.τ₀
    α     = params.α
    T_ref = params.T_ref
    γ     = params.γ
    β     = params.β

    R = 8.314

    q∞ = q∞ * (1 / (1 - γ * q_H2O))
    b  = b₀ * exp(-ΔH₀ / (R * gas.T)) * (1 + β * q_H2O)
    τ  = τ₀ + α * (1 - T_ref / gas.T)

    p_CO2_safe = max(eps(eltype(gas.p_CO2)), gas.p_CO2)

    q_star = q∞ * b * gas.p_CO2 / (1 + (b * p_CO2_safe) ^ τ) ^ (1/τ)

    return q_star
end

function Toth_isotherm_CO2_modified_H2O_Arvind(gas, q_H2O, params)
    q∞    = params.q∞
    b₀    = params.b₀
    ΔH₀   = params.ΔH₀
    τ₀    = params.τ₀
    α     = params.α
    T_ref = params.T_ref
    γ     = params.γ
    β     = params.β

    R = 8.314

    p_CO2_safe = _softplus(gas.p_CO2)

    q∞ = q∞ * (1 / (1 - γ * q_H2O))
    b  = b₀ * exp(-ΔH₀ / (R * T_ref) * (T_ref / gas.T - 1)) * _softplus(1 + β * q_H2O)
    τ  = τ₀ + α * (1 - T_ref / gas.T)

    q_star = q∞ * b * gas.p_CO2 / (1 + (b * p_CO2_safe) ^ τ) ^ (1/τ)

    return q_star
end
