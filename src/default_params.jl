# Default parameter factory functions

# --------------- PHYSICAL PARAMS ----------------------
PHYS_PARAMS_DEFAULT(T=Float64) = PhysicalParams(T;
                                Cₚ_CO2   = 37.520,
                                Cₚ_H2O   = 34.2,
                                λ_ev     = 2350.0,
                                Cₚ_steam = 2.08
                            )

# ----------------- COST PARAMS ------------------------
COST_PARAMS_DEFAULT(T=Float64) = CostParams(T;
                        efficiency_blower = 0.5,
                        adiabatic_index = 1.4,
                        eta_VP = 0.7,
                        annual_operating_hours = 8000,
                        target_CO2_ton_per_year = 1000)

# ----------------- COLUMN PARAMS ----------------------
COL_PARAMS_ARVIND(T=Float64) = ColumnParams(T;
                        Rᵢ = 0.04,
                        Rₒ = 0.041,
                        L  = 0.01,
                        h_L = 3.0,
                        h_W = 26.0
                    )

COL_PARAMS_YOUNG(T=Float64) = ColumnParams(T;
                        Rᵢ = 0.04,
                        Rₒ = 0.041,
                        L  = 0.01,
                        h_L = 14.0,
                        h_W = 22000
                    )

COL_PARAMS_LARGESCALE(T=Float64) = ColumnParams(T;
                        Rᵢ = 0.62,
                        Rₒ = 0.64,
                        L  = 0.02,
                        h_L = 236,
                        h_W = 22000
                    )

COL_PARAMS_HAGHPANAH(T=Float64) = ColumnParams(T;
                            Rᵢ = 0.1445,
                            Rₒ = 0.1620,
                            L  = 1,
                            h_L = 8.6,
                            h_W = 2.5
                        )

# --------------- ISOTHERM PARAMS ----------------------
ISOTHERM_PARAMS_LEWATIT(T=Float64) = IsothermParams(T;
                            # DRY SITE PARAMS
                            q∞     = 4.86,
                            b₀     = 2.85e-21, # Pa-1
                            ΔH₀    = -117798,
                            τ₀     = 0.209,
                            α      = 0.523,
                            T_ref  = 298.15,

                            # WET SITE PARAMS
                            b₀_wet  = 1.230e-18, # Pa-1
                            ΔH₀_wet = -203687,
                            τ₀_wet  = 0.053,
                            α_wet   = 0.053,
                            q∞_wet  = 9.035,
                            A       = 1.532,

                            # GAB PARAMS
                            qₘ    = 3.63,
                            C     = 47110,
                            D     = 0.023744,
                            F     = 57706,
                            G     = -47.814,

                            # Stampi model
                            γ = - 0.137,
                            β = 5.612
                        )
ISOTHERM_PARAMS_NbOFFIVE(T=Float64) = IsothermParams(T;
                            # DRY SITE PARAMS
                            q∞     = 2.22,
                            b₀     = 0.17567, # Pa-1
                            ΔH₀    = -50000,
                            τ₀     = 1.166,
                            α      = -0.4937,
                            T_ref  = 273,

                            # GAB PARAMS
                            Cg    = 36.81,
                            K_ads = 0.4063,
                            Cm    = 6.159,

                            # GAB PARAMS 2
                            qₘ    = 3.63,
                            C     = 47110,
                            D     = 0.023744,
                            F     = 57706,
                            G     = -47.814,

                            # Stampi model
                            γ = 0.0,
                            β = -0.11707
                        )

# --------------- SORBENT PARAMS ----------------------
SORB_PARAMS_LEWATIT(T=Float64; activation=1) = SorbentParams(T;
                        ε_bed = 0.4,
                        ε_total = 0.54,
                        dₚ = 0.00052,
                        ρ_bed = 528.0,
                        ρ_particle = 880.0,
                        k_CO2 = 0.003,
                        k_H2O = 0.0086,
                        k_N2 = 0.0,
                        ΔH_CO2 = -70000.0,
                        ΔH_H2O = -46000.0,
                        ΔH_N2 = 0.0,
                        Dₘ = 1.3e-5,
                        C_solid = 1580.0,
                        isotherm_params = ISOTHERM_PARAMS_LEWATIT(T),
                        q_star_H2O = GAB_isotherm_H2O_Tfunction_Resins,
                        q_star_CO2 = (T, p_CO2, q_H2O, params) -> activation * Toth_WADST_isotherm_CO2_wet(T, p_CO2, q_H2O, params)
                    )

SORB_PARAMS_NbOFFIVE(T=Float64; activation=1) = SorbentParams(T;
                        ε_bed = 0.4,
                        ε_total = 0.64,
                        dₚ = 0.0075,
                        ρ_bed = 704.0,
                        ρ_particle = 1173.6,
                        k_CO2 = 2.0e-4,
                        k_H2O = 2.0e-1,
                        k_N2 = 0.0,
                        ΔH_CO2 = -50000.0,
                        ΔH_H2O = -45036.0,
                        ΔH_N2 = 0.0,
                        Dₘ = 2.05e-5,
                        C_solid = 1000.0,
                        isotherm_params = ISOTHERM_PARAMS_NbOFFIVE(T),
                        q_star_H2O = GAB_isotherm_H2O,
                        q_star_CO2 = (T, p_CO2, q_H2O, params) -> activation * Toth_isotherm_CO2_modified_H2O_Arvind(T, p_CO2, q_H2O, params)
                    )

# --------------- PREDEFINED ISOTHERMS ----------------------

Psat_H2O(T) = 611.21 * exp((18.678 - T / 234.5) * T / (T + 273.15 - 16.01))

# Smooth positivity (AD/stiff-friendly)
@inline function _softplus(x, β::T=convert(typeof(x), 500)) where {T}
    z = clamp(β*x, oftype(x, -50), oftype(x, 50))
    return z > zero(z) ? x + inv(β)*log1p(exp(-z)) : inv(β)*log1p(exp(z))
end

@inline safeexp(x) = exp(clamp(x, oftype(x,-700), oftype(x,700)))

function Toth_WADST_isotherm_CO2_wet(T, p_CO2, q_H2O, params)
    # params
    q∞, b₀, ΔH₀, τ₀, α, T_ref      = params.q∞, params.b₀, params.ΔH₀, params.τ₀, params.α, params.T_ref
    b₀_wet, ΔH₀_wet, τ₀_wet, α_wet = params.b₀_wet, params.ΔH₀_wet, params.τ₀_wet, params.α_wet
    q∞_wet, A                      = params.q∞_wet, params.A
    R = oftype(T, 8.314)

    b_dry = b₀     * safeexp(-ΔH₀     / (R*T))
    b_wet = b₀_wet * safeexp(-ΔH₀_wet / (R*T))

    τ_dry = τ₀     + α     * (1 - T_ref / T)
    τ_wet = τ₀_wet + α_wet * (1 - T_ref / T)

    # Smooth-positive inputs
    p_CO2p = _softplus(p_CO2)
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

function Toth_isotherm_CO2_modified_H2O_Stampi(T, p_CO2, q_H2O, params)
    q∞     = params.q∞
    b₀     = params.b₀ # Pa^-1
    ΔH₀    = params.ΔH₀
    τ₀     = params.τ₀
    α      = params.α
    T_ref  = params.T_ref
    γ = params.γ
    β = params.β

    R = 8.314

    q∞ = q∞ * (1 / (1 - γ * q_H2O))
    b = b₀ * exp(-ΔH₀ / (R * T)) * (1 + β * q_H2O)
    τ = τ₀ + α * (1 - T_ref / T)

    p_CO2_safe = max(eps(eltype(p_CO2)), p_CO2)

    q_star = q∞ * b * p_CO2 / (1 + (b * p_CO2_safe) ^ τ) ^ (1/τ)

    return q_star
end

function Toth_isotherm_CO2_modified_H2O_Arvind(T, p_CO2, q_H2O, params)
    q∞     = params.q∞
    b₀     = params.b₀ # Pa^-1
    ΔH₀    = params.ΔH₀
    τ₀     = params.τ₀
    α      = params.α
    T_ref  = params.T_ref
    γ = params.γ
    β = params.β

    R = 8.314

    p_CO2_safe = _softplus(p_CO2)

    q∞ = q∞ * (1 / (1 - γ * q_H2O))
    b = b₀ * exp(-ΔH₀ / (R * T_ref) * (T_ref / T - 1)) * _softplus(1 + β * q_H2O)
    τ = τ₀ + α * (1 - T_ref / T)

    q_star = q∞ * b * p_CO2 / (1 + (b * p_CO2_safe) ^ τ) ^ (1/τ)

    return q_star
end

function GAB_isotherm_H2O(T, p_H2O, params)
    Cm = params.Cm
    Cg = params.Cg
    K_ads = params.K_ads

    p_H2O_safe = _softplus(p_H2O)
    x = p_H2O_safe / Psat_H2O(T-273)
    x = clamp(x, 0, 1)

    q_H2O = Cm * Cg * K_ads * x /((1 - K_ads * x) * (1 + (Cg - 1) * K_ads * x))

    return q_H2O
end

function GAB_isotherm_H2O_Tfunction_Resins(T, p_H2O, params)
    qm = params.qₘ
    C = params.C
    D = params.D
    F = params.F
    G = params.G
    R = 8.314

    E1 = C - exp(D * T)
    E2_9 = F + G * T
    E10 = -44.38 * T + 57220
    c = exp((E1-E10)/ R / T)
    k = exp((E2_9 - E10)/ R / T)

    p_H2O_safe = _softplus(p_H2O)
    x = p_H2O_safe / Psat_H2O(T - 273)
    x = clamp(x, 0, 1)

    q_star = qm * k * c * x / ((1 - k * x) * (1 + (c - 1) * k * x))

    q_star
end

function q_star_zeolite(u, data)
    b₀_CO2 = 8.65e-7
    b₀_N2  = 2.5e-6
    d₀_CO2 = 2.63e-10
    d₀_N2  = 0.0
    ΔUb_CO2 = -36641.21
    ΔUb_N2  = -1.58e4
    ΔUd_CO2 = -35690.66
    ΔUd_N2  = 0.0
    q_sb_CO2 = 3.09
    q_sb_N2 = 5.84
    q_sd_CO2 = 2.54
    q_sd_N2  = 0.0
    R = 8.314

    T = u[data.iT]
    c_CO2 = u[data.iCO2]
    c_N2 = u[data.iN2]

    b_CO2 = b₀_CO2 * exp(-ΔUb_CO2 / (R * T))
    b_N2 = b₀_N2 * exp(-ΔUb_N2 / (R * T))
    d_CO2 = d₀_CO2 * exp(-ΔUd_CO2 / (R * T))
    d_N2 = d₀_N2 * exp(-ΔUd_N2 / (R * T))

    sum_bc = b_CO2 * c_CO2 + b_N2 * c_N2
    sum_dc = d_CO2 * c_CO2 + d_N2 * c_N2
    q_star_CO2 = q_sb_CO2 * b_CO2 * c_CO2 / (1 + sum_bc) + q_sd_CO2 * d_CO2 * c_CO2 / (1 + sum_dc)
    q_star_N2 = q_sb_N2 * b_N2 * c_N2 / (1 + sum_bc) + q_sd_N2 * d_N2 * c_N2 / (1 + sum_dc)

    return (; q_star_H2O=0, q_star_CO2, q_star_N2)
end