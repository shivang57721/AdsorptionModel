# GAB isotherm models for H2O

function GAB_isotherm_H2O(gas, params)
    Cm    = params.Cm
    Cg    = params.Cg
    K_ads = params.K_ads

    p_H2O_safe = _softplus(gas.p_H2O)
    x = p_H2O_safe / Psat_H2O(gas.T - 273)
    x = clamp(x, 0, 1)

    q_H2O = Cm * Cg * K_ads * x / ((1 - K_ads * x) * (1 + (Cg - 1) * K_ads * x))

    return q_H2O
end

function GAB_isotherm_H2O_Tfunction_Resins(gas, params)
    qm = params.qₘ
    C  = params.C
    D  = params.D
    F  = params.F
    G  = params.G
    R  = 8.314

    E1   = C - exp(D * gas.T)
    E2_9 = F + G * gas.T
    E10  = -44.38 * gas.T + 57220
    c = exp((E1  - E10) / R / gas.T)
    k = exp((E2_9 - E10) / R / gas.T)

    p_H2O_safe = _softplus(gas.p_H2O)
    x = p_H2O_safe / Psat_H2O(gas.T - 273)
    x = clamp(x, 0, 1)

    q_star = qm * k * c * x / ((1 - k * x) * (1 + (c - 1) * k * x))

    q_star
end
