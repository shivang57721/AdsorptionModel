# Competitive dual-site Langmuir isotherm for zeolite (CO2 and N2)

# Shared helper: computes (q_CO2, q_N2) in mol/kg from concentrations in mol/m³
function _zeolite_dual_langmuir(T, c_CO2, c_N2)
    b₀_CO2  = 8.65e-7;  b₀_N2  = 2.5e-6
    d₀_CO2  = 2.63e-10; # d₀_N2 = 0
    ΔUb_CO2 = -36641.21; ΔUb_N2 = -1.58e4
    ΔUd_CO2 = -35690.66 # ΔUd_N2 = 0
    q_sb_CO2 = 3.09; q_sb_N2 = 5.84
    q_sd_CO2 = 2.54  # q_sd_N2 = 0
    R = oftype(T, 8.314)

    b_CO2 = b₀_CO2 * exp(-ΔUb_CO2 / (R * T))
    b_N2  = b₀_N2  * exp(-ΔUb_N2  / (R * T))
    d_CO2 = d₀_CO2 * exp(-ΔUd_CO2 / (R * T))

    sum_bc = b_CO2 * c_CO2 + b_N2 * c_N2

    q_CO2 = q_sb_CO2 * b_CO2 * c_CO2 / (1 + sum_bc) +
            q_sd_CO2 * d_CO2 * c_CO2 / (1 + d_CO2 * c_CO2)
    q_N2  = q_sb_N2  * b_N2  * c_N2  / (1 + sum_bc)

    return q_CO2, q_N2
end

function q_star_CO2_zeolite(gas, params)
    _zeolite_dual_langmuir(gas.T, gas.c_CO2, gas.c_N2)[1]
end

function q_star_N2_zeolite(gas, params)
    _zeolite_dual_langmuir(gas.T, gas.c_CO2, gas.c_N2)[2]
end
