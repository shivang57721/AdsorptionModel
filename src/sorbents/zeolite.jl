# Zeolite 13X sorbent parameters

SORB_PARAMS_ZEOLITE(T=Float64) = SorbentParams(T;
                        ε_bed = 0.37,
                        ε_total = 0.37,
                        dₚ = 2.0e-3,
                        ρ_bed = (1 - 0.37) * 1130,
                        ρ_particle = 1130,
                        k_CO2 = (c_CO2, q_star) -> c_CO2 / q_star * 15 * 0.35 * (1.6e-5 / 3) / 1e-6,
                        k_H2O = 0.0,
                        k_N2  = (c_N2, q_star) -> c_N2 / q_star * 15 * 0.35 * (1.6e-5 / 3) / 1e-6,
                        ΔH_CO2 = −36641,
                        ΔH_H2O = 0.0,
                        ΔH_N2  = -15800,
                        Dₘ = 1.6e-5,
                        C_solid = 1070.0,
                        isotherm_params = ISOTHERM_PARAMS_LEWATIT(T),
                        q_star_H2O = (gas, params) -> zero(gas.p_H2O),
                        q_star_CO2 = q_star_CO2_zeolite,
                        q_star_N2  = q_star_N2_zeolite
                    )