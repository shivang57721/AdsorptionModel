# Lewatit VP OC 1065 sorbent parameters

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
                        q_star_CO2 = (gas, q_H2O, params) -> activation * Toth_WADST_isotherm_CO2_wet(gas, q_H2O, params),
                        q_star_N2  = (gas, params) -> zero(gas.p_N2)
                    )
