# NbOFFIVE-1-Ni sorbent parameters

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
                        q_star_CO2 = (gas, q_H2O, params) -> activation * Toth_isotherm_CO2_modified_H2O_Arvind(gas, q_H2O, params),
                        q_star_N2  = (gas, params) -> zero(gas.p_N2)
                    )
