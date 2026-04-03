# APDES sorbent parameters

function ISOTHERM_PARAMS_APDES(T=Float64)
    IsothermParams(T;
                            # DRY SITE PARAMS
                            q∞     = 2.38,
                            b₀     = 0.0707,
                            ΔH₀    = -5.7047e4,
                            τ₀     = 0.4148,
                            α      = -1.606,
                            T_ref  = 296,

                            # GAB PARAMS
                            Cm    = 36.48,
                            Cg    = 0.1489,
                            K_ads = 0.5751,

                            # Stampi model
                            γ = 0.0016,
                            β = 59.1
                        )
end

SORB_PARAMS_APDES(T=Float64) = SorbentParams(T;
                        ε_bed = 0.092,
                        ε_total = 0.965,
                        dₚ = 0.0075,
                        ρ_bed = 55.4,
                        ρ_particle = 61.0,
                        k_CO2 = 2.0e-4,
                        k_H2O = 2.0e-3,
                        k_N2 = 0.0,
                        ΔH_CO2 = -57000,
                        ΔH_H2O = -49000,
                        ΔH_N2 = 0.0,
                        Dₘ = 4.3e-6,
                        C_solid = 2070,
                        isotherm_params = ISOTHERM_PARAMS_APDES(T),
                        q_star_H2O = GAB_isotherm_H2O,
                        q_star_CO2 = Toth_isotherm_CO2_modified_H2O_Stampi,
                        q_star_N2  = (gas, params) -> zero(gas.p_N2)
                    )
