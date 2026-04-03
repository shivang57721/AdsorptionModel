# Default parameter factory functions for physical, column, and cost parameters

# --------------- PHYSICAL PARAMS ----------------------
PHYS_PARAMS_DEFAULT(T=Float64) = PhysicalParams(T;
                                R        = 8.314,
                                γ₁       = 0.7,
                                γ₂       = 0.5,
                                μ        = 1.789e-5,
                                Cₚ_CO2   = 37.520,
                                Cₚ_H2O   = 34.2,
                                Cₚ_N2    = 29.171,
                                Cₚ_wall  = 4.0e6,
                                λ_ev     = 2350.0,
                                Cₚ_steam = 2.08,
                                MW_CO2   = 44.0,
                                MW_N2    = 28.0,
                                MW_H2O   = 18.0
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

COL_PARAMS_STAMPI(T=Float64) = ColumnParams(T;
                        Rᵢ = 0.08,
                        Rₒ = 0.082,
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
