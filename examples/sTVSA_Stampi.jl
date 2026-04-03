using AdsorptionModel
using Plots
using Revise

# ------------------------------------------------------------
# 1. Choose column and sorbent parameters
# ------------------------------------------------------------
col_params  = COL_PARAMS_STAMPI()
sorb_params = SORB_PARAMS_APDES()

phys_params = PhysicalParams(;
                    R        = 8.314,
                    γ₁       = 0.7,
                    γ₂       = 0.5,
                    μ        = 1.789e-5,
                    Cₚ_CO2   = 42.46,
                    Cₚ_H2O   = 73.1,
                    Cₚ_N2    = 29.171,
                    Cₚ_wall  = 4.0e6,
                    λ_ev     = 2400,
                    Cₚ_steam = 1.89,
                    MW_CO2   = 44.0,
                    MW_N2    = 28.0,
                    MW_H2O   = 18.0
                )

# ------------------------------------------------------------
# 2. Run a simple sTVSA simulation by specifying the step durations
# ------------------------------------------------------------

@time output = simulate_process(
    N = 10,
    col_params = col_params,
    sorb_params = sorb_params,
    phys_params = phys_params,

    # --- Step durations (in seconds) ---
    duration_adsorption = 13_772,
    duration_heating    = 704,
    duration_desorption = 30_000,

    max_duration_preheating = 0,
    max_duration_cooling = 0,

    # --- Operating temperatures ---
    T_amb_adsorption  = 293.15,
    T_feed_adsorption = 293.15,
    T_amb_desorption  = 95+273,
    T_feed_desorption = 95+273,

    # --- Inlet velocities ---
    u_feed_adsorption = 50e-6 / (π * col_params.Rᵢ^2),
    u_feed_desorption = 25e-6 / (π * col_params.Rᵢ^2),

    # --- Outlet pressures ---
    P_out_adsorption = 1e5,   # atmospheric pressure
    P_out_desorption = 50e2,  # vacuum level

    # --- Humidity ---
    y_H2O_adsorption = 0.0115,     # typical ambient humidity
    y_H2O_desorption = 1.0,        # steam

    # --- Logging ---
    enable_logging = true,
    reltol=1e-10,

    # --- Plotting ---
    enable_plotting = true,
    plotter = Plots
)

output.plots[7]