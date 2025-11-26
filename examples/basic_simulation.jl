using AdsorptionModel

# ------------------------------------------------------------
# 1. Choose column and sorbent parameters
# ------------------------------------------------------------
col_params  = COL_PARAMS_ARVIND()
sorb_params = SORB_PARAMS_NbOFFIVE()

# ------------------------------------------------------------
# 2. Run a simple TVSA simulation by specifying the step durations
# ------------------------------------------------------------

output = simulate_process(
    col_params = col_params,
    sorb_params = sorb_params,

    # --- Step durations (in seconds) ---
    duration_adsorption = 8_000,
    duration_heating    = 2_400,
    duration_desorption = 20_000,

    # --- Operating temperatures ---
    T_amb_adsorption  = 288.15,
    T_feed_adsorption = 288.15,
    T_amb_desorption  = 373.15,
    T_feed_desorption = 373.15,

    # --- Inlet velocities ---
    u_feed_adsorption = 7.06e-2,     # m/s
    u_feed_desorption = 0.0,         # no inflow during desorption

    # --- Outlet pressures ---
    P_out_adsorption = 1.01325e5,   # atmospheric pressure
    P_out_desorption = 0.2e5,       # vacuum level

    # --- Humidity ---
    y_H2O_adsorption = 0.00935,     # typical ambient humidity
    y_H2O_desorption = 0.0,         # dry sweep gas / vacuum

    # --- Logging ---
    enable_logging = true,
)