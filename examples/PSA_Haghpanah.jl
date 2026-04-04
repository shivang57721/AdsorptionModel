using AdsorptionModel
using Plots
using Revise

# ------------------------------------------------------------
# 1. Choose column and sorbent parameters
# ------------------------------------------------------------
col_params  = COL_PARAMS_HAGHPANAH()
sorb_params = SORB_PARAMS_ZEOLITE()

# ------------------------------------------------------------
# 2. Run a PSA simulation (post-combustion capture)
# ------------------------------------------------------------

@time output = simulate_process(PSA();
    N = 10,
    col_params  = col_params,
    sorb_params = sorb_params,

    steps = Dict(
        Adsorption => StepConfig(FixedDuration(15.0);
            u_feed     = 1.0,
            T_feed     = 298.15,
            T_amb      = 298.15, 
            P_out      = 1e5,
            y_CO2_feed = 0.15),

        Blowdown => StepConfig(FixedDuration(30.0);
            T_amb = 298.15,
            P_out = 0.2e5),

        Evacuation => StepConfig(FixedDuration(40.0);
            T_amb = 298.15,
            P_out = 0.1e5),

        Pressurization => StepConfig(FixedDuration(15.0);
            T_amb      = 298.15,
            P_out      = 1e5,
            y_CO2_feed = 0.15),
    ),

    steady_state_tol = 1e-3,
    max_cycles = 10,
    enable_logging = true,
    enable_plotting = true,
    plotter = Plots,
)

output.plots[1]
