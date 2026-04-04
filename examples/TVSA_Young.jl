using AdsorptionModel
using Plots
using Revise

col_params  = COL_PARAMS_YOUNG()
sorb_params = SORB_PARAMS_LEWATIT()

@time output = simulate_process(TVSA();
    N = 10,
    col_params = col_params,
    sorb_params = sorb_params,

    steps = Dict(
        # --- Adsorption: air flows through column ---
        Adsorption => StepConfig(FixedDuration(8_000.0);
            u_feed     = 7.06e-2,       # m/s
            T_feed     = 288.15,
            T_amb      = 288.15,
            P_out      = 1.01325e5,     # atmospheric
            y_CO2_feed = 0.0004,
            y_H2O_feed = 0.00935),

        # --- Heating: wall heating, no flow ---
        Heating => StepConfig(FixedDuration(2_400.0);
            T_amb  = 373.15,
            P_out  = 0.2e5),

        # --- Desorption: vacuum, no inflow ---
        Desorption => StepConfig(FixedDuration(20_000.0);
            T_amb  = 373.15,
            T_feed = 373.15,
            P_out  = 0.2e5),
    ),

    enable_logging = true,
    enable_plotting = true,
    plotter = Plots,
)

output.plots[7]