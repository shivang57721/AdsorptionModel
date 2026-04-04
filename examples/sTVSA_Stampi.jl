using AdsorptionModel
using Plots
using Revise

# ------------------------------------------------------------
# 1. Choose column and sorbent parameters
# ------------------------------------------------------------
col_params  = COL_PARAMS_STAMPI()
sorb_params = SORB_PARAMS_APDES()

phys_params = PhysicalParams(;
                    Cₚ_CO2   = 42.46,
                    Cₚ_H2O   = 73.1,
                    Cₚ_N2    = 29.171,
                    λ_ev     = 2400,
                    Cₚ_steam = 1.89,
                )

# ------------------------------------------------------------
# 2. Run an sTVSA simulation with fixed step durations
#    (sTVSA adds a Preheating step to avoid water condensation)
# ------------------------------------------------------------

u_feed_ads = 50e-6 / (π * col_params.Rᵢ^2)
u_feed_des = 25e-6 / (π * col_params.Rᵢ^2)

@time output = simulate_process(STVSA();
    N = 10,
    col_params  = col_params,
    sorb_params = sorb_params,
    phys_params = phys_params,

    steps = Dict(
        Adsorption => StepConfig(FixedDuration(13_772.0);
            u_feed     = u_feed_ads,
            T_feed     = 293.15,
            T_amb      = 293.15,
            P_out      = 1e5,
            y_CO2_feed = 0.0004,
            y_H2O_feed = 0.0115),

        # Override preheating defaults: disable it for this case
        Preheating => StepConfig(FixedDuration(0.0);
            T_amb      = 95 + 273,
            P_out      = 50e2,
            y_H2O_feed = 1.0,
            ΔT_heat    = 5.0),

        Heating => StepConfig(FixedDuration(704.0);
            T_amb   = 95 + 273,
            P_out   = 50e2,
            y_H2O_feed = 1.0),

        Desorption => StepConfig(FixedDuration(20_000.0);
            u_feed     = u_feed_des,
            T_feed     = 95 + 273,
            T_amb      = 95 + 273,
            P_out      = 50e2,
            y_H2O_feed = 1.0),

        # Disable cooling for this case
        Cooling => StepConfig(FixedDuration(0.0);
            T_amb = 293.15,
            P_out = 50e2),
    ),

    enable_logging = true,
    enable_plotting = true,
    plotter = Plots,
)

output.plots[7]