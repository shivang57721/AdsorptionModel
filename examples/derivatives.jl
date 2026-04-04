# How to run Forward-mode AD through simulate_process
using ForwardDiff
using DifferentialEquations
using AdsorptionModel

col_params  = COL_PARAMS_YOUNG()
sorb_params = SORB_PARAMS_LEWATIT()

# Define a function that takes a parameter vector
function f(x)
    output = simulate_process(TVSA();
        # VERY IMPORTANT: Specify the Dual type and change the linear solver to be AD compatible
        T = eltype(x),
        solver = FBDF(linsolve=SparspakFactorization(reuse_symbolic=false)),

        col_params = col_params,
        sorb_params = sorb_params,

        steps = Dict(
            Adsorption => StepConfig(SaturationLimit(x[1]);
                u_feed     = x[6],
                T_feed     = 288.15,
                T_amb      = 288.15,
                P_out      = 1.01325e5,
                y_CO2_feed = 0.0004,
                y_H2O_feed = x[9]),

            Heating => StepConfig(HeatingUntilTarget(x[3]);
                T_amb   = x[4],
                P_out   = x[8],
                ΔT_heat = 5.0,
                y_H2O_feed = x[9]),

            Desorption => StepConfig(SaturationLimit(x[2]);
                u_feed     = x[7],
                T_feed     = x[5],
                T_amb      = x[4],
                P_out      = x[8],
                y_H2O_feed = x[9]),
        ),
    )

    return output.productivity_kg_per_h_per_m3_bed
end

g = ForwardDiff.gradient(f, [0.9, 0.3, 0.5, 373, 373, 7.06e-2, 7.06e-2/3, 0.2e5, 0.5])
