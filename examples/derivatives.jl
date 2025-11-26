# How to run Forward-mode AD through simulate_process
using ForwardDiff
using DifferentialEquations
using AdsorptionModel

col_params  = COL_PARAMS_YOUNG()
sorb_params = SORB_PARAMS_LEWATIT()

# Define a function that takes a parameter vector
function f(x)
    output = simulate_process(
        # VERY IMPORTANT: Specify the Dual type and change the linear solver to be AD compatible
        T = eltype(x),
        solver = FBDF(linsolve=SparspakFactorization(reuse_symbolic=false)), 

        col_params = col_params,
        sorb_params = sorb_params,

        # --- Step termination limits ---
        q_CO2_saturation_limit_adsorption = x[1],
        q_CO2_saturation_limit_desorption = x[2],
        extra_heating_ratio               = x[3],

        # --- Operating temperatures ---
        T_amb_desorption  = x[4],
        T_feed_desorption = x[5],

        # --- Inlet velocities ---
        u_feed_adsorption = x[6],
        u_feed_desorption = x[7],

        # --- Outlet pressures ---
        P_out_desorption = x[8],

        # --- Humidity ---
        y_H2O_desorption = x[9],
    )

    return output.productivity_kg_per_h_per_m3_bed
end

g = ForwardDiff.gradient(f, [0.9, 0.3, 0.5, 373, 373, 7.06e-2, 7.06e-2/3, 0.2e5, 0.5])
