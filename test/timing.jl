using AdsorptionModel
using DifferentialEquations
using BenchmarkTools
using Plots
using Revise

col_params  = COL_PARAMS_YOUNG()
sorb_params = SORB_PARAMS_LEWATIT()

function run_baseline(N)
    out = simulate_process(;
        N,
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
        enable_logging = false,
        reltol=1e-6
    )
    return (; N, out...)
end

run_baseline(10)

results = Dict(:N => [], :time => [], :bytes => [], :gctime => [])

function add_to_dict(out)
    push!(results[:N], out.value.N)
    push!(results[:time], out.time)
    push!(results[:bytes], out.bytes)
    push!(results[:gctime], out.gctime)
end

out = @btimed run_baseline(5)
add_to_dict(out)
out = @btimed run_baseline(10)
add_to_dict(out)
out = @btimed run_baseline(50)
add_to_dict(out)
out = @btimed run_baseline(100)
add_to_dict(out)
out = @btimed run_baseline(150)
add_to_dict(out)
out = @btimed run_baseline(200)
add_to_dict(out)

plot(results[:N], results[:time], marker='o', xlabel="N", ylabel="Computational time (s)", label=nothing)
plot(results[:N], results[:bytes] ./ 10^6, marker='o', xlabel="N", ylabel="Memory (MB)", label=nothing)
