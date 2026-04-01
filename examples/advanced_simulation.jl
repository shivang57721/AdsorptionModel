using AdsorptionModel
using JLD2
using Plots
using GridVisualize
using VoronoiFVM
using HDF5
using Revise

# ------------------------------------------------------------
# 1. Define model resolution and parameters
# ------------------------------------------------------------
N = 10   # number of finite volumes

phys_params = PHYS_PARAMS_DEFAULT()
col_params  = COL_PARAMS_YOUNG()
sorb_params = SORB_PARAMS_LEWATIT(activation=0.8) # 80% activation
cost_params = COST_PARAMS_DEFAULT()

# ------------------------------------------------------------
# 2. Run an sTVSA simulation
# ------------------------------------------------------------
# The simulate_process function returns a named tuple 'output' with key performance metrics.
# If save_solution=true, then the following information is additionally saved:
# - all_solutions: vector of solutions for each cycle
# - output: key performance indicators (productivity, purity, energy, etc.)
# - index_data: indices for species/state variables
# - sys: VoronoiFVM system struct (mesh + DAE)
# - cycle_steps: metadata for cycle structure
#
# Use keyword arguments to control terminations and operating conditions.

output = simulate_process(
    N = N,
    col_params = col_params,
    sorb_params = sorb_params,
    phys_params = phys_params,
    cost_params = cost_params,

    # Termination criteria
    q_CO2_saturation_limit_adsorption = 0.95,
    q_CO2_saturation_limit_desorption = 0.30,
    ΔT_preheat = 5,
    extra_heating_ratio = 0.20,
    ΔT_heat = 5,
    T_safe_cooling = 343,

    # Operating conditions
    T_amb_desorption  = 373,
    T_feed_desorption = 373,
    u_feed_adsorption = 7.06e-2,
    u_feed_desorption = 7.06e-2 / 3,
    P_out_desorption  = 0.2e5,
    y_H2O_desorption  = 1.0,

    # Logging & saving
    enable_logging = true,
    save_solution = true,
    save_filepath = "examples/simulation_result.h5",
)

# ------------------------------------------------------------
# Plot results for the last cycle at column exit
# ------------------------------------------------------------
h5open("examples/simulation_result.h5", "r") do fid
    output_saved = Dict(k => read_attribute(fid["output"], k) for k in keys(attributes(fid["output"])))

    step_order = ["Adsorption", "Preheating", "Heating", "Desorption", "Cooling", "Pressurization"]
    all_times = Float64[]
    co2_exit  = Float64[]
    step_boundaries = Float64[]

    iCO2 = read_attribute(fid["species"], "CO2")
    data = read(fid, "cycles")
    last_cycle = length(data)

    time_offset = 0.0
    for step_name in step_order
        path = "cycles/cycle_$(last_cycle)/$(step_name)"
        t    = read(fid, path * "/t")
        data = read(fid, path * "/data")   # [n_species, n_nodes, n_times]

        append!(all_times, t .+ time_offset)
        append!(co2_exit,  data[iCO2, end, :])
        time_offset += t[end]

        push!(step_boundaries, time_offset)
    end
    pop!(step_boundaries)

    p = plot(all_times, co2_exit, xlabel="Time (s)", ylabel="CO2 concentration [mol/m³]", label="Last cycle")
    vline!(p, step_boundaries, linestyle=:dash, color=:black, alpha=0.5, label="")
end

# ------------------------------------------------------------
# Plot spatial profiles (e.g., at the end of Adsorption)
# ------------------------------------------------------------

h5open("examples/simulation_result.h5", "r") do fid
    xs = read(fid, "grid/x_coords")
    iCO2 = read_attribute(fid["species"], "CO2")
    data = read(fid, "cycles")
    last_cycle = length(data)
    co2_profile = read(fid, "cycles/cycle_$last_cycle/Adsorption/data")[iCO2, :, end]
    plot(xs, co2_profile,
        xlabel="Column position [m]",
        ylabel="CO₂ concentration [mol/m³]",
        title="CO₂ Spatial Profile at End of Adsorption")
end

# ---------------------------------------
# NOTE: How to access a transient solution
# ---------------------------------------
# tsol[:,:,it] returns the solution for timestep i
# tsol[ispec,:,it] returns the solution for component ispec at timestep i
# tsol(t) returns a (linearly) interpolated solution value for t.
# tsol.t[it] is the corresponding time
# tsol[ispec,ix,it] refers to solution of component ispec at node ix at moment it