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
col_params  = COL_PARAMS_YOUNG()
sorb_params = SORB_PARAMS_LEWATIT(activation=0.8) # 80% activation

# ------------------------------------------------------------
# 2. Run a sTVSA simulation with saturation-based step termination
# ------------------------------------------------------------
# The simulate_process function returns a named tuple 'output' with key performance metrics.
# If save_solution=true, then the following information is additionally saved:
# - all_solutions: vector of solutions for each cycle
# - output: key performance indicators (productivity, purity, energy, etc.)
# - index_data: indices for species/state variables
# - cycle_steps: metadata for cycle structure
#
# Use StepDuration types to control when each step ends:
#   FixedDuration(t)          — runs for exactly t seconds
#   SaturationLimit(frac)     — runs until CO₂ loading reaches frac of equilibrium
#   HeatingUntilTarget(ratio) — runs until temperature interpolates to target

output = simulate_process(STVSA();
    N = N,
    col_params = col_params,
    sorb_params = sorb_params,

    steps = Dict(
        Adsorption => StepConfig(SaturationLimit(0.95);
            u_feed     = 7.06e-2,
            T_feed     = 288.15,
            T_amb      = 288.15,
            P_out      = 1.01325e5,
            y_CO2_feed = 0.0004,
            y_H2O_feed = 0.00935),

        Heating => StepConfig(HeatingUntilTarget(0.05);
            T_amb   = 373.0,
            P_out   = 0.2e5,
            ΔT_heat = 5.0,
            y_H2O_feed = 1.0),

        Desorption => StepConfig(SaturationLimit(0.05);
            u_feed     = 7.06e-2 / 3,
            T_feed     = 373.0,
            T_amb      = 373.0,
            P_out      = 0.2e5,
            y_H2O_feed = 1.0),
    ),

    enable_logging = true,
    save_solution = true,
    save_filepath = "examples/simulation_result.h5",

    enable_plotting=true,
    plotter=Plots
)

# ------------------------------------------------------------
# Plot results for the last cycle at column exit
# ------------------------------------------------------------
h5open("examples/simulation_result.h5", "r") do fid
    output_saved = Dict(k => read_attribute(fid["output"], k) for k in keys(attributes(fid["output"])))

    step_order = ["Adsorption", "Blowdown", "Preheating", "Heating", "Desorption", "Cooling", "Pressurization"]
    all_times = Float64[]
    qco2_exit  = Float64[]
    step_boundaries = Float64[]

    q_CO2 = read_attribute(fid["species"], "q_CO2")
    data = read(fid, "cycles")
    last_cycle = length(data)

    time_offset = 0.0
    for step_name in step_order
        path = "cycles/cycle_$(last_cycle)/$(step_name)"
        t    = read(fid, path * "/t")
        data = read(fid, path * "/data")   # [n_species, n_nodes, n_times]

        append!(all_times, t .+ time_offset)
        append!(qco2_exit,  data[q_CO2, end, :])
        time_offset += t[end]

        push!(step_boundaries, time_offset)
    end
    pop!(step_boundaries)

    p = plot(all_times, qco2_exit, xlabel="Time (s)", ylabel="q_CO2 [mol/kg]", label="Last cycle")
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
