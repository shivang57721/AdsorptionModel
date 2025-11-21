using AdsorptionModel
using JLD2
using Plots
using GridVisualize
using VoronoiFVM

# ------------------------------------------------------------
# 1. Define model resolution and parameters
# ------------------------------------------------------------
N = 10   # number of finite volumes

phys_params = PHYS_PARAMS_DEFAULT()
col_params  = COL_PARAMS_YOUNG()
sorb_params = SORB_PARAMS_LEWATIT()
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
    u_feed_desorption = 7.06e-2 / 2,
    P_out_desorption  = 0.2e5,
    y_H2O_desorption  = 1.0,

    # Logging & saving
    enable_logging = true,
    save_solution = true,
    save_filepath = "simulation_result.jld2",
)

# ------------------------------------------------------------
# 3. Extracting key performance indicators
# ------------------------------------------------------------
productivity         = output.productivity_kg_per_h_per_m3_bed
co2_dry_purity       = output.co2_dry_purity
specific_consumption = output.specific_energy_consumption_MJ_per_kg

# ------------------------------------------------------------
# 4. Load the saved solution
# ------------------------------------------------------------
simulation_result = load("simulation_result.jld2")

all_solutions = simulation_result["all_solutions"]
index_data    = simulation_result["index_data"]
sys           = simulation_result["sys"]
cycle_steps   = simulation_result["cycle_steps"]

# ------------------------------------------------------------
# 5. Plot results for the last cycle at column exit
# ------------------------------------------------------------
# Plot all variables at column exit (grid_idx = N)
plots = plot_grid_idx(
    all_solutions;
    cycle_steps,
    cycle_number = length(all_solutions),
    grid_idx = N,
    plotter = Plots,
)

# Example: Plot CO₂ gas-phase mole fraction and adsorbed-phase CO₂ loading
plot(plots[index_data.iCO2], title = "CO₂ Mole Fraction at Column Exit")
plot(plots[index_data.iq_CO2], title = "Adsorbed CO₂ Loading at Column Exit")

# ------------------------------------------------------------
# 6. Plot spatial profiles (e.g., at the end of Adsorption)
# ------------------------------------------------------------
# Choose the adsorption step of the last cycle
sol_ads = all_solutions[end][Adsorption]

# Option 1: Get grid and co2 profile and plot directly
xs = coordinates(sys.grid)[1,:]
co2_profile = sol_ads[index_data.iCO2, :, end]

plot(xs, co2_profile,
     xlabel="Column position [m]",
     ylabel="CO₂ concentration [mol/m³]",
     title="CO₂ Spatial Profile at End of Adsorption")

# Option 2: Use scalarplot from GridVisualize.jl
scalarplot(
    sys,
    sol_ads[end];
    species = index_data.iCO2,
    Plotter = Plots,
    title = "CO₂ Spatial Profile at End of Adsorption",
)
xlabel!("Column position [m]")
ylabel!("CO₂ concentration [mol/m³]")