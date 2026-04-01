"""
    export_to_hdf5(filepath; all_solutions, cycle_steps, index_data, sys, output)

Export simulation results to an HDF5 file readable by Python via h5py.

# File structure
    grid/
        x_coords              # spatial node positions [m]
    species/                  # attributes: species name → 1-based index
    cycle_steps/
        <StepName>/           # attributes: all OperatingParameters fields
    output/                   # attributes: all post-processed KPIs
    cycles/
        cycle_<N>/
            <StepName>/
                t             # time array [s]
                data          # Float64[n_species, n_nodes, n_timesteps]

# Python usage (h5py)
    import h5py, numpy as np
    with h5py.File("solution.h5", "r") as f:
        x    = f["grid/x_coords"][:]
        iCO2 = f["species"].attrs["CO2"]           # 1-based index
        last = f["cycles/cycle_3/Adsorption"]
        t    = last["t"][:]
        data = last["data"][:]                      # shape: (8, n_nodes, n_times)
        co2_exit    = data[iCO2-1, -1, :]          # column exit vs time
        co2_profile = data[iCO2-1, :, -1]          # spatial profile at last time
"""
function export_to_hdf5(filepath::String; all_solutions, cycle_steps, index_data, sys, output)

    h5open(filepath, "w") do fid

        # --- Grid coordinates ---
        fid["grid/x_coords"] = coordinates(sys.grid)[1, :]

        # --- Species name → 1-based index ---
        g_species = create_group(fid, "species")
        for (name, idx) in index_data.species
            attributes(g_species)[name] = idx
        end

        # --- Operating parameters for each step ---
        g_steps = create_group(fid, "cycle_steps")
        for step in cycle_steps
            g_step = create_group(g_steps, string(step.step_name))
            for field in fieldnames(OperatingParameters)
                val = getfield(step, field)
                attributes(g_step)[string(field)] = field == :step_name ? string(val) : Float64(val)
            end
        end

        # --- Post-processed KPIs ---
        g_output = create_group(fid, "output")
        for (k, v) in pairs(output)
            try
                attributes(g_output)[string(k)] = v
            catch
                # skip non-serialisable entries
            end
        end

        # --- Solution data: one group per cycle, one subgroup per step ---
        g_cycles = create_group(fid, "cycles")
        for (cycle_idx, cycle_solutions) in enumerate(all_solutions)
            g_cycle = create_group(g_cycles, "cycle_$(cycle_idx)")
            for step in cycle_steps
                sol = cycle_solutions[step.step_name]
                g_step = create_group(g_cycle, string(step.step_name))

                g_step["t"] = collect(Float64, sol.t)

                n_species, n_nodes = size(sol.u[1])
                n_times = length(sol.t)
                data = Array{Float64}(undef, n_species, n_nodes, n_times)
                for (it, snap) in enumerate(sol.u)
                    data[:, :, it] = snap
                end
                g_step["data"] = data
            end
        end

    end

    return nothing
end
