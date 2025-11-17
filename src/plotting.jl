# tsol[ispec,ix,it] refers to solution of component ispec at node ix at moment it

function plot_grid_idx(all_solutions; plotter, cycle_steps, cycle_number, grid_idx, title="", kwargs...)
    # Get the cycle number
    cycle_solutions = all_solutions[cycle_number]
    
    # Initialize arrays to collect data across all steps
    all_times = Float64[]
    all_values = Matrix{Float64}[]
    time_offset = 0.0
    
    # Process each step in order
    for step in cycle_steps
        step_name = step.step_name
        sol = cycle_solutions[step_name]
        
        # Extract times and values for this step
        step_times = sol.t
        step_values = sol[:, grid_idx, :]
        
        # Adjust times to be continuous across steps
        adjusted_times = step_times .+ time_offset
        
        # Store data
        append!(all_times, adjusted_times)
        push!(all_values, step_values)
        
        # Update time offset for next step
        time_offset += sol.t[end]
    end
    
    all_values = hcat(all_values...)

    # Require a plot function on the plotter module
    if !isdefined(plotter, :plot)
        error("plotter must provide a `plot` function, e.g. pass `plotter=Plots`")
    end

    # Create the plot using the provided plotter module
    species_dict = Dict(1 => "Pressure [Pa]", 2 => "Concentration of N₂", 3 => "Concetration of CO₂", 4 => "Concetartion of H₂O",
                        5 => "Temperature [K]", 6 => "Temperature of wall [K]", 7 => "q_CO2", 8 => "q_H2O")
    plots = []
    for (i, row) in pairs(eachrow(all_values))
        p = plotter.plot(all_times, row;
            label="Cycle number: $(cycle_number)", 
            xlabel="Time (s)", 
            ylabel="$(species_dict[i])",
            title=title,
            linewidth=1,
            kwargs...)
         # Add vertical lines to separate steps if supported
        time_marker = 0.0
        for (j, step) in enumerate(cycle_steps)
            if j < length(cycle_steps)
                time_marker += step.duration
                if isdefined(plotter, Symbol("vline!"))
                    plotter.vline!(p, [time_marker]; label="", linestyle=:dash, color=:black, alpha=0.5)
                end
            end
        end
        push!(plots, p)
    end
    
    return plots
end