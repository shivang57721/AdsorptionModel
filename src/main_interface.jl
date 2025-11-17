# ---------------------------------------
# ACCESSING A TRANSIENT SOLUTION
# ---------------------------------------
# tsol[:,:,it] returns the solution for timestep i
# tsol[ispec,:,it] returns the solution for component ispec at timestep i
# tsol(t) returns a (linearly) interpolated solution value for t.
# tsol.t[it] is the corresponding time
# tsol[ispec,ix,it] refers to solution of component ispec at node ix at moment it

function simulate_DCB(; T::Type=Float64, N::Int=10, op_params, col_params=nothing, sorb_params=nothing, phys_params=nothing, 
                        enable_logging=true, solver=FBDF(), kwargs...)
    if col_params === nothing
        col_params = COL_PARAMS_YOUNG(T)
        sorb_params = SORB_PARAMS_LEWATIT(T)
        phys_params = PHYS_PARAMS_LEWATIT(T)
    end
    sys, data, inival = initialize_system(; T, N, op_params, col_params, sorb_params, phys_params)
    problem = ODEProblem(sys, inival, (0, op_params.duration))
    sol = solve(problem, solver; kwargs...)
    return reshape(sol, sys), data, sys
end

function run_simulation(; T::Type=Float64, N=10, cycle_steps::Dictionary{StepType, OperatingParameters},
        max_cycles=8, col_params::ColumnParams, sorb_params::SorbentParams, phys_params::PhysicalParams,
        steady_state_tol=1e-2, enable_logging=true, solver=FBDF(linsolve=KLUFactorization()), kwargs...)
    
    if enable_logging
        logger = ConsoleLogger(stderr, Logging.Info) 
        global_logger(logger)
    else
        global_logger(NullLogger())  # disables logging
    end

    sys, data, inival = initialize_system(; T, N, op_params=cycle_steps[Adsorption], col_params, sorb_params, phys_params)

    # --- Define callbacks ---

    # Preheating callback
    callback_preheating = nothing
    preheating_params = cycle_steps[Preheating]
    if (preheating_params.y_H2O_feed == 0) || (cycle_steps[Desorption].u_feed == 0)
        preheating_params.duration = 0
    elseif Preheating ∈ keys(cycle_steps)
        T_target = find_zero(T -> T_targ(T; P_heat=preheating_params.P_out, y_H2O=preheating_params.y_H2O_feed), 373) + preheating_params.ΔT_heat
        stop_condition(u, t, integrator) = minimum(@view reshape(u, sys)[data.iT, :]) ≥ T_target
        callback_preheating = DiscreteCallback(stop_condition, terminate!)
    end

    # Cooling callback
    callback_cooling = nothing
    if Cooling ∈ keys(cycle_steps) && cycle_steps[Cooling].T_safe != NaN
        callback_cooling = DiscreteCallback((u, _, _) -> maximum(@view reshape(u, sys)[data.iT, :]) ≤ cycle_steps[Cooling].T_safe, terminate!)
    end

    callbacks = Dict(
        Adsorption=>nothing,
        Preheating=>callback_preheating,
        Heating=>nothing,
        Desorption=>nothing,
        Cooling=>callback_cooling,
        Pressurization=>nothing
    )

    # --- Simulation Loop for Multiple Cycles ---
    all_solutions = Vector{Dict{StepType,TransientSolution}}()

    u_current = copy(inival)
    u_previous = similar(u_current)

    for cycle in 1:max_cycles
        if cycle > 2
            # steady-state check
            rel_diff = maximum(norm.(u_current .- u_previous) ./ norm.(u_previous))

            @info "Cycle $(cycle-1) relative difference: $rel_diff"
            if rel_diff < steady_state_tol
                @info "Steady state reached after $(cycle-1) cycles"
                break
            end
        end

        @info "Running cycle $cycle"
        push!(all_solutions, Dict{StepType,TransientSolution}())

        copyto!(u_previous, u_current)

        for step_params in cycle_steps
            @info "Running step $(step_params.step_name)"

            # update step parameters in-place
            step_params.P_out_start = step_params.step_name == Pressurization ? u_current[data.ip, 1] : u_current[data.ip, end]
            copy_params!(data.step_params, step_params)

            problem = ODEProblem(sys, u_current, (0, step_params.duration))

            # solver choice: keep FBDF by default but allow alternative via param
            step_sol = solve(problem, solver; callback = callbacks[step_params.step_name], dense=false, kwargs...)
            step_params.duration = step_sol.t[end]
            @info "Duration of step: $(step_params.duration)"

            all_solutions[cycle][step_params.step_name] = reshape(step_sol, sys)

            # update state for next step
            copyto!(u_current, step_sol.u[end])
        end
    end

    return all_solutions, data, sys
end

function get_cycle_params(; T::Type=Float64, 
                            duration_adsorption, duration_desorption, duration_heating, 
                            T_amb_adsorption=293, T_feed_adsorption=293,
                            T_amb_desorption, T_feed_desorption, 
                            u_feed_adsorption, u_feed_desorption, 
                            P_out_adsorption=1e5, P_out_desorption, 
                            ΔT_heat=5, max_duration_preheating=15*3600,
                            T_safe_cooling=343, max_duration_cooling=3600*5,
                            y_H2O_air=0.0115, y_CO2_air=0.0004, y_H2O_desorption)
    adsorption = OperatingParameters(T;
            step_name = Adsorption,
            u_feed = u_feed_adsorption,
            T_feed = T_feed_adsorption,
            y_CO2_feed = y_CO2_air,
            y_H2O_feed = y_H2O_air,
            T_amb = T_amb_adsorption,
            P_out = P_out_adsorption,
            duration = duration_adsorption)

    preheating = OperatingParameters(T;
            step_name = Preheating,
            T_amb = T_amb_desorption,
            P_out = P_out_desorption,
            ΔT_heat = ΔT_heat,
            y_H2O_feed = y_H2O_desorption,
            duration = max_duration_preheating)

    heating = OperatingParameters(T;
            step_name = Heating,
            T_amb = T_amb_desorption,
            P_out = P_out_desorption,
            y_H2O_feed = y_H2O_desorption,
            duration = duration_heating)
            
    desorption = OperatingParameters(T;
            step_name = Desorption,
            u_feed = u_feed_desorption,
            T_feed = T_feed_desorption,
            y_CO2_feed = 0.0,
            y_H2O_feed = y_H2O_desorption,
            T_amb = T_amb_desorption,
            P_out = P_out_desorption,
            duration = duration_desorption)

    cooling = OperatingParameters(T;
            step_name = Cooling,
            T_amb = T_amb_adsorption,
            P_out = P_out_desorption,
            T_safe = T_safe_cooling,
            duration = max_duration_cooling)

    pressurization = OperatingParameters(T;
                    step_name = Pressurization,
                    T_amb = T_amb_adsorption,
                    P_out = P_out_adsorption,
                    T_feed = T_feed_adsorption,
                    y_CO2_feed = y_CO2_air,
                    y_H2O_feed = y_H2O_air,
                    duration = 60)

    cycle_steps = Dictionary{StepType, OperatingParameters}(
                            [Adsorption, Preheating, Heating, Desorption, Cooling, Pressurization],
                            [adsorption, preheating, heating, desorption, cooling, pressurization])
                            
    return cycle_steps
end

function simulate_process(; T::Type=Float64, N=10,
                            duration_adsorption, duration_desorption, duration_heating, 
                            T_amb_adsorption=293, T_feed_adsorption=293,
                            T_amb_desorption, T_feed_desorption, 
                            u_feed_adsorption, u_feed_desorption, 
                            P_out_adsorption=1e5, P_out_desorption, 
                            ΔT_heat=5, max_duration_preheating=15*3600,
                            T_safe_cooling=343, max_duration_cooling=3600*5,
                            y_H2O_air=0.0115, y_CO2_air=0.0004, y_H2O_desorption,
                            enable_logging=false, enable_plotting=false, plotter=nothing, steady_state_tol=0.005, max_cycles=6,
                            col_params::ColumnParams,
                            sorb_params::SorbentParams,
                            phys_params::PhysicalParams=PHYS_PARAMS_DEFAULT(), 
                            cost_params::CostParams=COST_PARAMS_DEFAULT(),
                            solver = FBDF(linsolve=KLUFactorization()), kwargs...)
    if enable_logging
        logger = ConsoleLogger(stderr, Logging.Info) 
        global_logger(logger)
    else
        global_logger(NullLogger())  # disables logging
    end

    cycle_steps = get_cycle_params(;T,  duration_adsorption, duration_desorption, duration_heating, 
                                        T_amb_adsorption, T_feed_adsorption,
                                        T_amb_desorption, T_feed_desorption, 
                                        u_feed_adsorption, u_feed_desorption, 
                                        P_out_adsorption, P_out_desorption, 
                                        ΔT_heat, max_duration_preheating,
                                        T_safe_cooling, max_duration_cooling,
                                        y_H2O_air, y_CO2_air, y_H2O_desorption)

    all_solutions, index_data, sys = all_solutions, index_data, sys = run_simulation(
            ; T, N,
            cycle_steps,
            col_params,
            sorb_params,
            phys_params,
            steady_state_tol,
            max_cycles,
            enable_logging,
            solver,
            kwargs...
            )
            
    output = post_process(all_solutions, index_data, sys, cycle_steps, cost_params)

    @info "=" ^ 50
    @info "SIMULATION RESULTS"
    @info "=" ^ 50

    for (k,v) in pairs(output)
        @info "$k: $v"
    end

    if enable_plotting
        # cycle_number = length(all_solutions)
        for cycle_number in 1:length(all_solutions)
            if !isdir("plots/cycle_$cycle_number")
                mkdir("plots/cycle_$cycle_number")
            end
            plots = plot_grid_idx(all_solutions; plotter, cycle_steps, cycle_number, grid_idx=10)
            for (name, i) in index_data.species
                plotter.savefig(plots[i], "plots/cycle_$cycle_number/$name.png")
            end
        end
    end

    if output.co2_dry_purity < 0 || output.specific_energy_consumption_MJ_per_kg < 0
        throw(ErrorException("Unknown bug with model, computed solution was negative."))
    end

    return output
end

function is_feasible(T_amb_desorption, T_feed_desorption, y_H2O_desorption, P_out_desorption; ΔT_heat = 5)
    if y_H2O_desorption == 0
        return T_amb_desorption ≥ T_feed_desorption
    end
    T_saturation = find_zero(T -> T_targ(T; P_heat=P_out_desorption, y_H2O=y_H2O_desorption), 373)

    return (T_feed_desorption ≥ T_saturation) &&
            (T_amb_desorption ≥ T_feed_desorption) &&
            (T_amb_desorption ≥ T_saturation + ΔT_heat)
end