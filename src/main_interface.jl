# ---------------------------------------
# ACCESSING A TRANSIENT SOLUTION
# ---------------------------------------
# tsol[:,:,it] returns the solution for timestep i
# tsol[ispec,:,it] returns the solution for component ispec at timestep i
# tsol(t) returns a (linearly) interpolated solution value for t.
# tsol.t[it] is the corresponding time
# tsol[ispec,ix,it] refers to solution of component ispec at node ix at moment it

function simulate_DCB(; T::Type=Float64, N::Int=10, 
                        op_params::OperatingParameters,
                        col_params::ColumnParams,
                        sorb_params::SorbentParams,
                        phys_params::PhysicalParams=PHYS_PARAMS_DEFAULT(), 
                        solver=FBDF(linsolve=KLUFactorization()), kwargs...)
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

    # Adsorption and Desorption callbacks
    callback_adsorption = nothing
    callback_desorption = nothing
    adsorption_params = cycle_steps[Adsorption]
    desorption_params = cycle_steps[Desorption]
    if adsorption_params.q_CO2_saturation_limit !== NaN && desorption_params.q_CO2_saturation_limit !== NaN
        # Calculate q_star_CO2 from feed conditions
        p_H2O = adsorption_params.P_out * adsorption_params.c_H2O_feed / adsorption_params.c_total_feed
        p_CO2 = adsorption_params.P_out * adsorption_params.c_CO2_feed / adsorption_params.c_total_feed
        q_star_H2O = sorb_params.q_star_H2O(adsorption_params.T_feed, p_H2O, sorb_params.isotherm_params)
        q_star_CO2 = sorb_params.q_star_CO2(adsorption_params.T_feed, p_CO2, q_star_H2O, sorb_params.isotherm_params)

        callback_adsorption = ContinuousCallback(
                                (u,_,_) -> minimum(@view reshape(u, sys)[data.iq_CO2, :]) / q_star_CO2 - adsorption_params.q_CO2_saturation_limit,
                                terminate!)

        callback_desorption = ContinuousCallback(
                                (u,_,_) ->maximum(@view reshape(u, sys)[data.iq_CO2, :]) / q_star_CO2 - desorption_params.q_CO2_saturation_limit,
                                terminate!)
    end

    # Preheating callback
    callback_preheating = nothing
    preheating_params = cycle_steps[Preheating]
    if preheating_params.duration != 0
        T_target = find_zero(T -> T_targ(T; P_heat=preheating_params.P_out, y_H2O=preheating_params.y_H2O_feed), 373) + preheating_params.ΔT_heat
        callback_preheating = ContinuousCallback((u,_,_) -> max(0, T_target - minimum(@view reshape(u, sys)[data.iT, :])), terminate!)
    end

    # Heating callback
    callback_heating = nothing
    heating_params = cycle_steps[Heating]
    if heating_params.extra_heating_ratio !== NaN
        callback_heating = ContinuousCallback((u,t,integrator) -> begin
                                T_target = (1 - heating_params.extra_heating_ratio) * integrator.p.data.step_params.T_start + 
                                                heating_params.extra_heating_ratio * (heating_params.T_amb - heating_params.ΔT_heat)
                                minimum(@view reshape(u, sys)[data.iT, :]) - T_target
                            end, terminate!)
                            
    end

    # Cooling callback
    callback_cooling = nothing
    cooling_params = cycle_steps[Cooling]
    if cooling_params.T_safe_cooling != NaN
        callback_cooling = ContinuousCallback((u, _, _) -> max(0, maximum(@view reshape(u, sys)[data.iT, :]) - cooling_params.T_safe_cooling), terminate!)
    end

    callbacks = Dict(
        Adsorption=>callback_adsorption,
        Preheating=>callback_preheating,
        Heating=>callback_heating,
        Desorption=>callback_desorption,
        Cooling=>callback_cooling,
        Pressurization=>nothing
    )

    # --- Simulation Loop for Multiple Cycles ---
    all_solutions = Vector{Dict{StepType,TransientSolution}}()

    u_current = copy(inival)
    u_previous = similar(u_current)
    t_current = [step.duration for step in cycle_steps]
    t_previous = similar(t_current)

    for cycle in 1:max_cycles
        if cycle > 2
            # steady-state check
            u_rel_diff = maximum(norm.(u_current .- u_previous) ./ norm.(u_previous))
            t_rel_diff = norm(t_current .- t_previous) / norm(t_previous)

            @info "Cycle $(cycle-1) u relative difference: $u_rel_diff"
            @info "Cycle $(cycle-1) t relative difference: $t_rel_diff"
            if u_rel_diff < steady_state_tol && t_rel_diff < 0.01
                @info "Steady state reached after $(cycle-1) cycles"
                break
            end
        end

        @info "Running cycle $cycle"
        push!(all_solutions, Dict{StepType,TransientSolution}())

        copyto!(u_previous, u_current)
        t_previous = copy(t_current)

        for step_params in cycle_steps
            @info "Running step $(step_params.step_name)"
            # update step parameters
            step_params.P_out_start = step_params.step_name == Pressurization ? u_current[data.ip, 1] : u_current[data.ip, end]
            step_params.T_start     = minimum(u_current[data.iT, :])
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
        t_current = [step.duration for step in cycle_steps]
    end

    return all_solutions, data, sys
end

function get_cycle_params(; T::Type=Float64, 
                            duration_adsorption=10*3600, duration_desorption=10*3600, duration_heating=10*3600,
                            q_CO2_saturation_limit_adsorption=NaN,
                            q_CO2_saturation_limit_desorption=NaN,
                            extra_heating_ratio=NaN,
                            T_amb_adsorption=293, T_feed_adsorption=293,
                            T_amb_desorption, T_feed_desorption, 
                            u_feed_adsorption, u_feed_desorption, 
                            P_out_adsorption=1e5, P_out_desorption, 
                            ΔT_preheat=5, ΔT_heat=5,
                            max_duration_preheating=15*3600,
                            T_safe_cooling=343, max_duration_cooling=3600*5,
                            y_H2O_adsorption=0.0115, y_CO2_adsorption=0.0004, y_H2O_desorption)
    adsorption = OperatingParameters(T;
            step_name = Adsorption,
            u_feed = u_feed_adsorption,
            T_feed = T_feed_adsorption,
            y_CO2_feed = y_CO2_adsorption,
            y_H2O_feed = y_H2O_adsorption,
            T_amb = T_amb_adsorption,
            P_out = P_out_adsorption,
            duration = duration_adsorption,
            q_CO2_saturation_limit = q_CO2_saturation_limit_adsorption)

    preheating = OperatingParameters(T;
            step_name = Preheating,
            T_amb = T_amb_desorption,
            P_out = P_out_desorption,
            ΔT_heat = ΔT_preheat,
            y_H2O_feed = y_H2O_desorption,
            duration = max_duration_preheating)
     
    # Cancel preheating if there is no water in the purge feed
    if y_H2O_desorption < 1e-4 || u_feed_desorption == 0
        preheating.duration = 0
    end

    heating = OperatingParameters(T;
            step_name = Heating,
            T_amb = T_amb_desorption,
            P_out = P_out_desorption,
            y_H2O_feed = y_H2O_desorption,
            duration = duration_heating,
            ΔT_heat = ΔT_heat,
            extra_heating_ratio = extra_heating_ratio)
            
    desorption = OperatingParameters(T;
            step_name = Desorption,
            u_feed = u_feed_desorption,
            T_feed = T_feed_desorption,
            y_CO2_feed = 0.0,
            y_H2O_feed = y_H2O_desorption,
            T_amb = T_amb_desorption,
            P_out = P_out_desorption,
            duration = duration_desorption,
            q_CO2_saturation_limit = q_CO2_saturation_limit_desorption)

    cooling = OperatingParameters(T;
            step_name = Cooling,
            T_amb = T_amb_adsorption,
            P_out = P_out_desorption,
            T_safe_cooling = T_safe_cooling,
            duration = max_duration_cooling)

    pressurization = OperatingParameters(T;
                    step_name = Pressurization,
                    T_amb = T_amb_adsorption,
                    P_out = P_out_adsorption,
                    T_feed = T_feed_adsorption,
                    y_CO2_feed = y_CO2_adsorption,
                    y_H2O_feed = y_H2O_adsorption,
                    duration = 60)

    cycle_steps = Dictionary{StepType, OperatingParameters}(
                            [Adsorption, Preheating, Heating, Desorption, Cooling, Pressurization],
                            [adsorption, preheating, heating, desorption, cooling, pressurization])
                            
    return cycle_steps
end

function simulate_process(; T::Type=Float64, N=10,
                            duration_adsorption=10*3600, duration_desorption=10*3600, duration_heating=10*3600,
                            q_CO2_saturation_limit_adsorption=NaN,
                            q_CO2_saturation_limit_desorption=NaN,
                            extra_heating_ratio=NaN,
                            T_amb_adsorption=293, T_feed_adsorption=293,
                            T_amb_desorption, T_feed_desorption, 
                            u_feed_adsorption, u_feed_desorption, 
                            P_out_adsorption=1e5, P_out_desorption, 
                            ΔT_preheat=5, ΔT_heat=5,
                            max_duration_preheating=15*3600,
                            T_safe_cooling=343, max_duration_cooling=5*3600,
                            y_H2O_adsorption=0.0115, y_CO2_adsorption=0.0004, y_H2O_desorption,
                            enable_logging=false, save_solution=false, save_filepath=nothing,
                            steady_state_tol=0.005, max_cycles=6,
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
                                        q_CO2_saturation_limit_adsorption,
                                        q_CO2_saturation_limit_desorption,
                                        ΔT_heat, extra_heating_ratio,
                                        T_amb_adsorption, T_feed_adsorption,
                                        T_amb_desorption, T_feed_desorption, 
                                        u_feed_adsorption, u_feed_desorption, 
                                        P_out_adsorption, P_out_desorption, 
                                        ΔT_preheat, max_duration_preheating,
                                        T_safe_cooling, max_duration_cooling,
                                        y_H2O_adsorption, y_CO2_adsorption, y_H2O_desorption)

    all_solutions, index_data, sys = run_simulation(
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

    # Save the entire solution if requested
    if save_solution
        save(save_filepath, Dict("all_solutions"=>all_solutions, 
                                "cycle_steps"=>cycle_steps,
                                "index_data"=>index_data,
                                "sys"=>sys,
                                "output"=>output))
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
