# High-level simulation API

# ---------------------------------------
# ACCESSING A TRANSIENT SOLUTION
# ---------------------------------------
# tsol[:,:,it] returns the solution for timestep i
# tsol[ispec,:,it] returns the solution for component ispec at timestep i
# tsol(t) returns a (linearly) interpolated solution value for t.
# tsol.t[it] is the corresponding time
# tsol[ispec,ix,it] refers to solution of component ispec at node ix at moment it

# ---------------------------------------------------------------------------
# Single-step simulation (Dynamic Column Breakthrough)
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Core cycle loop (internal)
# ---------------------------------------------------------------------------

function _run_cycle_loop(; sys, data, inival, cycle_steps, callbacks,
                          max_cycles, steady_state_tol,
                          solver=FBDF(linsolve=KLUFactorization()), kwargs...)

    all_solutions = Vector{Dict{StepType,TransientSolution}}()

    u_current = copy(inival)
    u_previous = similar(u_current)
    t_current = [step.duration for step in cycle_steps]
    t_previous = similar(t_current)

    for cycle in 1:max_cycles
        if cycle > 2
            u_rel_diff = maximum(norm.(u_current .- u_previous) ./ (norm.(u_previous) .+ eps(eltype(u_previous))))
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
            step_params.P_out_start = step_params.step_name == Pressurization ? u_current[data.ip, 1] : u_current[data.ip, end]
            step_params.T_start     = minimum(u_current[data.iT, :])
            copy_params!(data.step_params, step_params)

            problem = ODEProblem(sys, u_current, (0, step_params.duration))

            cb = get(callbacks, step_params.step_name, nothing)
            step_sol = solve(problem, solver; callback=cb, dense=false, kwargs...)
            step_params.duration = step_sol.t[end]
            @info "Duration of step: $(step_params.duration)"

            all_solutions[cycle][step_params.step_name] = reshape(step_sol, sys)

            copyto!(u_current, step_sol.u[end])
        end
        t_current = [step.duration for step in cycle_steps]
    end

    return all_solutions
end

# ---------------------------------------------------------------------------
# Main entry point: simulate_process(process::ProcessType; ...)
# ---------------------------------------------------------------------------

function simulate_process(process::ProcessType;
        T::Type=Float64, N=10,
        steps,
        col_params::ColumnParams,
        sorb_params::SorbentParams,
        phys_params::PhysicalParams=PHYS_PARAMS_DEFAULT(),
        cost_params::CostParams=COST_PARAMS_DEFAULT(),
        steady_state_tol=0.005, max_cycles=6,
        enable_logging=false, save_solution=false, save_filepath=nothing,
        enable_plotting=false, plotter=nothing,
        solver=FBDF(linsolve=KLUFactorization()), kwargs...)

    if enable_logging
        global_logger(ConsoleLogger(stderr, Logging.Info))
    else
        global_logger(NullLogger())
    end

    # Build internal representations
    cycle_steps, all_configs = build_cycle_steps(process, steps)

    # Initialize PDE system from adsorption step
    sys, data, inival = initialize_system(; T, N,
        op_params=cycle_steps[Adsorption], col_params, sorb_params, phys_params)

    # Build callbacks from StepDuration types
    callbacks = build_callbacks(process, all_configs, sys, data, sorb_params)

    # Run the cycle loop
    all_solutions = _run_cycle_loop(;
        sys, data, inival, cycle_steps, callbacks,
        max_cycles, steady_state_tol, solver, kwargs...)

    # Post-process
    output = post_process(process, all_solutions, data, sys, cycle_steps, cost_params)

    @info "=" ^ 50
    @info "SIMULATION RESULTS"
    @info "=" ^ 50
    for (k, v) in pairs(output)
        @info "$k: $v"
    end

    if save_solution
        export_to_hdf5(save_filepath; all_solutions, cycle_steps, index_data=data, sys, output)
    end

    plots = nothing
    if enable_plotting
        cycle_number = length(all_solutions)
        plots = plot_grid_idx(all_solutions; plotter, cycle_steps, cycle_number, grid_idx=N)
    end

    return (; output..., plots)
end

# ---------------------------------------------------------------------------
# Feasibility check
# ---------------------------------------------------------------------------

function is_feasible(q_CO2_saturation_limit_adsorption, q_CO2_saturation_limit_desorption,
                    T_amb_desorption, T_feed_desorption, y_H2O_desorption, P_out_desorption; ΔT_heat = 5)
    if y_H2O_desorption == 0
        return T_amb_desorption ≥ T_feed_desorption
    end
    T_saturation = find_zero(T -> T_targ(T; P_heat=P_out_desorption, y_H2O=y_H2O_desorption), 373)

    return (T_feed_desorption ≥ T_saturation) &&
            (T_amb_desorption ≥ T_feed_desorption) &&
            (T_amb_desorption ≥ T_saturation + ΔT_heat) &&
            (q_CO2_saturation_limit_adsorption > q_CO2_saturation_limit_desorption)
end
