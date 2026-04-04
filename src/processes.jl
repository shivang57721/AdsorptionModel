# Process definitions: step sequences, validation, default steps, and callback generation

# ---------------------------------------------------------------------------
# Step sequences — the single source of truth for each process type
# ---------------------------------------------------------------------------

step_sequence(::TVSA)  = (Adsorption, Blowdown, Heating, Desorption, Cooling, Pressurization)
step_sequence(::STVSA) = (Adsorption, Blowdown, Preheating, Heating, Desorption, Cooling, Pressurization)
step_sequence(::PSA)   = (Pressurization, Adsorption, Blowdown, Evacuation)

required_steps(::TVSA)  = (Adsorption, Heating, Desorption)
required_steps(::STVSA) = (Adsorption, Heating, Desorption)
required_steps(::PSA)   = (Adsorption, Evacuation)

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

function validate_steps(process::ProcessType, user_steps)
    seq = step_sequence(process)
    req = required_steps(process)

    for step in req
        haskey(user_steps, step) || error("$process requires a $step step, but it was not provided.")
    end

    for step in keys(user_steps)
        step in seq || error("$process does not use a $step step. Valid steps: $(join(seq, ", "))")
    end

    for (step, config) in pairs(user_steps)
        _validate_step_duration(step, config)
    end
end

function _validate_step_duration(step::StepType, config::StepConfig)
    dur = config.duration
    if step in (Adsorption, Desorption)
        dur isa HeatingUntilTarget &&
            error("$step cannot use HeatingUntilTarget — use FixedDuration or SaturationLimit.")
    elseif step == Heating
        dur isa SaturationLimit &&
            error("Heating cannot use SaturationLimit — use FixedDuration or HeatingUntilTarget.")
    else
        dur isa FixedDuration ||
            error("$step must use FixedDuration (got $(typeof(dur))).")
    end
end

# ---------------------------------------------------------------------------
# Default step generation — fill in transition steps the user didn't specify
# ---------------------------------------------------------------------------

function fill_default_steps(process::ProcessType, user_steps)
    all_steps = Dict{StepType, StepConfig}(pairs(user_steps)...)
    _add_defaults!(process, all_steps)
    return all_steps
end

function _add_defaults!(::TVSA, steps)
    ads = steps[Adsorption]
    des = steps[Desorption]

    if !haskey(steps, Pressurization)
        steps[Pressurization] = StepConfig(FixedDuration(60.0);
            P_out=ads.P_out, T_amb=ads.T_amb, T_feed=ads.T_feed,
            y_CO2_feed=ads.y_CO2_feed, y_H2O_feed=ads.y_H2O_feed)
    end

    if !haskey(steps, Blowdown)
        dur = abs(ads.P_out - des.P_out) < 1.0 ? 0.0 : 30.0
        steps[Blowdown] = StepConfig(FixedDuration(dur);
            P_out=des.P_out, T_amb=ads.T_amb)
    end

    if !haskey(steps, Cooling)
        steps[Cooling] = StepConfig(FixedDuration(5*3600.0);
            P_out=des.P_out, T_amb=ads.T_amb, T_safe_cooling=343.0)
    end
end

function _add_defaults!(::STVSA, steps)
    ads = steps[Adsorption]
    des = steps[Desorption]

    if !haskey(steps, Pressurization)
        steps[Pressurization] = StepConfig(FixedDuration(60.0);
            P_out=ads.P_out, T_amb=ads.T_amb, T_feed=ads.T_feed,
            y_CO2_feed=ads.y_CO2_feed, y_H2O_feed=ads.y_H2O_feed)
    end

    if !haskey(steps, Blowdown)
        dur = abs(ads.P_out - des.P_out) < 1.0 ? 0.0 : 30.0
        steps[Blowdown] = StepConfig(FixedDuration(dur);
            P_out=des.P_out, T_amb=ads.T_amb)
    end

    if !haskey(steps, Preheating)
        steps[Preheating] = StepConfig(FixedDuration(15*3600.0);
            P_out=des.P_out, T_amb=des.T_amb,
            y_H2O_feed=des.y_H2O_feed, ΔT_heat=5.0)
    end

    if !haskey(steps, Cooling)
        steps[Cooling] = StepConfig(FixedDuration(5*3600.0);
            P_out=des.P_out, T_amb=ads.T_amb, T_safe_cooling=343.0)
    end
end

function _add_defaults!(::PSA, steps)
    ads = steps[Adsorption]
    evac = steps[Evacuation]

    if !haskey(steps, Pressurization)
        steps[Pressurization] = StepConfig(FixedDuration(15.0);
            P_out=ads.P_out, T_amb=ads.T_amb, T_feed=ads.T_feed,
            y_CO2_feed=ads.y_CO2_feed, y_H2O_feed=ads.y_H2O_feed)
    end

    if !haskey(steps, Blowdown)
        dur = abs(ads.P_out - evac.P_out) < 1.0 ? 0.0 : 30.0
        steps[Blowdown] = StepConfig(FixedDuration(dur);
            P_out=evac.P_out, T_amb=ads.T_amb)
    end
end

# ---------------------------------------------------------------------------
# Conversion: StepConfig -> OperatingParameters (internal)
# ---------------------------------------------------------------------------

function to_operating_params(step_type::StepType, config::StepConfig)
    OperatingParameters(Float64;
        step_name = step_type,
        u_feed = config.u_feed,
        T_feed = config.T_feed,
        y_CO2_feed = config.y_CO2_feed,
        y_H2O_feed = config.y_H2O_feed,
        T_amb = config.T_amb,
        P_out = config.P_out,
        duration = max_duration(config.duration),
        ΔT_heat = config.ΔT_heat,
        T_safe_cooling = config.T_safe_cooling,
        q_CO2_saturation_limit = config.duration isa SaturationLimit ? config.duration.limit : NaN,
        extra_heating_ratio = config.duration isa HeatingUntilTarget ? config.duration.extra_ratio : NaN,
        λ = 0.11)
end

function build_cycle_steps(process::ProcessType, user_steps)
    validate_steps(process, user_steps)
    all_configs = fill_default_steps(process, user_steps)
    seq = step_sequence(process)
    cycle_steps = Dictionary{StepType, OperatingParameters}(
        collect(seq),
        [to_operating_params(step, all_configs[step]) for step in seq])
    return cycle_steps, all_configs
end

# ---------------------------------------------------------------------------
# Callback generation — dispatched on StepDuration
# ---------------------------------------------------------------------------

function _compute_q_star_CO2(ads_config::StepConfig, sorb_params::SorbentParams)
    y_N2 = 1.0 - ads_config.y_CO2_feed - ads_config.y_H2O_feed
    c_total = ads_config.P_out / (8.314 * ads_config.T_feed)
    gas_feed = (
        T     = ads_config.T_feed,
        p     = ads_config.P_out,
        p_CO2 = ads_config.P_out * ads_config.y_CO2_feed,
        p_H2O = ads_config.P_out * ads_config.y_H2O_feed,
        p_N2  = ads_config.P_out * y_N2,
        c_CO2 = ads_config.y_CO2_feed * c_total,
        c_H2O = ads_config.y_H2O_feed * c_total,
        c_N2  = y_N2 * c_total,
    )
    q_star_H2O = sorb_params.q_star_H2O(gas_feed, sorb_params.isotherm_params)
    q_star_CO2 = sorb_params.q_star_CO2(gas_feed, q_star_H2O, sorb_params.isotherm_params)
    @show q_star_CO2
    return q_star_CO2
end

# FixedDuration: no early-termination callback
make_callback(::FixedDuration, ::StepType, sys, data, args...) = nothing

# SaturationLimit for Adsorption
function make_callback(dur::SaturationLimit, ::Val{Adsorption}, sys, data, q_star_CO2)
    ContinuousCallback(
        (u, _, _) -> minimum(@view reshape(u, sys)[data.iq_CO2, :]) / q_star_CO2 - dur.limit,
        terminate!)
end

# SaturationLimit for Desorption
function make_callback(dur::SaturationLimit, ::Val{Desorption}, sys, data, q_star_CO2)
    ContinuousCallback(
        (u, _, _) -> maximum(@view reshape(u, sys)[data.iq_CO2, :]) / q_star_CO2 - dur.limit,
        terminate!)
end

# HeatingUntilTarget for Heating
function make_callback(dur::HeatingUntilTarget, ::Val{Heating}, sys, data, config::StepConfig)
    T_amb = config.T_amb
    ΔT = config.ΔT_heat
    ratio = dur.extra_ratio
    ContinuousCallback(
        (u, t, integrator) -> begin
            T_start = integrator.p.data.step_params.T_start
            T_target = (1 - ratio) * T_start + ratio * (T_amb - ΔT)
            minimum(@view reshape(u, sys)[data.iT, :]) - T_target
        end, terminate!)
end

# Internal: Preheating callback (always applied for sTVSA)
function _preheating_callback(config::StepConfig, sys, data)
    max_duration(config.duration) ≈ 0 && return nothing
    T_target = find_zero(T -> T_targ(T; P_heat=config.P_out, y_H2O=config.y_H2O_feed), 373) + config.ΔT_heat
    ContinuousCallback(
        (u, _, _) -> max(0, T_target - minimum(@view reshape(u, sys)[data.iT, :])),
        terminate!)
end

# Internal: Cooling callback
function _cooling_callback(config::StepConfig, sys, data)
    isnan(config.T_safe_cooling) && return nothing
    T_safe = config.T_safe_cooling
    ContinuousCallback(
        (u, _, _) -> max(0, maximum(@view reshape(u, sys)[data.iT, :]) - T_safe),
        terminate!)
end

function build_callbacks(process::ProcessType, all_configs, sys, data, sorb_params)
    callbacks = Dict{StepType, Any}()
    seq = step_sequence(process)

    # Compute equilibrium loading for saturation-based callbacks
    ads_config = all_configs[Adsorption]
    q_star_CO2 = _compute_q_star_CO2(ads_config, sorb_params)

    for step in seq
        config = all_configs[step]
        dur = config.duration

        if step == Adsorption && dur isa SaturationLimit
            callbacks[step] = make_callback(dur, Val(Adsorption), sys, data, q_star_CO2)
        elseif step == Desorption && dur isa SaturationLimit
            callbacks[step] = make_callback(dur, Val(Desorption), sys, data, q_star_CO2)
        elseif step == Heating && dur isa HeatingUntilTarget
            callbacks[step] = make_callback(dur, Val(Heating), sys, data, config)
        elseif step == Preheating
            callbacks[step] = _preheating_callback(config, sys, data)
        elseif step == Cooling
            callbacks[step] = _cooling_callback(config, sys, data)
        else
            callbacks[step] = nothing
        end
    end

    return callbacks
end
