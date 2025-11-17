using Roots
using Statistics

# --- Trapezoidal integration utility ---
function trapz(ts, U)
    Δt = diff(ts)
    vals = (U[1:end-1, :] .+ U[2:end, :]) ./ 2
    if size(vals, 2) == 1
        return (Δt' * vals)[1, 1]
    else
        return vec(Δt' * vals)
    end
end

# --- Helper: Target temperature for water saturation ---
function T_targ(T; P_heat, y_H2O)
    P_heat * y_H2O - Psat_H2O(T - 273.15)
end


# =====================================================
# ===============  POST-PROCESS FUNCTION  ==============
# =====================================================

function post_process(all_solutions, data, sys, cycle_steps, cost_params)
    phys_params = data.phys_params
    col_params  = data.col_params
    sorb_params = data.sorb_params

    # --- Extract last-cycle step solutions ---
    adsorption      = all_solutions[end][Adsorption]
    preheating      = all_solutions[end][Preheating]
    heating         = all_solutions[end][Heating]
    desorption      = all_solutions[end][Desorption]

    adsorption_params = cycle_steps[Adsorption]
    desorption_params = cycle_steps[Desorption]

    cross_section = π * col_params.Rᵢ^2
    total_time = sum(step.duration for step in cycle_steps)

    # --- Flow rate integration setup ---
    tf = TestFunctionFactory(sys)
    T  = testfunction(tf, [1], [2])

    # --- Compute total outflow during desorption ---
    outflow = integrate(sys, T, desorption; rate = false)
    total_outflow = -trapz(outflow.t, reduce(hcat, outflow.u)')

    # =====================================================
    # =============  PURITY & PRODUCTIVITY  ================
    # =====================================================
    iCO2, iN2, iH2O = data.iCO2, data.iN2, data.iH2O

    purity_CO2      = total_outflow[iCO2] / (total_outflow[iCO2] + total_outflow[iN2] + total_outflow[iH2O])
    dry_purity_CO2  = total_outflow[iCO2] / (total_outflow[iCO2] + total_outflow[iN2])

    working_capacity = (
        total_outflow[iCO2] * cross_section -
        desorption_params.u_feed * cross_section * desorption_params.c_CO2_feed * desorption_params.duration
    ) / (col_params.L * cross_section * sorb_params.ρ_bed)  # mol/kg

    productivity_TPD = total_outflow[iCO2] * cross_section * phys_params.MW_CO2 * 1e-6 /
        (cross_section * col_params.L * (1 - sorb_params.ε_bed) * (total_time / 3600) / 24)

    n_cols_per_train = ceil(total_time / adsorption_params.duration)
    total_time_1cycle = n_cols_per_train * adsorption_params.duration / 3600
    productivity = working_capacity * phys_params.MW_CO2 / 1000 / total_time_1cycle * sorb_params.ρ_bed  # kg/h/m³_bed

    mass_sorbent = col_params.L * cross_section * sorb_params.ρ_bed
    CO2_prod_ton_per_year_single_train = working_capacity * 44 / 1000 / total_time_1cycle * cost_params.annual_operating_hours / 1000 * mass_sorbent * n_cols_per_train
    n_trains = ceil(cost_params.target_CO2_ton_per_year / CO2_prod_ton_per_year_single_train)

    # =====================================================
    # ===========  ELECTRICAL ENERGY CONSUMPTION  ==========
    # =====================================================
    DP_ads = trapz(adsorption.t, getindex.(adsorption.u, data.ip, 1) .- adsorption_params.P_out) / adsorption_params.duration
    DP_des = trapz(desorption.t, getindex.(desorption.u, data.ip, 1) .- desorption_params.P_out) / desorption_params.duration

    η_blow = cost_params.efficiency_blower
    γ = cost_params.adiabatic_index

    power_blower_ads = adsorption_params.u_feed * cross_section * DP_ads / (η_blow * 1000)
    power_blower_des = desorption_params.u_feed * cross_section * DP_des / (η_blow * 1000)

    molar_rate_out_des = cross_section *
        (total_outflow[iCO2] + total_outflow[iN2] + total_outflow[iH2O]) / desorption_params.duration

    power_vacuum_pump = (
        molar_rate_out_des * γ * phys_params.R * desorption[data.iT, end, end] /
        (cost_params.eta_VP * (γ - 1)) *
        ((adsorption_params.P_out / desorption_params.P_out) ^ ((γ - 1) / γ) - 1)
    ) / 1000  # kW

    total_electric_power = (
        power_blower_ads * adsorption_params.duration / 3600 +
        (power_blower_des + power_vacuum_pump) * desorption_params.duration / 3600
    ) / total_time_1cycle

    spec_electric_cons = total_electric_power / 1000 /
        (productivity / 3600 * mass_sorbent / sorb_params.ρ_bed)

    # =====================================================
    # =============  THERMAL ENERGY CONSUMPTION  ===========
    # =====================================================
    function wall_heat(step)
        q = -2π * col_params.h_L * col_params.Rᵢ *
            integrate(sys, (f, u, node, data) -> f[1] = u[data.iT] - u[data.iT_wall], step)[1, :]
        return trapz(step.t, q)
    end

    heat_preheat = wall_heat(preheating)
    heat_heating = wall_heat(heating)
    heat_des     = wall_heat(desorption)

    denom = productivity * total_time_1cycle * mass_sorbent / sorb_params.ρ_bed * 1e6
    spec_heat_preheat = heat_preheat / denom
    spec_heat_heating = heat_heating / denom
    spec_heat_des     = heat_des / denom

    # --- Steam and nitrogen heating ---
    mass_steam = desorption_params.y_H2O_feed * desorption_params.P_out *
        desorption_params.u_feed * cross_section / (phys_params.R * desorption_params.T_feed) *
        phys_params.MW_H2O / 1000 * desorption_params.duration

    moles_N2 = desorption_params.P_out * desorption_params.u_feed * cross_section *
        desorption_params.y_N2_feed / (phys_params.R * desorption_params.T_feed) * desorption_params.duration

    T_v = find_zero(T -> T_targ(T; P_heat = desorption_params.P_out, y_H2O = 1), 373)

    heat_steam = mass_steam * (
        phys_params.Cₚ_H2O * 1000 / phys_params.MW_H2O * (T_v - adsorption_params.T_feed) +
        phys_params.λ_ev * 1000 +
        phys_params.Cₚ_steam * (desorption_params.T_feed - T_v) * 1000
    ) / 1000  # kJ

    heat_N2 = moles_N2 * phys_params.Cₚ_N2 * (desorption_params.T_feed - adsorption_params.T_feed) / 1000
    heat_steam += heat_N2

    spec_heat_steam = heat_steam / 1000 / (productivity * total_time_1cycle * mass_sorbent / sorb_params.ρ_bed)  # MJ/kg

    # =====================================================
    # ==================  TOTAL ENERGY  ====================
    # =====================================================
    tot_spec_heat = spec_heat_des + spec_heat_preheat + spec_heat_heating + spec_heat_steam

    total_heat_cons = tot_spec_heat *
        (working_capacity * phys_params.MW_CO2 / 1000 / total_time_1cycle / 3600 * mass_sorbent) * 1000  # kW

    total_cons = total_electric_power + 0.5 * total_heat_cons
    spec_cons  = spec_electric_cons + 0.5 * tot_spec_heat  # MJ_el/kg

    membrane_energy = calculate_membrane_energy(total_outflow; data)

    # --- Return structured results ---
    (;  specific_energy_consumption_MJ_per_kg = spec_cons,
        co2_dry_purity = dry_purity_CO2,
        productivity_kg_per_h_per_m3_bed = productivity,
        
        # outflow totals (mol)
        co2_outflow_total_mol = total_outflow[iCO2],
        n2_outflow_total_mol  = total_outflow[iN2],
        h2o_outflow_total_mol = total_outflow[iH2O],
        co2_purity            = purity_CO2,

        # capacity & productivity
        working_capacity_mol_per_kg = working_capacity,
        productivity_tpd_per_m3_bed = productivity_TPD,

        # timing / cycles
        total_cycle_time_s     = total_time,
        effective_cycle_time_h = total_time_1cycle,
        cycles_completed       = length(all_solutions),
        adsorption_duration_s  = cycle_steps[Adsorption].duration,
        preheating_duration_s  = cycle_steps[Preheating].duration,
        heating_duration_s     = cycle_steps[Heating].duration,
        desorption_duration_s  = cycle_steps[Desorption].duration,
        cooling_duration_s     = cycle_steps[Cooling].duration,

        # pressure / flow
        avg_pressure_drop_ads_Pa = DP_ads,
        avg_pressure_drop_des_Pa = DP_des,
        molar_rate_out_des_mol_per_s = molar_rate_out_des,

        # electrical power (kW)
        power_blower_ads_kW = power_blower_ads,
        power_blower_des_kW = power_blower_des,
        power_vacuum_pump_kW = power_vacuum_pump,
        total_electric_power_kW = total_electric_power,
        specific_electric_energy_MJ_per_kg = spec_electric_cons,

        # thermal (energies and specific)
        heat_preheat_kJ = heat_preheat,
        heat_heating_kJ = heat_heating,
        heat_desorption_kJ = heat_des,
        specific_heat_preheat_MJ_per_kg = spec_heat_preheat,
        specific_heat_heating_MJ_per_kg = spec_heat_heating,
        specific_heat_desorption_MJ_per_kg = spec_heat_des,
        specific_heat_steam_MJ_per_kg = spec_heat_steam,
        total_specific_thermal_MJ_per_kg = tot_spec_heat,
        total_thermal_power_kW = total_heat_cons,

        # combined consumption
        total_energy_consumption_kW = total_cons,
        membrane_energy_MJ_per_kg = membrane_energy,

        # other details
        cross_section_m2 = cross_section,
        mass_sorbent_kg = mass_sorbent,
        mass_steam_kg = mass_steam,
        moles_n2_mol = moles_N2,
        n_cols_per_train = n_cols_per_train,
        n_trains = n_trains,
        column_length_m = col_params.L
        )
end

function calculate_membrane_energy(total_outflow, Tl=303, P=1e5; data)
    iCO2, iN2 = data.iCO2, data.iN2
    n_H2O_gas, n_H2O_liq = gasliq_H2O(total_outflow, Tl, P; data)
    n_total_gas = total_outflow[iCO2] + total_outflow[iN2] + n_H2O_gas
    x_CO2 =  total_outflow[iCO2] / n_total_gas
    membrane_energy = 0.15929655*(x_CO2^(-0.79870615)) # formula from membranes
    return membrane_energy
end

function gasliq_H2O(total_outflow, T, P; data)
    iCO2, iN2, iH2O = data.iCO2, data.iN2, data.iH2O
    n_total = total_outflow[iCO2] + total_outflow[iN2] + total_outflow[iH2O]
    y_H2O = total_outflow[iH2O] / n_total
    y_H2O_eq = Psat_H2O(T - 273.15) / P # maximum gaseous water molar fraction
    y_H2O_gas = min(y_H2O_eq, y_H2O)
    n_cond = (total_outflow[iH2O] - y_H2O_gas * n_total) / (1 - y_H2O_gas)
    n_H2O_gas = total_outflow[iH2O] - n_cond
    n_H2O_liq = n_cond
    return n_H2O_gas, n_H2O_liq
end