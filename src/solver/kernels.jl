# Kernel functions for the VoronoiFVM PDE system

_eval_k(k::Real,     _, _) = k
_eval_k(k::Function, c, q_star) = k(c, q_star + 1e-10)

# ---------------------------------------------------------------------------
# Velocity models
# ---------------------------------------------------------------------------

darcy_velocity(u, edge, data) = begin
    @inbounds begin
        return - 1 / data.darcy_k * (u[data.ip, 2] - u[data.ip, 1])
    end
end

ergun_velocity(u, edge, data) = begin
    @inbounds dpdz = (u[data.ip, 2] - u[data.ip, 1]) / edgelength(edge)
    b = data.ergun_b; c = data.ergun_c
    vh = project(edge, -(2c * dpdz) / (b + sqrt(b^2 + 4c * sign(dpdz) * dpdz)))
    return vh
end

# ---------------------------------------------------------------------------
# Flux functions
# ---------------------------------------------------------------------------

function flux_exponential(f, u, edge, data)
    ε_bed = data.sorb_params.ε_bed

    vh = data.velocity(u, edge, data)

    D_L = data.D_L(vh / edgelength(edge))
    bp, bm = (ε_bed * D_L) .* fbernoulli_pm(vh / (ε_bed * D_L))
    f[data.iN2] = bm * u[data.iN2, 1] - bp * u[data.iN2, 2]
    f[data.iCO2] = bm * u[data.iCO2, 1] - bp * u[data.iCO2, 2]
    f[data.iH2O] = bm * u[data.iH2O, 1] - bp * u[data.iH2O, 2]

    # --- Thermal Flux (Energy Balance) ---
    C_gas_1 = u[data.iCO2,1]*data.phys_params.Cₚ_CO2 + u[data.iH2O,1]*data.phys_params.Cₚ_H2O + u[data.iN2,1]*data.phys_params.Cₚ_N2
    C_gas_2 = u[data.iCO2,2]*data.phys_params.Cₚ_CO2 + u[data.iH2O,2]*data.phys_params.Cₚ_H2O + u[data.iN2,2]*data.phys_params.Cₚ_N2
    C_gas_edge = (C_gas_1 + C_gas_2) / 2

    K_L = data.K_L(vh / edgelength(edge), C_gas_edge)
    bp, bm = (ε_bed * K_L) .* fbernoulli_pm(vh * C_gas_edge / (ε_bed * K_L))
    f[data.iT] = bm * u[data.iT, 1] - bp * u[data.iT, 2]
end

function flux_upwind(f, u, edge, data)
    vh = data.velocity(u, edge, data)
    ε_bed = data.sorb_params.ε_bed

    D_L = data.D_L(vh / edgelength(edge))
    f[data.iN2] = vh * u[data.iN2, 1] - ε_bed * D_L * (u[data.iN2,2] - u[data.iN2,1])
    f[data.iCO2] = vh * u[data.iCO2, 1] - ε_bed * D_L * (u[data.iCO2,2] - u[data.iCO2,1])
    f[data.iH2O] = vh * u[data.iH2O, 1] - ε_bed * D_L * (u[data.iH2O,2] - u[data.iH2O,1])

    C_gas_1 = u[data.iCO2,1]*data.phys_params.Cₚ_CO2 + u[data.iH2O,1]*data.phys_params.Cₚ_H2O + u[data.iN2,1]*data.phys_params.Cₚ_N2
    C_gas_2 = u[data.iCO2,2]*data.phys_params.Cₚ_CO2 + u[data.iH2O,2]*data.phys_params.Cₚ_H2O + u[data.iN2,2]*data.phys_params.Cₚ_N2
    C_gas_edge = (C_gas_1 + C_gas_2) / 2

    K_L = data.K_L(vh / edgelength(edge), C_gas_edge)
    f[data.iT] = vh * C_gas_1 * u[data.iT, 1] - ε_bed * K_L * (u[data.iT, 2] - u[data.iT, 1])
end

# ---------------------------------------------------------------------------
# Storage, reaction, boundary conditions, source
# ---------------------------------------------------------------------------

function storage(y, u, node, data)
    @inbounds begin
        ε_total = data.sorb_params.ε_total
        C_solid = data.sorb_params.C_solid
        ρ_bed  = data.sorb_params.ρ_bed
        ΔH_CO2 = data.sorb_params.ΔH_CO2
        ΔH_H2O = data.sorb_params.ΔH_H2O
        ΔH_N2  = data.sorb_params.ΔH_N2

        N2 = u[data.iN2]
        CO2 = u[data.iCO2]
        H2O = u[data.iH2O]

        y[data.iN2]  = ε_total * N2
        y[data.iCO2] = ε_total * CO2
        y[data.iH2O] = ε_total * H2O

        C_gas = CO2 * data.phys_params.Cₚ_CO2 + H2O * data.phys_params.Cₚ_H2O + N2 * data.phys_params.Cₚ_N2
        C_ads = u[data.iq_CO2] * data.phys_params.Cₚ_CO2 + u[data.iq_H2O] * data.phys_params.Cₚ_H2O +
                u[data.iq_N2]  * data.phys_params.Cₚ_N2

        y[data.iT] = (ε_total * C_gas + ρ_bed * C_solid + ρ_bed * C_ads) * u[data.iT] - ε_total * u[data.ip]

        y[data.iT_wall] = u[data.iT_wall]

        # adsorbed species
        y[data.iq_CO2] = u[data.iq_CO2]
        y[data.iq_H2O] = u[data.iq_H2O]
        y[data.iq_N2]  = u[data.iq_N2]

        y[data.iCO2] += ρ_bed * u[data.iq_CO2]
        y[data.iH2O] += ρ_bed * u[data.iq_H2O]
        y[data.iN2]  += ρ_bed * u[data.iq_N2]
        y[data.iT] += ρ_bed * ΔH_CO2 * u[data.iq_CO2] + ρ_bed * ΔH_H2O * u[data.iq_H2O] +
                      ρ_bed * ΔH_N2  * u[data.iq_N2]
    end
end

function reaction(y, u, node, data)
    @inbounds begin
        col_params  = data.col_params
        sorb_params = data.sorb_params
        phys_params = data.phys_params
        step_params = data.step_params

        Rᵢ  = col_params.Rᵢ
        Rₒ  = col_params.Rₒ
        h_L = col_params.h_L
        h_W = col_params.h_W
        Cₚ_wall = phys_params.Cₚ_wall
        a_wall  = π * (Rₒ^2 - Rᵢ^2)

        # Pressure term
        c_total_node = u[data.iN2] + u[data.iCO2] + u[data.iH2O]
        y[data.ip] = u[data.ip] - c_total_node * phys_params.R * u[data.iT]

        # Temperature terms
        y[data.iT] = 2h_L/Rᵢ * (u[data.iT] - u[data.iT_wall])
        y[data.iT_wall] = - 2π/(Cₚ_wall * a_wall) * (h_L * Rᵢ * (u[data.iT] - u[data.iT_wall]) - h_W * Rₒ * (u[data.iT_wall] - step_params.T_amb))

        c_total = u[data.iCO2] + u[data.iN2] + u[data.iH2O]
        gas = (
            T     = u[data.iT],
            p     = u[data.ip],
            p_CO2 = u[data.ip] * u[data.iCO2] / c_total,
            p_H2O = u[data.ip] * u[data.iH2O] / c_total,
            p_N2  = u[data.ip] * u[data.iN2]  / c_total,
            c_CO2 = u[data.iCO2],
            c_H2O = u[data.iH2O],
            c_N2  = u[data.iN2],
        )

        q_star_H2O = sorb_params.q_star_H2O(gas, sorb_params.isotherm_params)
        q_star_CO2 = sorb_params.q_star_CO2(gas, q_star_H2O, sorb_params.isotherm_params)
        q_star_N2  = sorb_params.q_star_N2(gas, sorb_params.isotherm_params)

        y[data.iq_H2O] = - sorb_params.k_H2O * (q_star_H2O - u[data.iq_H2O])
        y[data.iq_CO2] = - _eval_k(sorb_params.k_CO2, u[data.iCO2], q_star_CO2) * (q_star_CO2 - u[data.iq_CO2])
        y[data.iq_N2]  = - _eval_k(sorb_params.k_N2,  u[data.iN2],  q_star_N2)  * (q_star_N2  - u[data.iq_N2])

        if NaN in u || NaN in y
            @show u[data.iCO2].value
            @show u[data.iN2].value
            @show _eval_k(sorb_params.k_CO2, u[data.iCO2], q_star_CO2).value
            @show _eval_k(sorb_params.k_N2,  u[data.iN2],  q_star_N2).value
            @show u[data.iq_CO2].value
            @show u[data.iq_N2].value
        end
    end
end

function bcondition(y, u, bnode, data)
    @inbounds begin
        params = data.step_params
        # inlet Neumann
        boundary_neumann!(y, u, bnode; species=data.iN2, region=data.Γ_in, value = params.u_feed * params.c_N2_feed)
        boundary_neumann!(y, u, bnode; species=data.iCO2, region=data.Γ_in, value = params.u_feed * params.c_CO2_feed)
        boundary_neumann!(y, u, bnode; species=data.iH2O, region=data.Γ_in, value = params.u_feed * params.c_H2O_feed)

        C_gas_feed = params.c_CO2_feed * data.phys_params.Cₚ_CO2 + params.c_H2O_feed * data.phys_params.Cₚ_H2O + params.c_N2_feed * data.phys_params.Cₚ_N2
        boundary_neumann!(y, u, bnode; species=data.iT, region=data.Γ_in, value = params.u_feed * C_gas_feed * params.T_feed)

        # Pressure Dirichlet
        if params.step_name == Adsorption
            boundary_dirichlet!(y, u, bnode; species=data.ip, region=data.Γ_out, value = params.P_out)
        elseif params.step_name == Pressurization
            P_out = params.P_out + (params.P_out_start - params.P_out) * exp(-params.λ * bnode.time)
            boundary_dirichlet!(y, u, bnode; species=data.ip, region=data.Γ_in, value = P_out)
        elseif params.step_name in (Blowdown, Preheating, Heating, Desorption)
            P_out = params.P_out + (params.P_out_start - params.P_out) * exp(-params.λ * bnode.time)
            boundary_dirichlet!(y, u, bnode; species=data.ip, region=data.Γ_out, value = P_out)
        elseif params.step_name == Evacuation
            P_out = params.P_out + (params.P_out_start - params.P_out) * exp(-params.λ * bnode.time)
            boundary_dirichlet!(y, u, bnode; species=data.ip, region=data.Γ_in, value = P_out)
        end
    end
end

function boutflow(y, u, edge, data)
    @inbounds begin
        params = data.step_params
        if params.step_name == Cooling
            return nothing
        end

        vh = data.velocity(u, edge, data)
        if params.step_name == Pressurization
            if vh < 0 return nothing end
            if outflownode(edge) == data.Γ_in
                y[data.iN2] = -vh * params.c_N2_feed
                y[data.iCO2] = -vh * params.c_CO2_feed
                y[data.iH2O] = -vh * params.c_H2O_feed

                C_gas_feed = params.c_CO2_feed * data.phys_params.Cₚ_CO2 + params.c_H2O_feed * data.phys_params.Cₚ_H2O + params.c_N2_feed * data.phys_params.Cₚ_N2
                y[data.iT] = -vh * C_gas_feed * params.T_feed
            end
        elseif params.step_name == Evacuation
            if vh > 0 return nothing end
            if outflownode(edge) == data.Γ_in
                nodeidx = outflownode(edge)
                y[data.iN2]  = -vh * u[data.iN2, nodeidx]
                y[data.iCO2] = -vh * u[data.iCO2, nodeidx]
                y[data.iH2O] = -vh * u[data.iH2O, nodeidx]

                C_gas = u[data.iCO2, nodeidx] * data.phys_params.Cₚ_CO2 + u[data.iH2O, nodeidx] * data.phys_params.Cₚ_H2O + u[data.iN2, nodeidx] * data.phys_params.Cₚ_N2
                y[data.iT] = -vh * C_gas * u[data.iT, nodeidx]
            end
        else
            if vh < 0 return nothing end
            if outflownode(edge) == data.Γ_out
                nodeidx = outflownode(edge)
                y[data.iN2]  = -vh * u[data.iN2, nodeidx]
                y[data.iCO2] = -vh * u[data.iCO2, nodeidx]
                y[data.iH2O] = -vh * u[data.iH2O, nodeidx]

                C_gas = u[data.iCO2, nodeidx] * data.phys_params.Cₚ_CO2 + u[data.iH2O, nodeidx] * data.phys_params.Cₚ_H2O + u[data.iN2, nodeidx] * data.phys_params.Cₚ_N2
                y[data.iT] = -vh * C_gas * u[data.iT, nodeidx]
            end
        end
    end
end

function source(f, node, data)
    return
end
