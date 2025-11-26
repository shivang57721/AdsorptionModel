# Simulation of adsorption equations with VoronoiFVM.jl

module AdsorptionModel

export
    PhysicalParams,
    ColumnParams,
    SorbentParams,
    IsothermParams,
    CostParams,
    StepType, Adsorption, Preheating, Heating, Desorption, Cooling, Pressurization,
    OperatingParameters,
    initialize_system,
    run_simulation,
    get_cycle_params,
    simulate_process,
    is_feasible,
    simulate_DCB,
    post_process,
    plot_grid_idx,
    PHYS_PARAMS_DEFAULT,
    COST_PARAMS_DEFAULT,
    COL_PARAMS_ARVIND,
    COL_PARAMS_YOUNG,
    COL_PARAMS_LARGESCALE,
    SORB_PARAMS_LEWATIT,
    SORB_PARAMS_NbOFFIVE

@enum StepType Adsorption Preheating Heating Desorption Cooling Pressurization

using DifferentialEquations
using VoronoiFVM
using LinearAlgebra
using Logging
using Dictionaries
using JLD2
using ExtendableGrids
include("params.jl")
include("default_params.jl")
include("post_processing.jl")
include("plotting.jl")
include("main_interface.jl")

using ArrayInterface
using ForwardDiff
using SparseArrays
ArrayInterface.issingular(A::SparseMatrixCSC{<:ForwardDiff.Dual}) = false

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

function flux_exponential(f, u, edge, data)
    @inbounds begin
        ε_bed = data.sorb_params.ε_bed
        Δz = edgelength(edge)
        vh = data.velocity(u, edge, data)

        # totals
        c1_N2  = u[data.iN2, 1]
        c1_CO2 = u[data.iCO2, 1]
        c1_H2O = u[data.iH2O, 1]
        c2_N2  = u[data.iN2, 2]
        c2_CO2 = u[data.iCO2, 2]
        c2_H2O = u[data.iH2O, 2]

        c_total_1 = c1_N2 + c1_CO2 + c1_H2O
        c_total_2 = c2_N2 + c2_CO2 + c2_H2O
        c_total_edge = (c_total_1 + c_total_2) * 0.5

        # mole fractions
        yN2_1 = c1_N2 / c_total_1
        yCO2_1 = c1_CO2 / c_total_1
        yH2O_1 = c1_H2O / c_total_1

        yN2_2 = c2_N2 / c_total_2
        yCO2_2 = c2_CO2 / c_total_2
        yH2O_2 = c2_H2O / c_total_2

        # dispersion
        D_L = data.D_L(vh / Δz)
        bp, bm = ε_bed * D_L * c_total_edge .* fbernoulli_pm(vh / (ε_bed * D_L))

        f[data.iN2] = bm * yN2_1 - bp * yN2_2
        f[data.iCO2] = bm * yCO2_1 - bp * yCO2_2
        f[data.iH2O] = bm * yH2O_1 - bp * yH2O_2

        # Thermal flux
        Cgas1 = c1_CO2 * data.phys_params.Cₚ_CO2 + c1_H2O * data.phys_params.Cₚ_H2O + c1_N2 * data.phys_params.Cₚ_N2
        Cgas2 = c2_CO2 * data.phys_params.Cₚ_CO2 + c2_H2O * data.phys_params.Cₚ_H2O + c2_N2 * data.phys_params.Cₚ_N2
        Cgas_edge = (Cgas1 + Cgas2) * 0.5

        K_L = data.K_L(vh / Δz, Cgas_edge)
        bpT, bmT = (ε_bed * K_L) .* fbernoulli_pm(vh * Cgas_edge / (ε_bed * K_L))
        f[data.iT] = bmT * u[data.iT, 1] - bpT * u[data.iT, 2]
    end
end

function flux_exponential_test(f, u, edge, data)
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

function storage(y, u, node, data)
    @inbounds begin
        ε_total = data.sorb_params.ε_total
        C_solid = data.sorb_params.C_solid
        ρ_bed  = data.sorb_params.ρ_bed
        ΔH_CO2 = data.sorb_params.ΔH_CO2
        ΔH_H2O = data.sorb_params.ΔH_H2O
   
        N2 = u[data.iN2]
        CO2 = u[data.iCO2]
        H2O = u[data.iH2O]

        y[data.iN2]  = ε_total * N2
        y[data.iCO2] = ε_total * CO2
        y[data.iH2O] = ε_total * H2O

        C_gas = CO2 * data.phys_params.Cₚ_CO2 + H2O * data.phys_params.Cₚ_H2O + N2 * data.phys_params.Cₚ_N2
        C_ads = u[data.iq_CO2] * data.phys_params.Cₚ_CO2 + u[data.iq_H2O] * data.phys_params.Cₚ_H2O

        y[data.iT] = (ε_total * C_gas + ρ_bed * C_solid + ρ_bed * C_ads) * u[data.iT] - ε_total * u[data.ip]

        y[data.iT_wall] = u[data.iT_wall]

        # adsorbed species
        y[data.iq_CO2] = u[data.iq_CO2]
        y[data.iq_H2O] = u[data.iq_H2O]

        y[data.iCO2] += ρ_bed * u[data.iq_CO2]
        y[data.iH2O] += ρ_bed * u[data.iq_H2O]
        y[data.iT] += ρ_bed * ΔH_CO2 * u[data.iq_CO2] + ρ_bed * ΔH_H2O * u[data.iq_H2O]
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

        # Parital pressure of H2O and CO2 for calculating q_star
        p_H2O = u[data.ip] * u[data.iH2O] / (u[data.iCO2] + u[data.iN2] + u[data.iH2O])
        p_CO2 = u[data.ip] * u[data.iCO2] / (u[data.iCO2] + u[data.iN2] + u[data.iH2O])

        q_star_H2O = sorb_params.q_star_H2O(u[data.iT], p_H2O, sorb_params.isotherm_params)
        q_star_CO2 = sorb_params.q_star_CO2(u[data.iT], p_CO2, q_star_H2O, sorb_params.isotherm_params)

        k_CO2 = sorb_params.k_CO2
        k_H2O = sorb_params.k_H2O

        y[data.iq_H2O] = - k_H2O * (q_star_H2O - u[data.iq_H2O])
        y[data.iq_CO2] = - k_CO2 * (q_star_CO2 - u[data.iq_CO2])
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
            P_out = params.P_out + (params.P_out_start - params.P_out) * exp(-0.11 * bnode.time)
            boundary_dirichlet!(y, u, bnode; species=data.ip, region=data.Γ_in, value = P_out)
        elseif params.step_name in (Preheating, Heating, Desorption)
            P_out = params.P_out + (params.P_out_start - params.P_out) * exp(-0.11 * bnode.time)
            boundary_dirichlet!(y, u, bnode; species=data.ip, region=data.Γ_out, value = params.P_out)
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
        if vh < 0
            return nothing
        end

        if params.step_name == Pressurization
            if outflownode(edge) == data.Γ_in
                y[data.iN2] = -vh * params.c_N2_feed
                y[data.iCO2] = -vh * params.c_CO2_feed
                y[data.iH2O] = -vh * params.c_H2O_feed

                C_gas_feed = params.c_CO2_feed * data.phys_params.Cₚ_CO2 + params.c_H2O_feed * data.phys_params.Cₚ_H2O + params.c_N2_feed * data.phys_params.Cₚ_N2
                y[data.iT] = -vh * C_gas_feed * params.T_feed
            end
        else
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

struct AdsorptionData{T, DL, KL, VEL}
    ip::Int; iN2::Int; iCO2::Int; iH2O::Int; iT::Int; iT_wall::Int; iq_CO2::Int; iq_H2O::Int
    Γ_in::Int; Γ_out::Int
    species::Dict{String,Int}
    step_params::OperatingParameters
    col_params::ColumnParams
    sorb_params::SorbentParams
    phys_params::PhysicalParams
    darcy_k::T
    D_L::DL
    K_L::KL
    velocity::VEL
    ergun_b::T
    ergun_c::T
end

function AdsorptionData(; T::Type, step_params::OperatingParameters, col_params::ColumnParams,
                        sorb_params::SorbentParams, phys_params::PhysicalParams, velocity)
    # fully concrete closures (no captured mutable structs)
    D_L = let γ₁=phys_params.γ₁, γ₂=phys_params.γ₂, Dₘ=sorb_params.Dₘ, dₚ=sorb_params.dₚ, ε=sorb_params.ε_bed
        u -> γ₁*Dₘ + γ₂*dₚ*u/ε
    end
    K_L = let D_L=D_L
        (u, C_gas) -> D_L(u) * C_gas
    end
    darcy_k = T(150) * phys_params.μ * (1 - sorb_params.ε_bed)^2 / (sorb_params.ε_bed^3 * sorb_params.dₚ^2)

    # Precompute Ergun constants (use feed properties)
    ρ_gas = step_params.P_out / (phys_params.R * step_params.T_amb) *
            (step_params.y_CO2_feed * phys_params.MW_CO2 + step_params.y_H2O_feed * phys_params.MW_H2O + step_params.y_N2_feed * phys_params.MW_N2) * 1e-3
    ε = sorb_params.ε_bed; dₚ = sorb_params.dₚ; μ = phys_params.μ
    b = T(150) * μ * (1 - ε) / (dₚ * T(1.75) * ρ_gas)
    c = ε^3 * dₚ / (T(1.75) * ρ_gas * (1 - ε))

    return AdsorptionData{T, typeof(D_L), typeof(K_L), typeof(velocity)}(1,2,3,4,5,6,7,8,
        1,2,
        Dict("P"=>1,"N2"=>2,"CO2"=>3,"H2O"=>4,"T"=>5,"T_wall"=>6,"q_CO2"=>7,"q_H2O"=>8),
        step_params, col_params, sorb_params, phys_params,
        darcy_k, D_L, K_L, velocity, b, c)
end

# -----------------------------------------------------------------------------
# Initialization
# -----------------------------------------------------------------------------
function initialize_system(; T::Type=Float64, N::Int,
        op_params::OperatingParameters, col_params::ColumnParams, sorb_params::SorbentParams, phys_params::PhysicalParams,
        velocity=ergun_velocity, flux=flux_exponential, reaction=reaction, bcondition=bcondition, boutflow=boutflow, source=source)

    data = AdsorptionData(; T, step_params=deepcopy(op_params), col_params, sorb_params, phys_params, velocity)
    X = range(0, stop=col_params.L, length=N)
    grid = VoronoiFVM.Grid(X)

    sys = VoronoiFVM.System(grid;
        storage, flux, reaction, source, bcondition, boutflow,
        data, outflowboundaries=[data.Γ_in, data.Γ_out], species=values(data.species), matrixtype=:sparse,
        unknown_storage=:dense, assembly=:edgewise,
        valuetype=T
    )

    inival = unknowns(sys)
    inival[data.ip, :]      .= op_params.P_out
    inival[data.iT, :]      .= op_params.T_feed
    inival[data.iT_wall, :] .= op_params.T_amb
    inival[data.iN2, :]     .= op_params.c_total_feed
    inival[data.iCO2, :]    .= zero(T)
    inival[data.iH2O, :]    .= zero(T)
    inival[data.iq_CO2, :]  .= zero(T)
    inival[data.iq_H2O, :]  .= zero(T)

    return (; sys, data, inival, grid)
end

# Ok so basically this doesn't really work - I think best bet is to keep it as before and
# only `activate` the baseline when running AD?
# const BASELINE = 20 * 3600.0

# function VoronoiFVM.eval_rhs!(du, u, state, t)
#     T = state.data.step_params.step_name ∈ [Adsorption, Desorption, Heating] ? state.data.step_params.duration / BASELINE : 1.0
#     VoronoiFVM._eval_res_jac!(state, u, T * t)
#     du .= T .* -VoronoiFVM.dofs(state.residual)
#     state.history.nf += 1
#     return nothing
# end
end # module