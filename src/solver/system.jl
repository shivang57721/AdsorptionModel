# AdsorptionData struct and system initialization

struct AdsorptionData{T, DL, KL, VEL, OP, CP, SP, PP}
    ip::Int; iN2::Int; iCO2::Int; iH2O::Int; iT::Int; iT_wall::Int; iq_CO2::Int; iq_H2O::Int; iq_N2::Int
    Γ_in::Int; Γ_out::Int
    species::Dict{String,Int}
    step_params::OP
    col_params::CP
    sorb_params::SP
    phys_params::PP
    darcy_k::T
    D_L::DL
    K_L::KL
    velocity::VEL
    ergun_b::T
    ergun_c::T
end

function AdsorptionData(; step_params::OperatingParameters, col_params::ColumnParams,
                        sorb_params::SorbentParams, phys_params::PhysicalParams, velocity)
    D_L = let γ₁=phys_params.γ₁, γ₂=phys_params.γ₂, Dₘ=sorb_params.Dₘ, dₚ=sorb_params.dₚ, ε=sorb_params.ε_bed
        u -> γ₁*Dₘ + γ₂*dₚ*u/ε
    end
    K_L = let D_L=D_L
        (u, C_gas) -> D_L(u) * C_gas
    end
    darcy_k = 150 * phys_params.μ * (1 - sorb_params.ε_bed)^2 / (sorb_params.ε_bed^2 * sorb_params.dₚ^2)

    ρ_gas = step_params.P_out / (phys_params.R * step_params.T_amb) *
            (step_params.y_CO2_feed * phys_params.MW_CO2 + step_params.y_H2O_feed * phys_params.MW_H2O + step_params.y_N2_feed * phys_params.MW_N2) * 1e-3
    ε = sorb_params.ε_bed; dₚ = sorb_params.dₚ; μ = phys_params.μ
    b = 150 * μ * (1 - ε) / (dₚ * 1.75 * ρ_gas)
    c = ε^3 * dₚ / (1.75 * ρ_gas * (1 - ε))

    return AdsorptionData{typeof(b), typeof(D_L), typeof(K_L), typeof(velocity),
                          typeof(step_params), typeof(col_params), typeof(sorb_params), typeof(phys_params)}(
        1,2,3,4,5,6,7,8,9,
        1,2,
        Dict("P"=>1,"N2"=>2,"CO2"=>3,"H2O"=>4,"T"=>5,"T_wall"=>6,"q_CO2"=>7,"q_H2O"=>8,"q_N2"=>9),
        step_params, col_params, sorb_params, phys_params,
        darcy_k, D_L, K_L, velocity, b, c)
end

function initialize_system(; T::Type=Float64, N::Int,
        op_params::OperatingParameters, col_params::ColumnParams, sorb_params::SorbentParams, phys_params::PhysicalParams,
        velocity=darcy_velocity, flux=flux_exponential, reaction=reaction, bcondition=bcondition, boutflow=boutflow, source=source)

    data = AdsorptionData(; step_params=deepcopy(op_params), col_params, sorb_params, phys_params, velocity)
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
    inival[data.iq_N2,  :]  .= zero(T)

    return (; sys, data, inival, grid)
end
