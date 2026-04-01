using AdsorptionModel
using VoronoiFVM
using Plots
using GridVisualize
using DifferentialEquations
using LinearAlgebra
using ForwardDiff

AM = AdsorptionModel

phys_params = PHYS_PARAMS_DEFAULT()

col_params = ColumnParams(Float64;
    Rᵢ = 0.4,
    Rₒ = 0.41,
    L = 1,
    h_L = 14,
    h_W = 8
)

sorb_params = SorbentParams(;
                        ε_bed = 0.4,
                        ε_total = 0.54,
                        dₚ = 0.00052,
                        ρ_bed = 528.0,
                        ρ_particle = 880.0,
                        k_CO2 = 0.003,
                        k_H2O = 0.0086,
                        k_N2 = 0.0,
                        ΔH_CO2 = -70000.0,
                        ΔH_H2O = -46000.0,
                        ΔH_N2 = 0.0,
                        Dₘ = 1.3e-5,
                        C_solid = 1580.0,
                        isotherm_params = AdsorptionModel.ISOTHERM_PARAMS_LEWATIT(),
                        q_star_H2O = (T, p_H2O, params) -> 1.0,
                        q_star_CO2 = (T, p_CO2, q_H2O, params) -> 1.0
                    )

# Extract parameters
R = phys_params.R
γ₁ = phys_params.γ₁
γ₂ = phys_params.γ₂
μ = phys_params.μ
Cₚ_CO2 = phys_params.Cₚ_CO2
Cₚ_H2O = phys_params.Cₚ_H2O
Cₚ_N2 = phys_params.Cₚ_N2
Cₚ_wall = phys_params.Cₚ_wall
MW_CO2 = phys_params.MW_CO2
MW_N2 = phys_params.MW_N2

Rᵢ = col_params.Rᵢ
Rₒ = col_params.Rₒ
h_L = col_params.h_L
h_W = col_params.h_W
a_wall =  π * (Rₒ^2 - Rᵢ^2)

ε_bed = sorb_params.ε_bed
dₚ = sorb_params.dₚ
ε_total = sorb_params.ε_total
C_solid = sorb_params.C_solid
ρ_bed = sorb_params.ρ_bed
k_CO2 = sorb_params.k_CO2
k_H2O = sorb_params.k_H2O
ΔH_CO2 = sorb_params.ΔH_CO2
ΔH_H2O = sorb_params.ΔH_H2O

# The manufactured solutions
c_CO2(x,t) = sin(π * x) * exp(-t)
c_N2(x,t) = cos(π * x)^2 * exp(-t)
c_H2O(x,t) = 0.5 * sin(2π*x)^2 * exp(-0.5 * t)
q_H2O(x,t) = (x - x^2) * exp(-2t)
q_CO2(x,t) = 0.5 + x * exp(-t)
T(x,t) = 300 + sin(2π * x + 0.5) * exp(-t)
T_wall(x,t) = 300 + cos(2π * x - 0.5) * exp(-t)
p(x,t) = 1e5 - 1000.0 * x

# Dict of exact solutions
fexact = Dict("P"=> p, "N2"=>c_N2, "CO2"=>c_CO2, "H2O"=>c_H2O, "T"=>T, "T_wall"=>T_wall, "q_CO2"=>q_CO2, "q_H2O"=>q_H2O)

# Define velocity with Ergun equation
∂p∂t(x,t) = ForwardDiff.derivative(t -> p(x,t), t)
∂p∂x(x,t) = ForwardDiff.derivative(x -> p(x,t), x)
ρ_gas = p(1,0) / (phys_params.R * T(0,0)) * (1.0 * phys_params.MW_N2) * 1e-3
b = 150 * μ * (1 - ε_bed) / (dₚ * 1.75 * ρ_gas)
c = ε_bed^3 * dₚ / (1.75 * ρ_gas * (1 - ε_bed))
u(x,t) = -2c * ∂p∂x(x,t) / (b + sqrt(b^2 + 4c * sign(∂p∂x(x,t)) * ∂p∂x(x,t)))

# Define operating parameters for adsorption
op_params = OperatingParameters(Float64;
                step_name = Adsorption,
                u_feed = u(0,0),
                T_feed = T(0,0),
                T_amb = T(0,0),
                P_out = p(1,0),
                y_CO2_feed = 0.0,
                y_H2O_feed = 0.0,
                duration = 1)

# Derivatives w.r.t t
∂c_CO2∂t(x,t) = ForwardDiff.derivative(t -> c_CO2(x,t), t)
∂c_N2∂t(x,t)  = ForwardDiff.derivative(t -> c_N2(x,t), t)
∂c_H2O∂t(x,t) = ForwardDiff.derivative(t -> c_H2O(x,t), t)
∂q_CO2∂t(x,t) = ForwardDiff.derivative(t -> q_CO2(x,t), t)
∂q_H2O∂t(x,t) = ForwardDiff.derivative(t -> q_H2O(x,t), t)
∂T∂t(x,t)     = ForwardDiff.derivative(t -> T(x,t), t)
∂T_wall∂t(x,t)= ForwardDiff.derivative(t -> T_wall(x,t), t)

# Derivatives w.r.t x
∂c_CO2∂x(x,t) = ForwardDiff.derivative(x -> c_CO2(x,t), x)
∂c_N2∂x(x,t)  = ForwardDiff.derivative(x -> c_N2(x,t), x)
∂c_H2O∂x(x,t) = ForwardDiff.derivative(x -> c_H2O(x,t), x)
∂T∂x(x,t)     = ForwardDiff.derivative(x -> T(x,t), x)

# Second derivatives w.r.t x
∂²c_CO2∂x²(x,t) = ForwardDiff.derivative(x -> ∂c_CO2∂x(x,t), x)
∂²c_N2∂x²(x,t)  = ForwardDiff.derivative(x -> ∂c_N2∂x(x,t), x)
∂²c_H2O∂x²(x,t) = ForwardDiff.derivative(x -> ∂c_H2O∂x(x,t), x)
∂²T∂x²(x,t)     = ForwardDiff.derivative(x -> ∂T∂x(x,t), x)

∂uc_CO2∂x(x,t) = ForwardDiff.derivative(x -> u(x,t) * c_CO2(x,t), x)
∂uc_N2∂x(x,t) = ForwardDiff.derivative(x -> u(x,t) * c_N2(x,t), x)
∂uc_H2O∂x(x,t) = ForwardDiff.derivative(x -> u(x,t) * c_H2O(x,t), x)
∂uT∂x(x,t) = ForwardDiff.derivative(x -> u(x,t) * T(x,t), x)

function source_MMStest(f, node, data)
    x = node[1]
    t = node.time
    D_L = data.D_L(u(x,t))

    f[data.ip] = p(x,t) - (c_CO2(x,t) + c_H2O(x,t) + c_N2(x,t)) * R * T(x, t)
    
    f[data.iCO2] = ε_total * ∂c_CO2∂t(x,t) + ∂uc_CO2∂x(x,t) + ρ_bed * ∂q_CO2∂t(x,t) - ε_bed * D_L * ∂²c_CO2∂x²(x,t)
    f[data.iN2] = ε_total * ∂c_N2∂t(x,t) + ∂uc_N2∂x(x,t) - ε_bed * D_L * ∂²c_N2∂x²(x,t)
    f[data.iH2O] = ε_total * ∂c_H2O∂t(x,t) + ∂uc_H2O∂x(x,t) + ρ_bed * ∂q_H2O∂t(x,t) - ε_bed * D_L * ∂²c_H2O∂x²(x,t)

    f[data.iq_H2O] = ∂q_H2O∂t(x,t) - k_H2O * (1 - q_H2O(x,t))
    f[data.iq_CO2] = ∂q_CO2∂t(x,t) - k_CO2 * (1 - q_CO2(x,t))

    C_gas = c_CO2(x,t) * Cₚ_CO2 + c_H2O(x,t) * Cₚ_H2O + c_N2(x,t) * Cₚ_N2
    C_ads = q_CO2(x,t) * Cₚ_CO2 + q_H2O(x,t) * Cₚ_H2O
    K_L = data.K_L(u(x,t), C_gas)
    f[data.iT] = (ε_total * C_gas + ρ_bed * C_solid + ρ_bed * C_ads) * ∂T∂t(x,t) - ε_total * ∂p∂t(x,t) +
                    C_gas * ∂uT∂x(x,t) + ρ_bed * ΔH_CO2 * ∂q_CO2∂t(x,t) + ρ_bed * ΔH_H2O * ∂q_H2O∂t(x,t) +
                    2h_L / Rᵢ * (T(x,t) - T_wall(x,t)) - ε_bed * K_L * ∂²T∂x²(x,t)
    
    f[data.iT_wall] = ∂T_wall∂t(x,t) - 2π/(Cₚ_wall * a_wall) * (h_L * Rᵢ * (T(x,t) - T_wall(x,t)) - h_W * Rₒ * (T_wall(x,t) - data.step_params.T_amb))
end

function bcondition_MMStest(y, u, bnode, data)
    t = bnode.time
    
    boundary_dirichlet!(y, u, bnode; species=data.iCO2, region=data.Γ_in, value = c_CO2(0,t))
    boundary_dirichlet!(y, u, bnode; species=data.iH2O, region=data.Γ_in, value = c_H2O(0,t))
    boundary_dirichlet!(y, u, bnode; species=data.iN2, region=data.Γ_in, value = c_N2(0,t))
    boundary_dirichlet!(y, u, bnode; species=data.iT, region=data.Γ_in, value = T(0,t))
    
    boundary_dirichlet!(y, u, bnode; species=data.iCO2, region=data.Γ_out, value = c_CO2(1,t))
    boundary_dirichlet!(y, u, bnode; species=data.iH2O, region=data.Γ_out, value = c_H2O(1,t))
    boundary_dirichlet!(y, u, bnode; species=data.iN2, region=data.Γ_out, value = c_N2(1,t))
    boundary_dirichlet!(y, u, bnode; species=data.iT, region=data.Γ_out, value = T(1,t))

    boundary_dirichlet!(y, u, bnode; species=data.ip, region=data.Γ_out, value = p(1,t))
end

function boutflow_MMStest(y, u, edge, data)
    return
end

function system(; op_params, N=10, kwargs...)
    sys, data, inival, grid = initialize_system(;N, op_params, col_params, sorb_params, phys_params,
                                            source=source_MMStest, bcondition=bcondition_MMStest, boutflow=boutflow_MMStest)
    
    inival[data.ip, :]      .= map(x -> p(x,0), grid)
    inival[data.iT, :]      .= map(x -> T(x,0), grid)
    inival[data.iT_wall, :] .= map(x -> T_wall(x,0), grid)
    inival[data.iN2, :]     .= map(x -> c_N2(x,0), grid)
    inival[data.iCO2, :]    .= map(x -> c_CO2(x,0), grid)
    inival[data.iH2O, :]    .= map(x -> c_H2O(x,0), grid)
    inival[data.iq_CO2, :]  .= map(x -> q_CO2(x,0), grid)
    inival[data.iq_H2O, :]  .= map(x -> q_H2O(x,0), grid)

    problem = ODEProblem(sys,inival,(0,op_params.duration))
    odesol = solve(problem, FBDF(); kwargs...)
    sol = reshape(odesol,sys)

    return sys, data, grid, sol
end


function plot_approximation(op_params;species="CO2", kwargs...)
    sys, data, grid, sol = system(;op_params, kwargs...)
    sol_exact = map(x -> fexact[species](x, sol.t[end]), grid)

    vis = GridVisualizer(legend = :rt, Plotter=Plots)
    scalarplot!(vis, sys, sol[end]; species=data.species[species], color = :green, markershape = :circle, markersize = 3, markevery = 1, label = "numerical")
    scalarplot!(vis, grid, sol_exact, color = :red, clear = false, markershape = :none, label = "exact")
    xlabel!("z")
    ylabel!("$species(z,T)")
    # title!("Numerical approximation of $species at final time T")
    reveal(vis)
end

# Convergence test space
function convtest_space(op_params)
    Ns = Int.([10 * 2.0^i for i in 2:8])
    L2 = Float64[]
    H1 = Float64[]

    for N in Ns
        @show N
        sys, data, grid, sol = system(;N, op_params, reltol=1e-10, abstol=1e-12)
        sol_exact = unknowns(sys)
        for (species, idx) in data.species
            sol_exact[idx, :] .= map(x -> fexact[species](x, sol.t[end]), grid)
        end
        push!(H1, h1norm(sys, sol[end] .- sol_exact))
        push!(L2, l2norm(sys, sol[end] .- sol_exact))
    end

    vis = GridVisualizer(
        xscale = :log, yscale = :log, legend = :lt, Plotter=Plots,
    )
    H = 1 ./ Ns
    # scalarplot!(vis, H, H1, color = :green, label = "H1")
    scalarplot!(vis, H, L2, color = :green, label = "L2", clear = false, linestyle = :dash)

    scalarplot!(vis, H, H * 20, color = :red, linestyle = :solid, label = "O(h)", clear = false)
    scalarplot!(vis, H, H .^ 2 * 800, color = :red, linestyle = :dash, label = "O(h^2)", clear = false)

    xlabel!("Mesh size h")
    ylabel!("L² error at final time T")
    # title!("Spatial Convergence (MMS)")
    reveal(vis)
end
    
# Convergence test time
function convtest_time(op_params)
    N = 2000
    dts = [0.1 * 2.0^(-i) for i in 4:10]

    L2 = Float64[]
    for dt in dts
        @show dt
        sys, data, grid, sol = system(;N=Int64(N), op_params, adaptive=false, dt)
        sol_exact = unknowns(sys)
        for (species, idx) in data.species
            sol_exact[idx, :] .= map(x -> fexact[species](x, sol.t[end]), grid)
        end
        push!(L2, l2norm(sys, sol[end] .- sol_exact))
    end

    vis = GridVisualizer(
        xscale = :log, yscale = :log, legend = :lt, Plotter=Plots,
    )
    H = dts
    scalarplot!(vis, H, L2, color = :green, label = "L2", clear = false, linestyle = :dash)

    scalarplot!(vis, H, H * 1000, color = :red, linestyle = :solid, label = "O(h)", clear = false)
    scalarplot!(vis, H, H .^ 2 * 10000, color = :red, linestyle = :dash, label = "O(h^2)", clear = false)

    xlabel!("Time step Δt")
    ylabel!("L² error at final time T")
    # title!("Temporal Convergence (MMS)")
    reveal(vis)
end

convtest_space(op_params)
savefig("test/convergence_test_h.png")

convtest_time(op_params)
savefig("test/convergence_test_dt.png")

plot_approximation(op_params; species="CO2", N=100, reltol=1e-9)
savefig("test/CO2_MMS_test_numerical_approximation.png")