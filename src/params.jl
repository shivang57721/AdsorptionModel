struct PhysicalParams{T<:Real}
    R::T
    γ₁::T
    γ₂::T
    μ::T
    Cₚ_CO2::T
    Cₚ_H2O::T
    Cₚ_N2::T
    Cₚ_wall::T
    λ_ev::T
    Cₚ_steam::T
    MW_CO2::T
    MW_N2::T
    MW_H2O::T
end

PhysicalParams(::Type{T}=Float64; R=8.314, γ₁=0.7, γ₂=0.5, μ=1.789e-5,
               Cₚ_CO2=37.520, Cₚ_H2O=34.2, Cₚ_N2=29.171, Cₚ_wall=4.0e6,
               λ_ev=2350.0, Cₚ_steam=2.08, MW_CO2=44.0, MW_N2=28.0, MW_H2O=18.0) where {T<:Real} =
    PhysicalParams{T}(T(R), T(γ₁), T(γ₂), T(μ), T(Cₚ_CO2), T(Cₚ_H2O), T(Cₚ_N2),
                      T(Cₚ_wall), T(λ_ev), T(Cₚ_steam), T(MW_CO2), T(MW_N2), T(MW_H2O))

struct IsothermParams{T<:Real}
    q∞::T; b₀::T; T_ref::T; ΔH₀::T; τ₀::T; α::T; # dry params
    b₀_wet::T; ΔH₀_wet::T; τ₀_wet::T; α_wet::T; q∞_wet::T; A::T; # wet params
    qₘ::T; C::T; D::T; F::T; G::T;
    γ::T; β::T
end

IsothermParams(::Type{T}=Float64; q∞=NaN, b₀=NaN, T_ref=NaN, ΔH₀=NaN, τ₀=NaN, α=NaN,
                b₀_wet=NaN, ΔH₀_wet=NaN, τ₀_wet=NaN, α_wet=NaN, q∞_wet=NaN, A=NaN,
                qₘ=NaN, C=NaN, D=NaN, F=NaN, G=NaN,
                γ=NaN, β=NaN) where {T<:Real} =
    IsothermParams{T}(T(q∞), T(b₀), T(T_ref), T(ΔH₀), T(τ₀), T(α),
                        T(b₀_wet), T(ΔH₀_wet), T(τ₀_wet), T(α_wet), T(q∞_wet), T(A),
                      T(qₘ), T(C), T(D), T(F), T(G),
                      T(γ), T(β))

struct ColumnParams{T<:Real}
    Rᵢ::T; Rₒ::T; L::T; h_L::T; h_W::T
end

ColumnParams(::Type{T}=Float64; Rᵢ, Rₒ, L, h_L, h_W) where {T<:Real} =
    ColumnParams{T}(T(Rᵢ), T(Rₒ), T(L), T(h_L), T(h_W))

struct SorbentParams{T<:Real}
    ε_bed::T; ε_total::T; dₚ::T; ρ_bed::T; ρ_particle::T
    k_CO2::T; k_H2O::T; k_N2::T
    ΔH_CO2::T; ΔH_H2O::T; ΔH_N2::T
    Dₘ::T; C_solid::T
    isotherm_params::IsothermParams
    q_star_CO2::Function
    q_star_H2O::Function
end

SorbentParams(::Type{T}=Float64; ε_bed, ε_total, dₚ, ρ_bed, ρ_particle,
              k_CO2, k_H2O, k_N2, ΔH_CO2, ΔH_H2O, ΔH_N2, Dₘ, C_solid,
              q_star_CO2::Function, q_star_H2O::Function, isotherm_params::IsothermParams) where {T<:Real} =
    SorbentParams{T}(T(ε_bed), T(ε_total), T(dₚ), T(ρ_bed), T(ρ_particle),
                     T(k_CO2), T(k_H2O), T(k_N2),
                     T(ΔH_CO2), T(ΔH_H2O), T(ΔH_N2),
                     T(Dₘ), T(C_solid), 
                     isotherm_params, q_star_CO2, q_star_H2O)

struct CostParams{T<:Real}
    efficiency_blower::T
    adiabatic_index::T
    eta_VP::T
    annual_operating_hours::T
    target_CO2_ton_per_year::T
end

CostParams(::Type{T}=Float64; efficiency_blower=0.6, adiabatic_index=1.4, eta_VP=0.7,
           annual_operating_hours=8000, target_CO2_ton_per_year=1000) where {T<:Real} =
    CostParams{T}(T(efficiency_blower), T(adiabatic_index), T(eta_VP),
                  T(annual_operating_hours), T(target_CO2_ton_per_year))

mutable struct OperatingParameters{T<:Real}
    u_feed::T; T_feed::T; T_amb::T; T_start::T
    y_CO2_feed::T; y_H2O_feed::T
    duration::T; step_name::StepType
    P_out::T; P_out_start::T
    y_N2_feed::T; c_total_feed::T
    c_N2_feed::T; c_CO2_feed::T; c_H2O_feed::T
    ΔT_preheating::T; T_safe_cooling::T
    q_CO2_saturation_limit::T; extra_heating_ratio::T
end

OperatingParameters(::Type{T}; u_feed=0, T_feed=293, y_CO2_feed=0, y_H2O_feed=0,
                    T_amb, P_out, duration, step_name::StepType, 
                    ΔT_preheating=NaN, T_safe_cooling=NaN, q_CO2_saturation_limit=NaN, extra_heating_ratio=NaN) where {T<:Real} = begin
    y_N2_feed = T(1) - T(y_CO2_feed) - T(y_H2O_feed)
    c_total_feed = T(P_out) / (T(8.314) * T(T_feed))
    c_N2_feed = y_N2_feed * c_total_feed
    c_CO2_feed = T(y_CO2_feed) * c_total_feed
    c_H2O_feed = T(y_H2O_feed) * c_total_feed
    OperatingParameters{T}(
        T(u_feed), T(T_feed), T(T_amb), T(T_amb),
        T(y_CO2_feed), T(y_H2O_feed),
        T(duration), step_name,
        T(P_out), T(P_out),
        y_N2_feed, c_total_feed,
        c_N2_feed, c_CO2_feed, c_H2O_feed,
        T(ΔT_preheating), T(T_safe_cooling),
        T(q_CO2_saturation_limit), T(extra_heating_ratio)
    )
end

function copy_params!(dest::OperatingParameters, src::OperatingParameters)
    for name in fieldnames(OperatingParameters)
        setfield!(dest, name, getfield(src, name))
    end
    return nothing
end