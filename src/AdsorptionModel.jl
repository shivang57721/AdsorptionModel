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
    export_to_hdf5,
    PHYS_PARAMS_DEFAULT,
    COST_PARAMS_DEFAULT,
    COL_PARAMS_ARVIND,
    COL_PARAMS_YOUNG,
    COL_PARAMS_LARGESCALE,
    COL_PARAMS_HAGHPANAH,
    COL_PARAMS_STAMPI,
    SORB_PARAMS_LEWATIT,
    SORB_PARAMS_NbOFFIVE,
    SORB_PARAMS_ZEOLITE,
    SORB_PARAMS_APDES

@enum StepType Adsorption Blowdown Preheating Heating Desorption Cooling Pressurization

using DifferentialEquations
using VoronoiFVM
using LinearAlgebra
using Logging
using Dictionaries
using JLD2
using HDF5
using ExtendableGrids
using ArrayInterface
using ForwardDiff
using SparseArrays
ArrayInterface.issingular(A::SparseMatrixCSC{<:ForwardDiff.Dual}) = false

include("params.jl")

# Isotherms
include("isotherms/helpers.jl")
include("isotherms/toth.jl")
include("isotherms/gab.jl")
include("isotherms/langmuir.jl")

# Sorbent presets
include("sorbents/lewatit.jl")
include("sorbents/nbofffive.jl")
include("sorbents/zeolite.jl")
include("sorbents/APDES.jl")

# Column / physical / cost presets
include("presets.jl")

# PDE solver kernels and system initialization
include("solver/kernels.jl")
include("solver/system.jl")

# Post-processing, plotting, export, simulation API
include("post_processing.jl")
include("plotting.jl")
include("hdf5_export.jl")
include("simulation.jl")

end # module
