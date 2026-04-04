# Simulation of adsorption equations with VoronoiFVM.jl

module AdsorptionModel

export
    # Parameter types
    PhysicalParams,
    ColumnParams,
    SorbentParams,
    IsothermParams,
    CostParams,
    OperatingParameters,

    # Process types
    ProcessType, TVSA, STVSA, PSA, TSA,

    # Step types
    StepType, Adsorption, Blowdown, Preheating, Heating, Desorption, Cooling, Pressurization, Evacuation,

    # Step duration types
    StepDuration, FixedDuration, SaturationLimit, HeatingUntilTarget, max_duration,

    # Step configuration
    StepConfig, step_sequence,

    # Simulation API
    simulate_process,
    simulate_DCB,
    is_feasible,
    post_process,
    plot_grid_idx,
    export_to_hdf5,

    # Presets
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

@enum StepType Adsorption Blowdown Preheating Heating Desorption Cooling Pressurization Evacuation

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

# Post-processing, plotting, export
include("post_processing.jl")
include("plotting.jl")
include("hdf5_export.jl")

# Process definitions and simulation API
include("processes.jl")
include("simulation.jl")

end # module
