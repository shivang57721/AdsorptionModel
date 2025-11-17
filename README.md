*A Julia package for simulating sTVSA, TVSA, and TCSA adsorption processes for Direct Air Capture (DAC).*

# Overview

AdsorptionModel.jl provides a fast, flexible, and robust numerical framework for simulating 1D adsorption columns used in DAC processes, including:

- sTVSA — steam-temperature–vacuum swing adsorption

- TVSA — temperature–vacuum swing adsorption

- TCSA — temperature–concentration swing adsorption

The model uses:

- Finite Volume Method (FVM) for spatial discretization with VoronoiFVM.jl

- DifferentialEquations.jl as the DAE solver backend

- LDF (Linear Driving Force) model for interphase mass transfer

- Pressure drop  described by Ergun's equation

- Modular structure for sorbent, cycle, and process parameters

This package is designed for high-performance simulation and optimization of DAC processes, with an emphasis on speed, numerical stability, and extensibility.

# Installation

To install the package directly from GitHub, open the julia REPL in your project directory and run
```
] add https://github.com/shivang57721/AdsorptionModel.git
```

# Examples
Example scripts are provided in the `examples/` directory.

# Testing
TODO: To run the package test suite
```
] test AdsorptionModel
```
