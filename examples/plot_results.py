"""
Plotting simulation results from AdsorptionModel HDF5 output.

Requires: h5py, numpy, matplotlib
    pip install h5py numpy matplotlib
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt

HDF5_FILE = "examples/simulation_result.h5"

STEP_ORDER = ["Adsorption", "Blowdown", "Preheating", "Heating", "Desorption", "Cooling", "Pressurization"]

# NOTE: Julia writes arrays in column-major (Fortran) order.
# h5py reads in row-major (C) order, so a Julia array of shape
# [n_species, n_nodes, n_times] is read in Python as [n_times, n_nodes, n_species].
# All indexing below accounts for this transposition.

# ------------------------------------------------------------
# Plot CO2 at column exit for the last cycle (all steps)
# ------------------------------------------------------------

with h5py.File(HDF5_FILE, "r") as f:
    # Species indices are 1-based (Julia convention); convert to 0-based for Python
    iCO2 = f["species"].attrs["CO2"] - 1

    last_cycle = max(int(k.replace("cycle_", "")) for k in f["cycles"].keys())

    all_times = []
    co2_exit  = []
    step_boundaries = []
    time_offset = 0.0

    for step_name in STEP_ORDER:
        path = f"cycles/cycle_{last_cycle}/{step_name}"
        t    = f[path]["t"][:]                 # [n_times]
        data = f[path]["data"][:]              # [n_times, n_nodes, n_species]

        all_times.extend(t + time_offset)
        co2_exit.extend(data[:, -1, iCO2])    # last node, CO2 species
        time_offset += t[-1]
        step_boundaries.append(time_offset)

    step_boundaries.pop()  # remove boundary after last step

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(all_times, co2_exit, label="Last cycle")
for t_boundary in step_boundaries:
    ax.axvline(t_boundary, linestyle="--", color="black", alpha=0.5)
ax.set_xlabel("Time (s)")
ax.set_ylabel("CO₂ concentration [mol/m³]")
ax.set_title("CO₂ at Column Exit — Last Cycle")
ax.legend()
plt.tight_layout()
plt.show()
plt.savefig("co2_at_exit.png")

# ------------------------------------------------------------
# Plot CO2 spatial profile at end of Adsorption (last cycle)
# ------------------------------------------------------------

with h5py.File(HDF5_FILE, "r") as f:
    iCO2 = f["species"].attrs["CO2"] - 1

    last_cycle = max(int(k.replace("cycle_", "")) for k in f["cycles"].keys())

    x    = f["grid/x_coords"][:]                                           # [n_nodes]
    data = f[f"cycles/cycle_{last_cycle}/Adsorption/data"][:]              # [n_times, n_nodes, n_species]

    co2_profile = data[-1, :, iCO2]   # last timestep, all nodes, CO2

fig, ax = plt.subplots()
ax.plot(x, co2_profile)
ax.set_xlabel("Column position [m]")
ax.set_ylabel("CO₂ concentration [mol/m³]")
ax.set_title("CO₂ Spatial Profile at End of Adsorption")
plt.tight_layout()
plt.show()
plt.savefig("x_co2_end_of_ads.png")

# ------------------------------------------------------------
# Read output KPIs
# ------------------------------------------------------------

with h5py.File(HDF5_FILE, "r") as f:
    output = dict(f["output"].attrs)

print("Output:")
for k, v in output.items():
    print(f"  {k}: {v}")
