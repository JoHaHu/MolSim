---
Group: D
Members: Julius Kramer, Tim Scholl, Johannes Hupe
PR: https://github.com/JoHaHu/MolSim/pull/71
---

# Report Group-D

## Tasks - Worksheet 5

### Task 1 - Simulation of a membrane

- Extend Molecule Class: Modified to store neighboring molecules and initialized on a rectangular grid, setting direct and diagonal neighbors.
- Harmonic Potential: Implemented interactions for direct and diagonal neighbors using harmonic potential formulas.
- Self-Penetration Prevention: Applied truncated Lennard-Jones potential to ensure only repulsive forces are active.
- Ran the simulation with specified parameters.

### Task 2 - Parallelization

- Parallelized the most time-consuming parts of the algorithm using OpenMP.
- Implemented two parallelization strategies selectable via the input file.
- Ensured compatibility without OpenMP using precompiler statements (#ifdef OPENMP).
- Measured speedup with varying amount of threads.

### Task 3 - Rayleigh-Taylor instability in 3D

- Extended the 2D Rayleigh-Taylor instability simulation to a 3D scenario.
- Applied periodic boundaries on the additional third dimension.
- Conducted the experiment with specified parameters.

### Contest 2

- 

### Task 4 - Nano-scale flow simulation

- Fixed Wall Particles: Enabled fixing positions of particles in the outer cuboids, ensuring they do not move but still exert forces on fluid particles.
- Thermostat Extension: Extended the thermostat to ignore total fluid velocity, calculating temperature using velocity deviations from the mean.
- Velocity Adjustment: For each particle, subtracted the average velocity, scaled the deviation, and then added back the average velocity to obtain new velocities.
- Profile Computation: Implemented a component to compute density and velocity profiles along the x-axis by subdividing it into bins and calculating averages per bin.

We studied various influences on simulation profiles:
=====================================================

Gravity factor:
- The higher the gravity factor, the denser the stable formation (with reflecting boundaries).

Mass:
- Higher mass of wall molecules leads to the fluid sticking better to it.
- Vice-versa, lower mass of wall molecules or higher mass of fluid molecules leads to less interaction effects with the wall molecules.

σ or ϵ of the molecules:
- Higher σ expands the effective diameter of the particles.
- Higher ϵ leads to higher inter-particle attraction which increases the density of the liquid.

Removing walls:
- Particles flowing out in various directions

Adding walls:
- All reflecting boundaries neatly shows the "sagging" of the particles.
- Settling into stable formation.

Old thermostat:
- Maintains uniform temperature by scaling velocities directly, allowing potential fluid drift.

New thermostat:
- Accurately controls temperature by removing average fluid velocity before scaling, proventing net drift.

No thermostat:
- Temperature fluctuates significantly without control, and velocity distribution deviates over time. Potential for energy build-up and instability.

Task 4 extra thermostat:
- Maintains temperature by scaling only the x- and z-components of velocities and adjusts the y-component for mean flow,
resulting in realistic flow dynamics without affecting the temperature.
This approach differs by focusing on directional components and excluding flow speed from temperature calculations.
