---
Group: D
Members: Julius Kramer, Tim Scholl, Johannes Hupe
PR: https://github.com/JoHaHu/MolSim/pull/71
---

# Report Group-D

#### Pull Request: https://github.com/JoHaHu/MolSim/pull/71

## Tasks - Worksheet 5

### Media Link

In this worksheet, we once again have a couple of small and larger simulation video files (some of which are also in the repo) that we would also like to share via this OneDrive link for you to access them:
https://1drv.ms/f/s!Ar3cyyDhVL_oktssXaQhihLO8r66rw?e=dBheet

## Task 1 - Simulation of a membrane

- Extend Molecule Class: Modified to store neighboring molecules and initialized on a rectangular grid, setting direct
  and diagonal neighbors.
- Harmonic Potential: Implemented interactions for direct and diagonal neighbors using harmonic potential formulas.
- Self-Penetration Prevention: Applied truncated Lennard-Jones potential to ensure only repulsive forces are active.
- Ran the simulation with specified parameters.

## Task 2 - Parallelization

- Parallelized the most time-consuming parts of the algorithm using OpenMP.
- Implemented two parallelization strategies selectable via the input file.
- Ensured compatibility without OpenMP using precompiler statements (#ifdef OPENMP).
- Measured speedup with varying amount of threads.
- we use domain coloring to divide the domain in partition that can be iterated without locking

## Task 3 - Rayleigh-Taylor instability in 3D

- Extended the 2D Rayleigh-Taylor instability simulation to a 3D scenario.
- Most of the necessary features were already implemented or prepared for the extension to 3D
- Applied periodic boundaries on the additional third dimension. The upper and lower boundary is still reflecting so particles can interchange and mix and form the typical patterns.
- Conducted the experiment with specified parameters but with shorter end times for better runtimes.
- The long run with an end time of 100 took about 8 hours and 40 minutes on one of our laptops and would approximately have taken 4-5 hours on the 16-core Ryzen 9 (unfortunately the run was interrupted).
- The simulation shows the expected behavior that we saw in last worksheet's Rayleigh Taylor instability but in a large cuboid/die.
- After some time the denser fluid is settled below while the other one floats above.

### Contest 2

- Runtime on the cluster for the 1000 iterations of RT in 3D: **02:03:821** (mm:ss:ms)
- Runtime on the local machine for the 1000 iterations of RT in 3D: **00:43:933** (mm:ss:ms)

## Task 4 - Nano-scale flow simulation

- Fixed Wall Particles: Enabled fixing positions of particles in the outer cuboids, ensuring they do not move but still exert forces on fluid particles. This works by adding an additional attribute for cuboids
- In the position updating, only particles that are not fixed are considered and their positions updated. Therefore, the forces between wall and fluid particles still remain.
- Thermostat Extension: Extended the thermostat to ignore total fluid velocity, calculating temperature using velocity deviations from the mean.
- Velocity Adjustment: For each particle, subtracted the average velocity, scaled the deviation, and then added back the average velocity to obtain new velocities.
- Profile Computation: Implemented a component to compute density and velocity profiles along the x-axis by subdividing it into bins and calculating averages per bin.
- Iterations and amount of bins can be input in the XML with an optional tag to use only when needed. The corresponding csv files are automatically put out in the build folder.

=====================================================

### We studied various influences on simulation profiles:

**Gravity factor:**
- The higher the gravity factor, the denser the stable formation (with reflecting boundaries).
- The higher the gravity factor, the faster the flow (imaginary flow of a nanotube sequence) with periodic boundaries

**Mass:**
- Higher mass of wall molecules leads to the fluid sticking better to it.
- Vice-versa, lower mass of wall molecules or higher mass of fluid molecules leads to less interaction effects with the
  wall molecules.

**σ or ϵ of the molecules:**
- Higher σ expands the effective diameter of the particles.
- Higher ϵ leads to higher inter-particle attraction which increases the density of the liquid.

**Removing walls:**
- Particles flowing out in various directions

**Adding walls:**
- All reflecting boundaries neatly shows the "sagging" of the particles.
- Settling into stable formation.

**Old thermostat:**
- Maintains uniform temperature by scaling velocities directly, allowing potential fluid drift.

**New thermostat:**
- Accurately controls temperature by removing average fluid velocity before scaling, proventing net drift.

**No thermostat:**
- Temperature fluctuates significantly without control, and velocity distribution deviates over time. Potential for energy build-up and instability.

**Task 4 extra thermostat:**
- Maintains temperature by scaling only the x- and z-components of velocities and adjusts the y-component for mean flow,
  resulting in realistic flow dynamics without affecting the temperature.
  This approach differs by focusing on directional components and excluding flow speed from temperature calculations.
