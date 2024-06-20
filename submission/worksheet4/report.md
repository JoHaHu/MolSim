---
Group: D
Members: Julius Kramer, Tim Scholl, Johannes Hupe
PR: https://github.com/JoHaHu/MolSim/pull/50
---

# Report Group-D

## Tasks - Worksheet 4

### Task 1 - Thermostats

- implemented as specified on the worksheet
- five functions to integrate the functionality

### Task 2 - Simulation of the Rayleigh-Taylor instability

periodic boundary condition is done by applying a precalculated correction vector to the position of the particles.

The Rayleigh taylor instability looks like the one on the worksheet

### Task 3 - Simulation of a falling drop - Liquid

checkpointing implemented

The drop creates a wave that continues through the periodic boundary or bounces up on the reflecting boundary.

### Task 4 & 5 - Performance Measurement and Profiling

We completely rewrote the linked cells, the first working version already contained all the optimizations we could think
of.
We now use SoA.
core operations are vectorized with std::experimental::simd
optimized LJF to not use sqrt for norm, because it's squared afterward
added precomputation for constant factors
checked that all function in the core part are inlined
switched to 2D calculations, because all tasks are 2D this week.

The the profiler showed that most of the time was spent in calculate_force.
Intel advisor showed almost no suggestions.
One advise about unaligned memory access, is not problematic and correcting the alignment there decreases performance.
Intel vtune reported no frontend issues and only backend issues, which means we utilize the hardware very much.

Because of inlining we couldn't get better profile infos for subsequent functions.

But we expect that the most significant time is spent iteration pairwise + calculate LJF and then sorting the particles.

On a laptop with AVX512 for the rayleigh-taylor instability scenario we get:

- 700 MUPS/s with clang++
- 466 MUPS/s with icpx
- 606 MUPS/s with g++

Clang++ outperformed in all our testings, but g++ and icpx could benefit from pgo, but we didn't finish implementing it
for them

On the cluster only g++ build, because we used very recent features of the c++ standard and g++ was the only compiler
that we could get to work

- 236,07 MUPS/s with g++

