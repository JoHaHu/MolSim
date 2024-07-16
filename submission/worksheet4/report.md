---
Group: D
Members: Julius Kramer, Tim Scholl, Johannes Hupe
PR: https://github.com/JoHaHu/MolSim/pull/70
---

# Report Group-D

### Media

All the media is stored in the git repository for this week because the files were too large for moodle even after packing. Images are also packed in the submission

Additionally, here is a OneDrive Link with the video files: https://1drv.ms/f/s!Ar3cyyDhVL_okssq7odwUhiSI0WrTQ?e=BcDDnv


## Tasks - Worksheet 4

### Task 1 - Thermostats

- Implemented as specified on the worksheet
- Five functions to integrate the functionality
- Tests for each functionality (cooling, heating and keeping temperature)

### Task 2 - Simulation of the Rayleigh-Taylor instability

- Periodic boundary condition is done by applying a precalculated correction vector to the position of the particles.
- The Rayleigh-Taylor instability looks like the one on the worksheet

### Task 3 - Simulation of a falling drop - Liquid

- Checkpointing implemented according to the worksheet. We used a new XSD schema for that
- The drop creates a wave that continues through the periodic boundary or bounces up on the reflecting boundary.

### Task 4 & 5 - Performance Measurement and Profiling

- We completely rewrote the linked-cells algorithm, the first working version already contained all the optimizations we
  could think
  of.
  We now use SoA.
- Core operations are vectorized with std::experimental::simd.
  optimized LJF to not use sqrt for norm, because it's squared afterward.
- Added precomputation for constant factors.
- Checked that all function in the core part are inlined.
- Switched to 2D calculations, because all tasks are 2D this week.

- The profiler showed that most of the time was spent in calculate_boundary_force.
- Intel advisor showed almost no suggestions.
- One advice about unaligned memory access is not problematic and correcting the alignment decreased performance.
- Intel VTune reported no frontend issues and only backend issues, which means we utilize the hardware very much.

Because of inlining we couldn't get better profile infos for subsequent functions.
But we expect that the most significant time is spent iteration pairwise + calculate LJF and then sorting the particles.

On a laptop with AVX512 for the rayleigh-taylor instability scenario we get:

- 700 MUPS/s with clang++
- 466 MUPS/s with icpx
- 606 MUPS/s with g++

Clang++ outperformed in all our testings, g++ and icpx could benefit from pgo, but we didn't finish implementing it
for them.

On the cluster only g++ built, because we used very recent features of the c++ standard and g++ was the only compiler
that we could get to work

- 236,07 MUPS/s with g++

