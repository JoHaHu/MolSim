---
Group: D
Members: Julius Kramer, Tim Scholl, Johannes Hupe
PR: https://github.com/JoHaHu/MolSim/pull/1
---

# Report Group-D

## Tasks

### Task 1

- followed instructions
- decided on using suggested code sanitizer, static analysis and debugging tools as well as .py script for ParaView

### Task 2

- just followed instructions

### Task 3

- calculations based on the 3 formulas provided in the slide.
- implementation of calculation is in MolSim.cpp in the respective functions
- Optimizations we did: forces are calculated in pairs. This saves computation because it only needs to be calculated once.

### Task 4

- We suspect the celestial bodies to be the sun, jupiter, the earth and Halley's Comet (compare screenshot for reference). This is based on the masses of the objects
  and their trajectories in the simulation.
- We observed that the force vector of the sun and jupiter always point at each other as expected by Newtons laws, but
  get slightly offset based on the location of the other objects.

### Task 5

- The iterator patter looks like a good abstraction for the particle container.
- We implemented the iterator in ParticleContainer.cpp and ParticleContainer.h
- While implementing we encountered some compiling issues e.g. when using const iterators and some issues that we broke down to local setup differences
- The strategy pattern can be used to abstract over the different kinds of force calculations
- We have switched to using a vector for now to store particles.
- We also added the UML class diagram generation to the Doxyfile which creates a class diagram according to the project, this works by using Graphviz and Doxygen
