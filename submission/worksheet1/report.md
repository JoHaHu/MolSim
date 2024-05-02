---
Group: D
Members: Julius Kramer, Tim Scholl, Johannes Hupe
PR: https://github.com/JoHaHu/MolSim/pull/1
---

# Report Group-D

## Tasks

### Task 1

- followed instructions
- decide on using suggested code sanitizer, static analysis and debugging tools

### Task 2

- just followed instructions

### Task 3

- calculations based on the 3 formulas provided in the slide.
- Optimization: forces are calculated in pairs. This saves computation because it only needs to be calculated once.

### Task 4

- We suspect the object to be the sun, jupiter, the earth and Halley's Comet. This is based on the masses of the objects
  and the trajectories.
- We observed that the force vector of the sun and jupiter always point at each other as expected by Newtons laws, but
  get slightly offset based on the location of the other objects.

### Task 5

- The iterator patter looks like a good abstraction for the particle container.
- the strategy pattern can be used to abstract over the different kinds of force calculations
- We've switched to using a vector for now to store particles.
