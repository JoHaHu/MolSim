---
Group: D
Members: Julius Kramer, Tim Scholl, Johannes Hupe
PR: https://github.com/JoHaHu/MolSim/pull/1
---

# Report Group-D

## Tasks

### Task 1 - Unit Testing

#### Testing suite setup
- Integrated Googletest Suite via FetchContent into CMake so no further user installation is necessary
- For CMake integration, we introduced a new CMakeLists structure and split the config content up in several files
- Configuration for googletest inside own test/CMakeLists.txt file
- Tests are all runnable by using ctest (and also integrated in CI)

#### Previous implementation
- Implemented tests for previous week's functionality
- Tested ParticleContainer / iterator functionality with 6 unit tests
- Tested gravitational force calculation with manually pre-calculated expected results

#### New implementation - molecular collision
- New tests for Lennard-Jones forces, also with manually pre-calculated expected results
- Testing behavior of calculate_force function to confirm proper calculation
- Used mainly/only expects and not asserts because one non-matching value does not mean that everything else is wrong, too
- Expects and therefore continuation of testing can give valuable feedback about where exactly a mistake could have happened


### Task 2 - Continuous Integration
- Integrated several checks in GitHub with the provided continuous integration features

#### Each push is checked for the following:
- Code quality is checked with clang-tidy
- Unit tests of the tests folder are run to make sure the functionality is properly working
- Formatting of the code is checked with clang-tidy
- Code sanitizers with clang and gcc


### Task 3 - Logging

- Chose spdlog functions due to maintainability and debug functionality
- We intend to use macros in the future for zero runtime overhead, better performance
- Included spdlog in CMake project 
- Set up LoggerManager for global logger handling
- Implemented multi level logging functionality, where deemed necessary  


### Task 4 - Collision Simulation

- Created a function to generate particles arranged in a 3D grid
- Checked for and handled invalid data in the input cuboid parameters
- Calculated particle positions based on grid dimensions and spacing
- Initialized particle velocities with a predefined distribution.
- Stored the generated particles efficiently in a vector



