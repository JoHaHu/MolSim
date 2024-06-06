---
Group: D
Members: Julius Kramer, Tim Scholl, Johannes Hupe
PR: https://github.com/JoHaHu/MolSim/pull/1
---

# Report Group-D

## Tasks (Worksheet 3)

### Task 1 - XML input

#### XSD Tool
- Installed the CodeSynthesis tool XSD XML Schema to C++ compiler 4.2.0
- Created the custom XML schema SimulationInputSchema.xsd to define the structure of an input XML
- The corresponding .hxx and .cxx files, we generated with the tool to create the XML binding to our application logic
- The XML schema defines one string, one positive integer and five double values (all non-repeated) for the simulation (name, output frequency, t_end, delta_t, epsilon, sigma and the average brownian motion velocity)
- Further, there are three arrays and two doubles per each cuboid (repeatable values) that describe cuboid specifications (coordinate, number of particles, distance h, mass m and initial velocity)


#### Processing
- The parsing and storing of information for further processing happens in the XMLFileReader class
- For that storing purpose, we created the SimConfigXML class which stores the data in appropriate attributes
- The simulation-bound parameters (single value for one simulation) are each stored in a separate variable
- The cuboid-bound parameters (single value/array per cuboid) are stored in a tuple and all tuples in a vector
- This enables multiple cuboids with different specifications to be created and stored
- The config file can then be passed to the simulation as we did in last week's worksheet when parsing from .txt files
- That enables the previous functionality to remain functional and adding that extra ability of input


### Task 2 - Linked-cell algorithm

- Implemented a new particle container using the linked-cell algorithm, retaining the existing direct sum implementation.
- Ensured container boundaries align with scenario descriptions and handled non-divisible cell widths.
- Adapted XML format to specify domain size and cutoff radius.
- Created a container to iterate over boundary and halo particles, supporting their deletion.
- Allowed iteration over particle pairs while utilizing Newton’s third law.
- Performed the "Collision of two bodies" simulation with specified parameters.
- Compared runtime per iteration for naive and linked-cell implementations.
- Visualized runtime comparisons for different molecule counts (1000, 2000, 4000, 8000) in a plot, embedded in Doxygen documentation.
    - 1000: TODO add results
    - 2000: TODO add results
    - 4000: TODO add results
    - 8000: TODO add results


### Task 3 - Boundary conditions

- Implemented outflow boundary condition: particles are removed when they cross the domain boundary.
- Implemented reflecting boundary condition: particles reflect off the boundary using "counter"-particles based on the Lennard-Jones potential.
- Adapted XML input to specify boundary conditions for each domain boundary.
- Ensured particles remain inside the domain and maintained the qualitative behavior of particles.
- Performed unit tests for each boundary condition to ensure correctness and stability.
- The time step delta t must be small to ensure stability because larger steps can cause particles to cross boundaries before the reflecting forces are accurately applied.


### Task 4 - Simulation of a falling drop - Wall

- Extended the particle generator to generate 2D discs packed with molecules arranged in a grid.
    - Disc properties:
        - Center coordinates
        - Initial velocity
        - Radius in terms of the number of molecules along the radius
        - Meshwidth (distance between molecules)
- Adapted the XML input format to include the new disc properties.
- Conducted the simulation with a disc flying against a reflecting boundary using specified parameters.
    - Parameters: xcenter = {60, 25, 0}, v = {0, −10, 0}, m = 1.0, R = 15, h ≈ 1.1225, ϵ = 5, σ = 1, ∆t = 0.00005, tend = 10, rcutoff = 3.0, domain size = {120, 50, 1}
- Additionally implemented generators for 3D shapes: sphere, torus, and double helix.
- Simulated experiments to test the new implementations. (See videos)


