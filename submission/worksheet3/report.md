---
Group: D
Members: Julius Kramer, Tim Scholl, Johannes Hupe
PR: https://github.com/JoHaHu/MolSim/pull/50
---

# Report Group-D

## Tasks (Worksheet 3)

### Task 1 - XML input

#### XSD Tool

#### XSD Tool
- Installed the CodeSynthesis tool XSD XML Schema to C++ compiler 4.2.0
- Created the custom XML schema SimulationInputSchema.xsd to define the structure of an input XML
- The corresponding .hxx and .cxx files, we generated with the tool to create the XML binding to our application logic
- The XML schema defines several input variables next to the ones given on the worksheet, e.g. like different choices between types of bodies, particle loader types etc.
- Further, we have defined different bodies like the cuboid or disc but also added spheres, tori (torus) and double_helices that are described by coordinates, velocity and other body-specific values like radius etc.
- Functionality of previous weeks has been preserved by creating corresponding XML schema elements and XML files for the ws1 and ws2 gravitational and two-body collision simulations

#### Processing
- The parsing and storing of information for further processing happens in the XMLFileReader class
- For that storing purpose, we created the Config class (which replaces the old Config class) which stores the data in appropriate attributes
- The class contains all attributes for every type of body, simulation and force models as well as deciding enums that tell us which type of action/simulation is desired
- Only the given input is set in the variables while others stay empty, which is fine because they will not be called in that situation
- Global simulation-bound parameters (single value for one simulation) are stored in a header and settings element which has attributes for these parameters for a slimmer XML file
- Specific cuboid-bound parameters (single value/array per cuboid) are stored in a corresponding class (e.g. Disc) and multiple instances in a vector of that class
- Enables multiple cuboids with different specifications to be created and stored
- Config file can then be passed to the simulation methods that pick out the necessary info
- Enables the previous functionality to remain functional and adding that extra ability of input

### Task 2 - Linked-cell algorithm

- Implemented linked cells with a linear array of cells containing references to particles.
- All particles are stored sequential in a memory arena, that allows to keep static references and mark elements as
  iterating and skip while iterating
- Compared runtime per iteration for naive and linked-cell implementations.
- Visualized runtime comparisons for different molecule counts (1000, 2000, 4000, 8000) in a plot, embedded in Doxygen
  documentation.
- Results:
    - 1000: 0:9 1:22
    - 2000: 0:18 5:12
    - 4000: 0:38 20:53
    - 8000: To be conducted

### Task 3 - Boundary conditions

- Implemented outflow boundary condition: particles are removed when they cross the domain boundary.
- Implemented reflecting boundary condition: particles reflect off the boundary using "counter"-particles based on the
  Lennard-Jones potential.
- Adapted XML input to specify boundary conditions for each domain boundary.
- Ensured particles remain inside the domain and maintained the qualitative behavior of particles.
- Performed unit tests for each boundary condition to ensure correctness and stability.
- The time step delta t must be small to ensure stability because larger steps can cause particles to cross boundaries
  before the reflecting forces are accurately applied.

### Task 4 - Simulation of a falling drop - Wall

- Extended the particle generator to generate 2D discs packed with molecules arranged in a grid.
    - Disc properties:
        - Center coordinates
        - Initial velocity
        - Radius in terms of the number of molecules along the radius
        - Meshwidth (distance between molecules)
- Adapted the XML input format to include the new disc properties.
- Conducted the simulation with a disc flying against a reflecting boundary using specified parameters.
    - Parameters: xcenter = {60, 25, 0}, v = {0, −10, 0}, m = 1.0, R = 15, h ≈ 1.1225, ϵ = 5, σ = 1, ∆t = 0.00005,
      tend = 10, rcutoff = 3.0, domain size = {120, 50, 1}
- Additionally implemented generators for 3D shapes: sphere, torus, and double helix.
- Simulated experiments to test the new implementations. (See videos)


