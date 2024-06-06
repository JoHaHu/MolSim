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

#### New particle container
- to be completed


### Task 3 - Boundary conditions

- to be completed


### Task 4 - Simulation of a falling drop - Wall

- actual simulation task: to be executed and documented here



