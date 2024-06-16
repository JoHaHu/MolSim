#pragma once

#include "CelestialBody.h"
#include "Cuboid.h"
#include "Disc.h"
#include "DoubleHelix.h"
#include "Sphere.h"
#include "Torus.h"
#include "container/boundary.h"
#include "simulator/physics/ForceModel.h"
#include <array>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

/**
 * definition of an enum class that gives information about the desired ParticleLoader
 */
enum class ParticleContainerType {
  LinkedCells,
  Vector
};

/**
 * definition of an enum class that gives information about the desired simulated body type in Lennard Jones
 */
enum class BodyType {
  Cub,
  Dis,
  Sph,
  Tor,
  Hel
};

/**
 * All supported runtime configurations
 */
namespace config {

/**
  * All supported runtime configurations
  */
class Config {

 public:
  // Default constructor
  Config() = default;

  /**
   * the frequency of written output files
   */
  std::string base_name;

  /**
   * the end time of the simulation
   */
  double end_time{};

  /**
   * the frequency of written output files
   */
  int output_frequency{};

  /**
   * the output file name
   */
  std::string output_filename;

  /**
   * an enum value that gives information about the desired simulation type
   */
  simulator::physics::ForceModel simulation_type{};

  /**
   * an enum value that gives information about the desired simulation type
   */
  ParticleContainerType particle_loader_type{};

  /**
   * an enum value that gives information about the desired simulated body type in Lennard Jones
   */
  BodyType body_type{};

  /**
   * an array to store the domain size for the linked-cells algorithm
   */
  std::array<double, 3> domain_size{};

  std::array<BoundaryCondition, 6> boundary_conditions;

  /**
   * the total amount of celestial bodies in the gravitational planetary simulation
   */
  int total_bodies{};

  /**
   * the timestep size
   */
  double delta_t{};

  /**
   * the sigma value used in the calculation of Lennard-Jones forces
   */
  double sigma{};

  /**
   * the epsilon value used in the calculation of Lennard-Jones forces
   */
  double epsilon{};

  /**
   * the mass of one particle in a disc or cuboid
   */
  double mass_m{};

  /**
   * the distance h between particles (mesh width of the grid) in cuboid or disc simulation
   */
  double distance_h{};

  /**
   * average brownian motion velocity
   */
  double brownian_motion{};

  /**
   * the cutoff radius for the linked-cells algorithm
   */
  double cutoff_radius{};

  /**
   * a vector that can store multiple celestial bodies for simulation defined in the CelestialBody class
   */
  std::vector<CelestialBody> celestial_bodies;

  /**
   * a vector that can store multiple cuboids for simulation defined in the Cuboid class
   */
  std::vector<Cuboid> cuboids;

  /**
   * a vector that can store multiple discs for simulation defined in the Disc class
   */
  std::vector<Disc> discs;

  /**
   * a vector that can store multiple spheres for simulation defined in the Sphere class
   */
  std::vector<Sphere> spheres;

  /**
   * a vector that can store multiple tori for simulation defined in the Torus class
   */
  std::vector<Torus> tori;

  /**
   * a vector that can store multiple double helices for simulation defined in the DoubleHelix class
   */
  std::vector<DoubleHelix> double_helices;

  /**
   * a random seed that is necessary for simulation
   */
  int seed{};

  /**
   * the input filename
   */
  std::string input_filename;

  /**
   * prints the help message
   */
  static auto print_usage() -> void;

  /**
   * build a config from provided cmdline arguments.
   */
  static auto parse_config(int argc, char *argv[]) -> std::shared_ptr<Config>;
};

}// namespace config
