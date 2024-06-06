//
// Created by TimSc on 03.06.2024.
//

#ifndef SIM_CONFIG_XML_H
#define SIM_CONFIG_XML_H

#pragma once

#include "CelestialBody.h"
#include "Cuboid.h"
#include "Disc.h"
#include <array>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

/**
 * definition of an enum class that gives information about the desired simulation
 */
enum class ForceModel {
  Gravity,
  LennardJones
};

/**
 * definition of an enum class that gives information about the desired simulated body type in Lennard Jones
 */
enum class BodyType {
  Cub,
  Dis
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
  std::string base_name{};

  /**
   * the end time of the simulation
   */
  double end_time{};

  /**
   * the frequency of written output files
   */
  double output_frequency{};

  /**
   * the output file name
   */
  std::string output_filename{};

  /**
   * an enum value that gives information about the desired simulation type
   */
  ForceModel simulation_type{};

  /**
   * an enum value that gives information about the desired simulated body type in Lennard Jones
   */
  BodyType body_type{};

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
   * Static function to create the config by calling the respective constructors, setting the values and returning a shared_ptr to the config object
   */
  static auto store_config_values(ForceModel simulationType, BodyType bodyType, std::string &baseName, double endTime,
                                  double outputFrequency, std::string &outputFilename, int totalBodies, double deltaT,
                                  double inputSigma, double inputEpsilon, double massM, double distanceH,
                                  double averageBrownianMotion, std::vector<CelestialBody> celestialBodies,
                                  std::vector<Cuboid> cuboidVector, std::vector<Disc> discVector,
                                  int RngSeed) -> std::shared_ptr<Config>;

  /**
   * build a config from provided cmdline arguments.
   */
  static auto parse_config(int argc, char *argv[]) -> std::shared_ptr<Config>;
};

}// namespace config

#endif// SIM_CONFIG_XML_H
