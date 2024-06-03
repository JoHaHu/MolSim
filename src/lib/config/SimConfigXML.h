//
// Created by TimSc on 03.06.2024.
//

#ifndef SIM_CONFIG_XML_H
#define SIM_CONFIG_XML_H

#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

using Cuboid = std::tuple<std::array<int, 3>, std::array<int, 3>, double, double, std::array<int, 3>>;

/**
 * All supported runtime configurations
 * */
class SimConfigXML {
 public:
  /**
   * build a config from provided XML file arguments parsed by XSD.
   * */
  SimConfigXML(std::string name, int output_frequency, double t_end, double delta_t, double epsilon, double sigma, double av_brown_motion, std::vector<Cuboid> cuboids);

  // Function to create and store the configuration values
  static auto store_config_values(std::string name, int output_frequency, double t_end, double delta_t, double epsilon, double sigma, double av_brown_motion, std::vector<Cuboid> cuboids) -> std::shared_ptr<SimConfigXML>;

  /**
   * the frequency of written output files
   * */
  std::string base_name;
  /**
   * the frequency of written output files
   * */
  double output_frequency{};
  /**
   * start time of the simulation, not used right now
   * */
  double start_time{};
  /**
   * the end time of the simulation
   * */
  double end_time{};
  /**
   * the timestep size
   * */
  double delta_t{};
  /**
   * the epsilon value used in the calculation of Lennard-Jones forces
   * */
  double epsilon{};
  /**
   * the sigma value used in the calculation of Lennard-Jones forces
   * */
  double sigma{};
  /**
   * average brownian motion velocity (if it should be specifically specified, otherwise initialised with 0.1)
   * */
  double average_brownian_motion = 0.1;
  /**
   * a vector that can store cuboids (defined for this simulation parameter/config storing purpose only)
   * */
  std::vector<Cuboid> cuboids;
  /**
   * the output file name
   * */
  std::string output_filename;
  /**
   * the input filename
   * */
  std::string input_filename;
  /**
   * the seed for the rng
   * */
  int seed{};
};

#endif// SIM_CONFIG_XML_H