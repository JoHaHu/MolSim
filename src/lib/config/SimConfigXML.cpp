//
// Created by TimSc on 03.06.2024.
//

#include "SimConfigXML.h"
#include "spdlog/spdlog.h"
#include <utility>

/**
 * @brief Creates an object to store the simulation parameters given by the input XML and parsed by XSD.
 *
 *
 * @param name The base name.
 * @param output_frequency The write frequency of output files.
 * @param t_end The simulation end time
 * @param delta_t The time steps of t.
 * @param epsilon The epsilon constant used in calculation of Lennard-Jones forces.
 * @param sigma The sigma constant used in calculation of Lennard-Jones forces.
 * @param av_brown_motion The average brownian motion velocity used for calculation.
 * @param cuboids A vector of cuboids that should be simulated.
 *
 * @return std::shared_ptr<SimConfigXML> containing the configured parameters.
 *
 */

// Constructor definition
SimConfigXML::SimConfigXML(std::string name, int output_frequency, double t_end, double delta_t, double epsilon, double sigma, double av_brown_motion, std::vector<Cuboid> cuboids)
    : base_name(std::move(name)),
      output_frequency(output_frequency),
      end_time(t_end),
      delta_t(delta_t),
      epsilon(epsilon),
      sigma(sigma),
      average_brownian_motion(av_brown_motion),
      cuboids(std::move(cuboids)){};

// Static function to create and return a shared_ptr to SimConfigXML
auto SimConfigXML::store_config_values(std::string name, int output_frequency, double t_end, double delta_t, double epsilon, double sigma, double av_brown_motion, std::vector<Cuboid> cuboids) -> std::shared_ptr<SimConfigXML> {
  return std::make_shared<SimConfigXML>(std::move(name), output_frequency, t_end, delta_t, epsilon, sigma, av_brown_motion, std::move(cuboids));
}
