//
// Created by TimSc on 03.06.2024.
//

#include "Config.h"
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
 * @return std::shared_ptr<Config> containing the configured parameters.
 *
 */


auto Config::store_config_values(SimulationType simulationType, BodyType bodyType, std::string &baseName, double endTime, double outputFrequency,
                                       std::string &outputFilename, int totalBodies, double deltaT, double inputSigma, double inputEpsilon, double massM, double distanceH,
                                       double averageBrownianMotion, std::vector<CelestialBody> celestialBodies, std::vector<Cuboid> cuboidVector,
                                       std::vector<Disc> discVector, int RngSeed) -> std::shared_ptr<Config> {
  if (simulationType == SimulationType::LennardJones) {
    Config config = Config();
      config.body_type = bodyType;
      config.simulation_type = simulationType;
      config.base_name = baseName;
      config.end_time = endTime;
      config.output_frequency = outputFrequency;
      config.output_filename = outputFilename;
      config.delta_t = deltaT;
      config.sigma = inputSigma;
      config.epsilon = inputEpsilon;
      config.mass_m = massM;
      config.distance_h = distanceH;
      config.brownian_motion = averageBrownianMotion;
      config.seed = RngSeed;
      if (bodyType == BodyType::Cub) {
        config.cuboids = std::move(cuboidVector);
      }
      else {
        config.discs = std::move(discVector);
      }
      return std::make_shared<Config>(config);
  }
  if (simulationType == SimulationType::Gravity) {
    Config config = Config();
    config.body_type = bodyType;
    config.simulation_type = simulationType;
    config.base_name = baseName;
    config.end_time = endTime;
    config.output_frequency = outputFrequency;
    config.output_filename = outputFilename;
    config.total_bodies = totalBodies;
    config.celestial_bodies = std::move(celestialBodies);
    config.seed = RngSeed;
    return std::make_shared<Config>(config);
  }
  return nullptr;
}