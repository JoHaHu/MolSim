//
// Created by template from MolSim
//

#include "config/Config.h"
#include "simulator/Simulator.h"
#include "simulator/io/ParticleGenerator.h"
#include "simulator/io/VTKPlotter.h"
#include "simulator/io/xml_reader/XMLFileReader.h"
#include "simulator/physics/ForceModel.h"
#include "simulator/physics/Gravity.h"
#include "simulator/physics/LennardJones.h"
#include "utils/LoggerManager.h"
#include <spdlog/spdlog.h>
#include <vector>

/**
 * Main entry point of the MolSim particle simulation program.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @return 0 if successful, 1 on error.
 *
 * Orchestrates initialization, input validation, and simulation. Logs final output or errors.
 */
auto main(int argc, char *argv[]) -> int {

  auto config = config::Config::parse_config(argc, argv);
  LoggerManager::setup_logger(config);

  auto startTime = std::chrono::high_resolution_clock::now();

  auto particle_loader = simulator::io::ParticleGenerator(config);

  auto particles = particle_loader.load_particles();
  auto particle_container = ParticleContainer(particles);
  auto plotter = std::make_unique<simulator::io::VTKPlotter>(config);

  auto simulator = simulator::Simulator(particle_container, std::move(plotter), config);

  if (config->output_frequency == 0) {
    switch (config->simulation_type) {
      case ForceModel::Gravity:
        simulator.run<simulator::physics::Gravity, false>();
        break;
      case ForceModel::LennardJones:
        simulator.run<simulator::physics::LennardJones, false>();
        break;
    }
  } else {
    switch (config->simulation_type) {
      case ForceModel::Gravity:
        simulator.run<simulator::physics::Gravity, true>();
        break;
      case ForceModel::LennardJones:
        simulator.run<simulator::physics::LennardJones, true>();
        break;
    }
  }

  auto endTime = std::chrono::high_resolution_clock::now();
  auto durationMs = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
  auto minutes = std::chrono::duration_cast<std::chrono::minutes>(durationMs);
  durationMs -= minutes;
  auto seconds = std::chrono::duration_cast<std::chrono::seconds>(durationMs);
  durationMs -= seconds;
  auto milliseconds = durationMs.count();

  spdlog::info("Simulation completed in {:02d}:{:02d}:{:03d} (mm:ss:ms)", minutes.count(), seconds.count(), milliseconds);

  return 0;
}
