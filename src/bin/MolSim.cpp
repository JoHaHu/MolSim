#include "lib/simulator/io/FileReader.h"

#include <spdlog/spdlog.h>
#include <vector>

#include "lib/config/config.h"
#include "lib/simulator/Simulator.h"
#include "lib/simulator/io/VTKPlotter.h"
#include "lib/simulator/physics/Gravity.h"
#include "lib/simulator/physics/LennardJones.h"
#include "lib/utils/LoggerManager.h"

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

  std::vector<Particle> particles;
  try {
    spdlog::info("Reading file: {}", config->input_filename);
    simulator::io::FileReader::read_file(particles, config->input_filename);
  } catch (const std::exception &e) {
    spdlog::error("Failed to read file: {}", e.what());
    return 1;
  }

  auto plotter = std::make_unique<simulator::io::VTKPlotter>();
  auto particleContainer = ParticleContainer(particles);

  auto simulator = simulator::Simulator(std::move(particleContainer), std::move(plotter), config);

  if (config->io_interval == 0) {
    switch (config->task) {
      case simulator::Task::gravity:
        simulator.run<simulator::physics::Gravity, false>();
        break;
      case simulator::Task::collision:
        simulator.run<simulator::physics::LennardJones, false>();
        break;
    }
  } else {
    switch (config->task) {
      case simulator::Task::gravity:
        simulator.run<simulator::physics::Gravity, true>();
        break;
      case simulator::Task::collision:
        simulator.run<simulator::physics::LennardJones, true>();
        break;
    }
  }

  return 0;
}
