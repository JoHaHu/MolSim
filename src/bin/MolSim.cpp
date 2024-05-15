#include <spdlog/spdlog.h>
#include <vector>

#include "lib/config/config.h"
#include "lib/simulator/Simulator.h"
#include "lib/simulator/io/ParticleLoader.h"
#include "lib/simulator/io/VTKPlotter.h"
#include "lib/simulator/physics/ForceModel.h"
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

  auto particle_loader = simulator::io::ParticleLoader(config);
  auto [particle_container, force_model] = particle_loader.load_particles();

  auto plotter = std::make_unique<simulator::io::VTKPlotter>();

  auto simulator = simulator::Simulator(std::move(particle_container), std::move(plotter), config);

  if (config->io_interval == 0) {
    switch (force_model) {
      case simulator::physics::ForceModel::Gravity:
        simulator.run<simulator::physics::Gravity, false>();
        break;
      case simulator::physics::ForceModel::LennardJones:
        simulator.run<simulator::physics::LennardJones, false>();
        break;
    }
  } else {
    switch (force_model) {
      case simulator::physics::ForceModel::Gravity:
        simulator.run<simulator::physics::Gravity, true>();
        break;
      case simulator::physics::ForceModel::LennardJones:
        simulator.run<simulator::physics::LennardJones, true>();
        break;
    }
  }

  return 0;
}
