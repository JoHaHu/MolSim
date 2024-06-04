
#include <spdlog/spdlog.h>
#include <variant>

#include "config/config.h"
#include "simulator/Simulator.h"
#include "simulator/io/ParticleLoader.h"
#include "simulator/io/VTKPlotter.h"
#include "simulator/physics/ForceModel.h"
#include "simulator/physics/Gravity.h"
#include "simulator/physics/LennardJones.h"
#include "utils/LoggerManager.h"

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
  LoggerManager::setup_logger(*config);

  auto particle_loader = simulator::io::ParticleLoader(config);

  auto [particles, force_model] = particle_loader.load_particles();

  auto particle_container = container::linked_cell<container::index::simple_index>(
      {140, 90, 1},
      3.0,
      {container::boundary_condition::outflow,
       container::boundary_condition::outflow,
       container::boundary_condition::outflow,
       container::boundary_condition::outflow,
       container::boundary_condition::reflecting,
       container::boundary_condition::reflecting},
      particles.size(), 1.0);
  //  auto particle_container = std::vector<Particle>();
  auto unique_plotter = std::make_unique<simulator::io::VTKPlotter>(config);

  auto container = container::particle_container(std::move(particle_container));

  for (auto &p : particles) {
    container.insert(p);
  }

  simulator::physics::force_model physics;
  switch (force_model) {
    case simulator::physics::ForceModel::Gravity:
      physics = {simulator::physics::Gravity()};
      break;
    case simulator::physics::ForceModel::LennardJones:
      physics = {simulator::physics::LennardJones(3.0, 1, 5)};
      break;
  }

  auto simulator = simulator::Simulator(std::move(container), physics, std::move(unique_plotter), config);

  auto startTime = std::chrono::high_resolution_clock::now();

  if (config->io_interval == 0) {
    simulator.run<false>();
  } else {
    simulator.run<true>();
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
