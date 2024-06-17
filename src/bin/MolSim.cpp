
#include <spdlog/spdlog.h>
#include <variant>

#include "config/Config.h"
#include "simulator/Simulator.h"
#include "simulator/io/ParticleGenerator.h"
#include "simulator/io/VTKPlotter.h"
#include "simulator/io/xml_reader/XMLFileReader.h"
#include "simulator/physics/ForceModel.h"
#include "simulator/physics/LennardJones.h"
#include "utils/LoggerManager.h"

template<const size_t DIMENSIONS>
auto setup(std::shared_ptr<config::Config> config) -> auto {
  auto particle_loader = simulator::io::ParticleGenerator<DIMENSIONS>(config);

  auto particles_vector = particle_loader.load_particles();

  auto plotter = std::make_unique<simulator::io::VTKPlotter<DIMENSIONS>>(config);

  simulator::physics::ForceModel physics = config->simulation_type;
  switch (config->simulation_type) {
    case simulator::physics::ForceModel::Gravity:
      break;
    case simulator::physics::ForceModel::LennardJones:
      simulator::physics::lennard_jones::initialize_constants(config->epsilon, config->sigma, config->cutoff_radius);
      break;
  }

  container::ParticleContainer pc = container::ParticleContainer<DIMENSIONS>(Particles<DIMENSIONS>());

  switch (config->particle_container_type) {
    case ParticleContainerType::Vector:
      break;
    case ParticleContainerType::LinkedCells:

      std::array<double, DIMENSIONS> domain;
      std::array<BoundaryCondition, 2 * DIMENSIONS> boundary;

      for (int i = 0; i < DIMENSIONS; ++i) {
        domain[i] = config->domain_size[i];
        boundary[i] = config->boundary_conditions[i];
        boundary[DIMENSIONS + i] = config->boundary_conditions[DIMENSIONS + i];
      }

      auto lc = container::LinkedCell<DIMENSIONS>(
          domain,
          config->cutoff_radius,
          boundary,
          config->sigma);

      pc = container::ParticleContainer<DIMENSIONS>(std::move(lc));
      break;
  }

  for (auto p : particles_vector) {
    pc.insert(p);
  }

  pc.refresh();

  auto simulator = simulator::Simulator<DIMENSIONS>(std::move(pc), physics, std::move(plotter), config);
  return simulator;
}

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

  auto simulator = setup<2>(config);

  const auto startTime = std::chrono::high_resolution_clock::now();

  if (config->output_frequency == 0) {
    simulator.run<false>();
  } else {
    simulator.run<true>();
  }
  auto endTime = std::chrono::high_resolution_clock::now();
  auto durationMs = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
  auto rate = (double) durationMs.count() / (config->end_time / config->delta_t);
  auto minutes = std::chrono::duration_cast<std::chrono::minutes>(durationMs);
  durationMs -= minutes;
  auto seconds = std::chrono::duration_cast<std::chrono::seconds>(durationMs);
  durationMs -= seconds;
  auto milliseconds = durationMs.count();

  SPDLOG_INFO("Simulation completed in {:02d}:{:02d}:{:03d} (mm:ss:ms) | {:0.5f} ms/iteration | {:0.5f} MUPS/s", minutes.count(), seconds.count(), milliseconds, rate, 1000 / rate);

  return 0;
}
