
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

  auto particle_loader = simulator::io::ParticleGenerator(config);

  auto particles_vector = particle_loader.load_particles();

  auto plotter = std::make_unique<simulator::io::VTKPlotter>(config);

  simulator::physics::ForceModel physics = config->simulation_type;
  switch (config->simulation_type) {
    case simulator::physics::ForceModel::Gravity:
      break;
    case simulator::physics::ForceModel::LennardJones:
      simulator::physics::lennard_jones::initialize_constants(config->epsilon, config->sigma, config->cutoff_radius);
      break;
  }

  auto particles = Particles();

  for (auto p : particles_vector) {
    particles.insert_particle(p);
  }

  container::particle_container pc = container::particle_container(particles);

  switch (config->particle_loader_type) {
    case ParticleContainerType::Vector:
      break;
    case ParticleContainerType::LinkedCells:
      auto lc = container::LinkedCell(
          std::move(particles),
          config->domain_size,
          config->cutoff_radius,
          config->boundary_conditions,
          config->sigma);

      pc = container::particle_container(std::move(lc));
      break;
  }

  pc.refresh();
  auto simulator = simulator::Simulator(std::move(pc), physics, std::move(plotter), config);

  auto startTime = std::chrono::high_resolution_clock::now();

  if (config->output_frequency == 0) {
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

  SPDLOG_INFO("Simulation completed in {:02d}:{:02d}:{:03d} (mm:ss:ms)", minutes.count(), seconds.count(), milliseconds);

  return 0;
}
