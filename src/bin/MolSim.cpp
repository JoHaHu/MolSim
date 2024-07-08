
#include <spdlog/spdlog.h>
#include <variant>

#include "config/Config.h"
#include "container/LinkedCell.h"
#include "simulator/Simulator.h"
#include "simulator/io/ParticleGenerator.h"
#include "simulator/io/VTKPlotter.h"
#include "simulator/io/xml_reader/XMLFileReader.h"
#include "simulator/physics/Force.h"
#include "simulator/physics/ForceModel.h"
#include "simulator/physics/LennardJones.h"
#include "utils/LoggerManager.h"

template<const size_t DIMENSIONS>
auto setup(std::shared_ptr<config::Config> config) -> auto{
  auto particle_loader = simulator::io::ParticleGenerator<DIMENSIONS>(config);

  auto particles_vector = particle_loader.load_particles();

  auto plotter = std::make_unique<simulator::io::VTKPlotter<DIMENSIONS>>(config);

  auto checkpointer = Checkpointer<DIMENSIONS>();

  std::unique_ptr<container::Container<DIMENSIONS>> pc;
  switch (config->particle_container_type) {
    case ParticleContainerType::Vector:
      // TODO
      break;
    case ParticleContainerType::LinkedCells:

      std::array<double, DIMENSIONS> domain;
      std::array<BoundaryCondition, 2 * DIMENSIONS> boundary;

      for (size_t i = 0; i < DIMENSIONS; ++i) {
        domain[i] = config->domain_size[i];
        boundary[i] = config->boundary_conditions[i];
        boundary[DIMENSIONS + i] = config->boundary_conditions[DIMENSIONS + i];
      }

      switch (config->simulation_type) {
        case simulator::physics::ForceModel::LennardJones: {
          auto physics = simulator::physics::LennardJonesForce(config->cutoff_radius, config->epsilon, config->sigma);
          pc = std::make_unique<container::LinkedCell<simulator::physics::LennardJonesForce, DIMENSIONS>>(
              std::move(physics),
              domain,
              config->cutoff_radius,
              boundary,
              config->sigma);
          break;
        }
        default:
          SPDLOG_ERROR("Unsupported Force with Linked Cells");
          exit(1);
      }

      break;
  }

  for (auto p : particles_vector) {
    pc->insert(p);
  }

  for (auto &cp_file : config->input_checkpoints) {
    auto cp = checkpointer.load_checkpoint(cp_file);
    for (auto p : cp) {
      pc->insert(p);
    }
  }

  pc->refresh();

  auto simulator = simulator::Simulator<DIMENSIONS>(
      std::move(pc),
      std::move(plotter),
      config,
      checkpointer);
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

  // Needs to be initialized and get overwritten
  std::chrono::high_resolution_clock::time_point startTime;

  if (config->dimensions == 2) {

    auto simulator = setup<2>(config);
    startTime = std::chrono::high_resolution_clock::now();
    if (config->output_frequency == 0) {
      simulator.run<false>();
    } else {
      simulator.run<true>();
    }
  } else if (config->dimensions == 3) {

    auto simulator = setup<3>(config);
    startTime = std::chrono::high_resolution_clock::now();
    if (config->output_frequency == 0) {
      simulator.run<false>();
    } else {
      simulator.run<true>();
    }
  } else {
    SPDLOG_ERROR("Unsupported dimension");
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
