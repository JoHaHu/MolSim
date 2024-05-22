#pragma once

#include "Particle.h"
#include "config/config.h"
#include "container/container.h"
#include "simulator/io/Plotter.h"
#include "simulator/physics/ForceModel.h"
#include "utils/ArrayUtils.h"
#include "utils/variants.h"
#include <cmath>
#include <ranges>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <utility>

namespace simulator {

/**
 * The main Simulator class. can be configured by providing a config and a plotter. Some methods use a physics model provided at compile time.
 * */
class Simulator {
 private:
  container::particle_container particles;
  physics::force_model physics;
  std::unique_ptr<io::Plotter> plotter;
  std::shared_ptr<config::Config> config;

  double start_time;
  double end_time;
  double delta_t;

 public:
  /**
   * \param particles the particle container
   * \param plotter An instance plotter.
   * \param config the runtime configuration
   * */
  explicit Simulator(
      container::particle_container &&particles,
      physics::force_model physics,
      std::unique_ptr<io::Plotter> &&plotter,
      const std::shared_ptr<config::Config> &config);

  /**
  * @brief Runs the simulation.
  *
  * Executes the main simulation loop, updating particle positions, forces, and velocities
  * until the end time is reached. Periodically plots the particles based on the IO interval.
  */
  template<bool IO>
  auto run() -> void;

  /*! <p> Function for position calculation </p>
   *
   * calculates the position for all particles, takes no arguments and has no return value
   */
  auto calculate_position() -> void;

  /*! <p> Function for velocity calculation </p>
  * calculates the velocity for all particles, takes no arguments and has no return value
  */
  auto calculate_velocity() -> void;

  /**
  * @brief Calculates forces between particles.
  *
  * Resets forces for all particles, then calculates and updates forces for each particle pair.
  */

  auto calculate_force() -> void;
};

Simulator::Simulator(
    container::particle_container &&particles,
    physics::force_model physics,
    std::unique_ptr<io::Plotter> &&plotter,
    const std::shared_ptr<config::Config> &config)
    : particles(std::move(particles)),
      physics(physics),
      plotter(std::move(plotter)),
      config(config),
      start_time(config->start_time),
      end_time(config->end_time),
      delta_t(config->delta_t) {
}

/**
 * runs the main simulation loop. The physics concept must be specified to provide compile time transparency to allow compiler optimizations
 * */
template<bool IO>
auto Simulator::run() -> void {
  spdlog::info("Running simulation...");
  double current_time = start_time;
  int iteration = 0;
  auto interval = config->io_interval;

  calculate_force();

  while (current_time < end_time) {
    calculate_position();
    calculate_force();
    calculate_velocity();

    iteration++;
    if (iteration % interval == 0 && IO) {
      plotter->plotParticles(particles, iteration);
      spdlog::debug("Iteration {} plotted.", iteration);
    }

    spdlog::debug("Iteration {} finished.", iteration);

    current_time += delta_t;
  }
  spdlog::info("Output written. Terminating...");
}

inline void Simulator::calculate_position() {
  spdlog::debug("Updating positions for {} particles.", particles.size());
  for (auto &particle : particles) {
    particle.position = particle.position + delta_t * particle.velocity + pow(delta_t, 2) * (1 / (2 * particle.mass)) * particle.old_force;
    spdlog::trace("Particle position updated: ({}, {}, {})", particle.position[0], particle.position[1],
                  particle.position[2]);
  }
}

inline void Simulator::calculate_velocity() {
  spdlog::debug("Updating velocities for {} particles.", particles.size());
  for (auto &particle : particles) {
    particle.velocity = particle.velocity + delta_t * (1 / (2 * particle.mass)) * (particle.old_force + particle.force);
    spdlog::trace("Particle velocity updated: ({}, {}, {})", particle.velocity[0], particle.velocity[1],
                  particle.velocity[2]);
  }
}

inline void Simulator::calculate_force() {
  spdlog::debug("Starting force calculation for {} particles.", particles.size());
  for (auto &particle : particles) {
    particle.old_force = particle.force;
    particle.force = {0, 0, 0};
  }
  auto comb = particles.pairwise();
  std::visit(
      [this](auto &comb) { for (auto pair : comb) {

      Particle &particle1 = std::get<0>(pair);
      Particle &particle2 = std::get<1>(pair);
      const auto force = physics::calculate_force(physics, particle1, particle2);
      particle1.force = particle1.force + force;
      particle2.force = particle2.force - force;
      spdlog::trace("Force updated for particle pair: ({}, {}, {}) - ({}, {}, {})", particle1.force[0],
                    particle1.force[1], particle1.force[2], particle2.force[0], particle2.force[1],
                    particle2.force[2]);
    } },
      comb);

  spdlog::trace("Force calculation completed.");
}

}// namespace simulator