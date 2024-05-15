#pragma once

#include "lib/Particle.h"
#include "lib/ParticleContainer.h"
#include "lib/config/config.h"
#include "lib/simulator/io/Plotter.h"
#include "lib/utils/ArrayUtils.h"
#include <cmath>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <utility>

namespace simulator {

template<typename T>
/** <p> physics concept for force calculation </p>
    * calculates the force for all particles,
    * \param T
    * \param particle1
    * \param particle2
    */
concept Physics = requires(T type, Particle const &particle1, Particle const &particle2) {
  { T::calculate_force(particle1, particle2) } -> std::convertible_to<std::array<double, 3>>;
};

class Simulator {

 private:
  ParticleContainer particles;
  std::unique_ptr<io::Plotter> plotter;
  std::shared_ptr<config::Config> config;
  double start_time;
  double end_time;
  double delta_t;

 public:
  /**
   * \param particles
   * \param plotter
   * \param config
   * */
  explicit Simulator(
      ParticleContainer particles,
      std::unique_ptr<io::Plotter> plotter,
      const std::shared_ptr<config::Config> &config);

  template<Physics PY, bool IO>
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
  template<Physics PY>
  auto calculate_force() -> void;
};

Simulator::Simulator(
    ParticleContainer particles, std::unique_ptr<io::Plotter> plotter, const std::shared_ptr<config::Config> &config)
    : particles(std::move(particles)), plotter(std::move(plotter)), config(config), start_time(config->start_time), end_time(config->end_time), delta_t(config->delta_t) {
}

template<Physics PY, bool IO>
auto Simulator::run() -> void {
  spdlog::info("Running simulation...");
  double current_time = start_time;
  int iteration = 0;
  auto interval = config->io_interval;

  calculate_force<PY>();

  while (current_time < end_time) {
    calculate_position();
    calculate_force<PY>();
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

void Simulator::calculate_position() {
  spdlog::debug("Updating positions for {} particles.", particles.size());
  for (auto &particle : particles) {
    particle.position = particle.position + delta_t * particle.velocity + pow(delta_t, 2) * (1 / (2 * particle.mass)) * particle.old_force;
    spdlog::trace("Particle position updated: ({}, {}, {})", particle.position[0], particle.position[1], particle.position[2]);
  }
}

void Simulator::calculate_velocity() {
  spdlog::debug("Updating velocities for {} particles.", particles.size());
  for (auto &particle : particles) {
    particle.velocity = particle.velocity + delta_t * (1 / (2 * particle.mass)) * (particle.old_force + particle.force);
    spdlog::trace("Particle velocity updated: ({}, {}, {})", particle.velocity[0], particle.velocity[1], particle.velocity[2]);
  }
}

template<Physics PY>
void Simulator::calculate_force() {
  spdlog::debug("Starting force calculation for {} particles.", particles.size());
  for (auto &particle : particles) {
    particle.old_force = particle.force;
    particle.force = {0, 0, 0};
  }

  for (auto pair = particles.begin_pair(); pair != particles.end_pair(); pair++) {
    const auto [particle1, particle2] = *pair;

    const auto force = PY::calculate_force(particle1, particle2);

    particle1.force = particle1.force + force;
    particle2.force = particle2.force - force;

    spdlog::trace("Force updated for particle pair: ({}, {}, {}) - ({}, {}, {})", particle1.force[0], particle1.force[1], particle1.force[2], particle2.force[0], particle2.force[1], particle2.force[2]);
  }
  spdlog::trace("Force calculation completed.");
}

}// namespace simulator