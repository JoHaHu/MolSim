#pragma once

#include "Particle.h"
#include "config/Config.h"
#include "container/container.h"
#include "experimental/simd"
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
  container::ParticleContainer particles;
  physics::ForceModel physics;
  std::unique_ptr<io::Plotter> plotter;
  std::shared_ptr<config::Config> config;

  double end_time;
  double delta_t;
  unsigned long iteration = 0;

 public:
  /**
   * \param particles the particle container
   * \param plotter An instance plotter.
   * \param config the runtime configuration
   * */
  explicit Simulator(
      container::ParticleContainer &&particles,
      physics::ForceModel physics,
      std::unique_ptr<io::Plotter> &&plotter,
      const std::shared_ptr<config::Config> &config)
      : particles(std::move(particles)),
        physics(physics),
        plotter(std::move(plotter)),
        config(config),
        end_time(config->end_time),
        delta_t(config->delta_t) {};

  auto calculate_position_particle(
      double &position_x, double &position_y, double &position_z,
      double const &velocity_x, double const &velocity_y, double const &velocity_z,
      double const &old_force_x, double const &old_force_y, double const &old_force_z,
      double const &mass) const {

    const auto temp = pow(delta_t, 2) * (1 / (2 * mass));
    position_x += delta_t * velocity_x + temp * old_force_x;
    position_y += delta_t * velocity_y + temp * old_force_y;
    position_z += delta_t * velocity_z + temp * old_force_z;
  }
  auto calculate_velocity_particle(
      double &velocity_x, double &velocity_y, double &velocity_z,
      double const &force_x, double const &force_y, double const &force_z,
      double const &old_force_x, double const &old_force_y, double const &old_force_z,
      double const &mass) const {

    const auto temp = delta_t * (1 / (2 * mass));
    velocity_x += temp * (old_force_x + force_x);
    velocity_y += temp * (old_force_y + force_y);
    velocity_z += temp * (old_force_z + force_z);
  }

  template<typename F>
  auto calculate_force_particle_pair(F f, VectorizedParticle &p1, VectorizedParticle &p2, double_mask mask) {

    const auto force = f(p1, p2, mask);

    where(mask, p2.force[0]) = p2.force[0] - force[0];
    where(mask, p2.force[1]) = p2.force[1] - force[1];
    where(mask, p2.force[2]) = p2.force[2] - force[2];

    p1.force[0] += stdx::reduce(where(mask, force[0]), std::plus<>());
    p1.force[1] += stdx::reduce(where(mask, force[1]), std::plus<>());
    p1.force[2] += stdx::reduce(where(mask, force[2]), std::plus<>());
  }

  /*! <p> Function for position calculation </p>
   *
   * calculates the position for all particles, takes no arguments and has no return value
   */
  auto calculate_position() -> void {
    SPDLOG_DEBUG("Updating positions");
    particles.linear([this](Particles &p, size_t index) {
      calculate_position_particle(
          p.position_x[index], p.position_y[index], p.position_z[index],
          p.velocity_x[index], p.velocity_y[index], p.velocity_z[index],
          p.old_force_x[index], p.old_force_y[index], p.old_force_z[index],
          p.mass[index]);
    });
  }

  /*! <p> Function for velocity calculation </p>
  * calculates the velocity for all particles, takes no arguments and has no return value
  */
  auto calculate_velocity() -> void {
    SPDLOG_DEBUG("Updating velocities");

    particles.linear([this](Particles &p, size_t index) {
      calculate_velocity_particle(
          p.velocity_x[index], p.velocity_y[index], p.velocity_z[index],
          p.force_x[index], p.force_y[index], p.force_z[index],
          p.old_force_x[index], p.old_force_y[index], p.old_force_z[index],
          p.mass[index]);
    });
  }
  /**
  * @brief Calculates forces between particles.
  *
  * Resets forces for all particles, then calculates and updates forces for each particle pair.
  */

  auto calculate_force() -> void {
    SPDLOG_DEBUG("Starting force calculation");

    particles.swap_force();

    //    particles.boundary([this](Particles &p, auto index) {
    //      calculate_force_particle_pair(p, index);
    //    });
    particles.refresh();

    switch (physics) {
      case physics::ForceModel::Gravity:
        particles.pairwise([this](auto &p1, auto &p2, auto mask) {
          calculate_force_particle_pair(physics::gravity::calculate_force, p1, p2, mask);
        });
        break;
      case physics::ForceModel::LennardJones:
        particles.pairwise([this](auto &p1, auto &p2, auto mask) {
          calculate_force_particle_pair(physics::lennard_jones::calculate_force, p1, p2, mask);
        });
        break;
    }

    SPDLOG_TRACE("Force calculation completed.");
  };

  /**
    * @brief Runs the simulation.
    *
    * Executes the main simulation loop, updating particle positions, forces, and velocities
    * until the end time is reached. Periodically plots the particles based on the IO interval.
    */
  template<bool IO>
  auto run() -> void {
    spdlog::info("Running simulation...");
    double current_time = 0;
    iteration = 0;
    auto interval = config->output_frequency;

    calculate_force();
    // Plot initial position and forces

    if (IO) {
      plotter->plotParticles(particles, iteration);
      SPDLOG_DEBUG("Iteration {} plotted.", iteration);
    }

    while (current_time < end_time) {
      calculate_position();
      calculate_force();
      calculate_velocity();

      iteration++;
      if (IO && iteration % interval == 0) {
        plotter->plotParticles(particles, iteration);
        SPDLOG_DEBUG("Iteration {} plotted.", iteration);
      }

      SPDLOG_DEBUG("Iteration {} finished.", iteration);

      current_time += delta_t;
    }
    SPDLOG_INFO("Output written. Terminating...");
  }
};

}// namespace simulator