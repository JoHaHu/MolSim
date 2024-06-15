#pragma once

#include "Particle.h"
#include "config/Config.h"
#include "container/container.h"
#include "experimental/simd"
#include "simulator/io/Plotter.h"
#include "simulator/physics/ForceModel.h"
#include "simulator/physics/Gravity.h"
#include "simulator/physics/LennardJones.h"
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

template<const size_t DIMENSIONS>
class Simulator {
 private:
  container::ParticleContainer<DIMENSIONS> particles;
  physics::ForceModel physics;
  std::unique_ptr<io::Plotter<DIMENSIONS>> plotter;
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
      container::ParticleContainer<DIMENSIONS> &&particles,
      physics::ForceModel physics,
      std::unique_ptr<io::Plotter<DIMENSIONS>> &&plotter,
      const std::shared_ptr<config::Config> &config)
      : particles(std::move(particles)),
        physics(physics),
        plotter(std::move(plotter)),
        config(config),
        end_time(config->end_time),
        delta_t(config->delta_t){};

  auto calculate_position_particle(Particles<DIMENSIONS> &p, size_t index) const {

    const auto temp = pow(delta_t, 2) * (1 / (2 * p.mass[index]));
    for (int i = 0; i < DIMENSIONS; ++i) {
      p.positions[i][index] += delta_t * p.velocities[i][index] + temp * p.old_forces[i][index];
    }
  }
  auto calculate_velocity_particle(Particles<DIMENSIONS> &p, size_t index) const {

    const auto temp = delta_t * (1 / (2 * p.mass[index]));

    for (int i = 0; i < DIMENSIONS; ++i) {
      p.velocities[i][index] += temp * (p.old_forces[i][index] + p.forces[i][index]);
    }
  }

  template<typename F>
  auto calculate_force_particle_pair(F f, VectorizedParticle<DIMENSIONS> &p1, VectorizedParticle<DIMENSIONS> &p2, double_mask mask, std::array<double_v, DIMENSIONS> &correction) {

    std::array<double_v, DIMENSIONS> force{};
    f(p1, p2, mask, force, correction);

    for (int i = 0; i < DIMENSIONS; ++i) {
      where(mask, p2.force[i]) = p2.force[i] - force[i];
    }

    for (int i = 0; i < DIMENSIONS; ++i) {
      p1.force[i] += stdx::reduce(where(mask, force[i]), std::plus<>());
    }
  }

  /*! <p> Function for position calculation </p>
   *
   * calculates the position for all particles, takes no arguments and has no return value
   */
  auto calculate_position() -> void {
    SPDLOG_DEBUG("Updating positions");
    particles.linear([this](Particles<DIMENSIONS> &p, size_t index) {
      calculate_position_particle(p, index);
    });
  }

  /*! <p> Function for velocity calculation </p>
  * calculates the velocity for all particles, takes no arguments and has no return value
  */
  auto calculate_velocity() -> void {
    SPDLOG_DEBUG("Updating velocities");

    particles.linear([this](Particles<DIMENSIONS> &p, size_t index) {
      calculate_velocity_particle(p, index);
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

    switch (physics) {
      case physics::ForceModel::Gravity: {

        // TODO clarify what boundary forces should apply when simulating with gravity,
        //  since gravity will always pull towards the center it might be possible to have no boundary condition if we start from stationary
        particles.boundary(physics::lennard_jones::calculate_force);

        particles.refresh();
        particles.pairwise([this](auto &p1, auto &p2, auto mask, auto &correction) {
          calculate_force_particle_pair(physics::gravity::calculate_force_vectorized<DIMENSIONS>, p1, p2, mask, correction);
        });
        break;
      }
      case physics::ForceModel::LennardJones: {
        particles.boundary(physics::lennard_jones::calculate_force);
        particles.refresh();
        particles.pairwise([this](auto &p1, auto &p2, auto mask, auto &correction) {
          calculate_force_particle_pair(physics::lennard_jones::calculate_force_vectorized<DIMENSIONS>, p1, p2, mask, correction);
        });
        break;
      }
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
    particles.swap_force();
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