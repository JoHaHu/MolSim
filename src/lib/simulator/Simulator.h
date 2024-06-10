#pragma once

#include "Particle.h"
#include "config/Config.h"
#include "container/container.h"
#include "simulator/io/Plotter.h"
#include "simulator/physics/ForceModel.h"
#include "simulator/physics/Thermostat.h"
#include "utils/ArrayUtils.h"
#include "utils/variants.h"
#include <cmath>
#include <ranges>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <utility>

namespace simulator {

// TODO @TIM Diese Variablen mit Config bzw. XML verstellbar machen.
// TODO @IRGENDJEMAND Die constexpr's hier raus löschen, und die aus der Config verwernden.
static constexpr bool use_brownian_motion = true;
static constexpr double Tinit = 300;
static constexpr double Ttarget = 300;
static constexpr double deltaT = 10;
static constexpr int nthermostat = 100;

/**
 * The main Simulator class. can be configured by providing a config and a plotter. Some methods use a physics model provided at compile time.
 * */
class Simulator {
 private:
  container::particle_container particles;
  physics::force_model physics;
  std::unique_ptr<io::Plotter> plotter;
  std::shared_ptr<config::Config> config;
  Thermostat thermostat;

  double end_time;
  double delta_t;
  unsigned long iteration = 0;

  auto calculate_position_particle(Particle &particle) const {
    particle.position = particle.position + delta_t * particle.velocity + pow(delta_t, 2) * (1 / (2 * particle.mass)) * particle.old_force;
    SPDLOG_TRACE("Particle position updated: ({}, {}, {})", particle.position[0], particle.position[1],
                 particle.position[2]);
  }
  auto calculate_velocity_particle(Particle &particle) const {
    particle.velocity = particle.velocity + delta_t * (1 / (2 * particle.mass)) * (particle.old_force + particle.force);
    SPDLOG_TRACE("Particle velocity updated: ({}, {}, {})", particle.velocity[0], particle.velocity[1],
                 particle.velocity[2]);
  }

  void calculate_old_force_particle(Particle &particle) {
    particle.old_force = particle.force;
    particle.force = {0, 0, 0};
  }

  auto calculate_force_particle_pair(std::tuple<Particle &, Particle &> pair) {
    auto [particle1, particle2] = pair;

    const auto force = physics::calculate_force(physics, particle1, particle2);

    particle1.force = particle1.force + force;
    particle2.force = particle2.force - force;
    SPDLOG_TRACE("Force updated for particle pair: ({}, {}, {}) - ({}, {}, {})", particle1.force[0],
                 particle1.force[1], particle1.force[2], particle2.force[0], particle2.force[1],
                 particle2.force[2]);
  }

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
      const std::shared_ptr<config::Config> &config)
      : particles(std::move(particles)),
        physics(physics),
        plotter(std::move(plotter)),
        config(config),
        end_time(config->end_time),
        delta_t(config->delta_t),
        thermostat(Tinit, Ttarget, deltaT, nthermostat, config->seed){};

  /*! <p> Function for position calculation </p>
   *
   * calculates the position for all particles, takes no arguments and has no return value
   */
  auto calculate_position() -> void {
    SPDLOG_DEBUG("Updating positions");

    particles.linear([this](auto &p) { calculate_position_particle(p); });
  }

  /*! <p> Function for velocity calculation </p>
  * calculates the velocity for all particles, takes no arguments and has no return value
  */
  auto calculate_velocity() -> void {
    SPDLOG_DEBUG("Updating velocities");
    particles.linear([this](auto &p) { calculate_velocity_particle(p); });
  }
  /**
  * @brief Calculates forces between particles.
  *
  * Resets forces for all particles, then calculates and updates forces for each particle pair.
  */

  auto calculate_force() -> void {
    SPDLOG_DEBUG("Starting force calculation");
    particles.linear([this](auto &p) { calculate_old_force_particle(p); });
    particles.pairwise([this](auto p) { calculate_force_particle_pair(p); });
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

    thermostat.initializeVelocities(particles, use_brownian_motion, config->brownian_motion);
    calculate_force();

    // Plot initial position and forces
    if (IO) {
      plotter->plotParticles(particles, iteration);
      SPDLOG_DEBUG("Iteration {} plotted.", iteration);
    }

    while (current_time < end_time) {
      calculate_position();
      particles.boundary([this](auto p) { calculate_force_particle_pair(p); });
      particles.refresh();
      calculate_force();
      calculate_velocity();
      //thermostat.apply(particles); //TODO activate

      iteration++;
      if (IO && iteration % interval == 0) {
        plotter->plotParticles(particles, iteration);
        SPDLOG_DEBUG("Iteration {} plotted.", iteration);
      }

      SPDLOG_DEBUG("Iteration {} finished.", iteration);

      current_time += delta_t;
    }
    spdlog::info("Output written. Terminating...");
  }
};

}// namespace simulator