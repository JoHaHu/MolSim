#pragma once

#include "Particle.h"
#include "config/Config.h"
#include "container/Container.h"
#include "experimental/simd"
#include "simulator/io/Plotter.h"
#include "simulator/physics/ForceModel.h"
#include "simulator/physics/Gravity.h"
#include "simulator/physics/LennardJones.h"
#include "simulator/physics/Thermostat.h"

#include "simulator/io/checkpoint/Checkpointer.h"
#include "utils/ArrayUtils.h"
#include <cmath>
#include <omp.h>
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
  std::unique_ptr<container::Container<DIMENSIONS>> particles;
  std::unique_ptr<io::Plotter<DIMENSIONS>> plotter;
  std::shared_ptr<config::Config> config;
  Thermostat<DIMENSIONS> thermostat;
  Checkpointer<DIMENSIONS> checkpoint;

        double end_time;
        double delta_t;
        unsigned long iteration = 0;
        double gravity = 0.0;

 public:
  /**
   * \param particles the particle container
   * \param plotter An instance plotter.
   * \param config the runtime configuration
   * */
  explicit Simulator(std::unique_ptr<container::Container<DIMENSIONS>> &&particles,
                     std::unique_ptr<io::Plotter<DIMENSIONS>> &&plotter,
                     const std::shared_ptr<config::Config> &config,
                     Checkpointer<DIMENSIONS> checkpoint)
      : particles(std::move(particles)),
        plotter(std::move(plotter)),
        config(config),
        thermostat(config->temp_init, config->temp_target, config->max_temp_diff, config->seed),
        checkpoint(checkpoint),
        end_time(config->end_time),
        delta_t(config->delta_t),
        gravity(config->ljf_gravity){};

        auto inline calculate_position_particle(Particles<DIMENSIONS> &p, size_t index) const {
            const auto temp = pow(delta_t, 2) * (1 / (2 * p.mass[index]));
            if (p.fixed[index] == 0)[[likely]] {
                for (size_t i = 0; i < DIMENSIONS; ++i) {
                    p.positions[i][index] += delta_t * p.velocities[i][index] + temp * p.old_forces[i][index];
                }
            }
        }

        auto inline calculate_velocity_particle(Particles<DIMENSIONS> &p, size_t index) const {

            const auto temp = delta_t * (1 / (2 * p.mass[index]));

            for (size_t i = 0; i < DIMENSIONS; ++i) {
                p.velocities[i][index] += temp * (p.old_forces[i][index] + p.forces[i][index]);
            }
        }

  /*! <p> Function for position calculation </p>
   *
   * calculates the position for all particles, takes no arguments and has no return value
   */
  auto calculate_position() -> void {
    SPDLOG_DEBUG("Updating positions");
    particles->linear([this](Particles<DIMENSIONS> &p, size_t index) {
      calculate_position_particle(p, index);
    });
  }

        /*! <p> Function for velocity calculation </p>
        * calculates the velocity for all particles, takes no arguments and has no return value
        */
        auto calculate_velocity() -> void {
            SPDLOG_DEBUG("Updating velocities");

    particles->linear([this](Particles<DIMENSIONS> &p, size_t index) {
      calculate_velocity_particle(p, index);
    });
  }
  auto apply_gravity() -> void {
    SPDLOG_DEBUG("Applying gravity");

    particles->linear([this](Particles<DIMENSIONS> &p, size_t index) {
      p.forces[1][index] += simulator::physics::constant_gravity(p.mass[index], gravity);
    });
  }

        /**
        * @brief Calculates forces between particles.
        *
        * Resets forces for all particles, then calculates and updates forces for each particle pair.
        */

  auto calculate_force() -> void {
    SPDLOG_DEBUG("Starting force calculation");
    particles->boundary();
    particles->refresh();
    particles->pairwise(config->parallelized, config->vectorized);

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
    auto temp_interval = config->thermo_step;
    thermostat.initializeVelocities(*particles, config->use_brownian_motion, config->brownian_motion);
    calculate_force();
    // calculate twice to initialize old_force and force to have proper simulation and output
    particles->swap_force();
    calculate_force();

    // Plot initial position and forces
    if (IO) {
      plotter->plotParticles(*particles, iteration);
      SPDLOG_DEBUG("Iteration {} plotted.", iteration);
    }

    while (current_time < end_time) {
      calculate_position();
      particles->swap_force();
      calculate_force();

      if (gravity != 0) {
        apply_gravity();
      }
      if (temp_interval != 0 && iteration % temp_interval == 0) {
        thermostat.apply(*particles);
      }
      calculate_velocity();

      iteration++;
      if (IO && iteration % interval == 0) {
        plotter->plotParticles(*particles, iteration);
        SPDLOG_DEBUG("Iteration {} plotted.", iteration);
      }

                SPDLOG_DEBUG("Iteration {} finished.", iteration);

                current_time += delta_t;
            }

    if (config->output_checkpoint.has_value()) {
        checkpoint.save_checkpoint(*config->output_checkpoint, particles);
    }

            SPDLOG_INFO("Output written. Terminating...");
        }
    };

}// namespace simulator