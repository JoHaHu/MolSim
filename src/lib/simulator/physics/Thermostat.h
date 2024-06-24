#pragma once

#include "Particle.h"
#include "container/container.h"
#include "range/v3/algorithm.hpp"
#include "range/v3/view/iota.hpp"
#include "spdlog/spdlog.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include <bits/random.h>
#include <cmath>
#include <limits>

template<const size_t DIMENSIONS>
class Thermostat {
 public:
  /**
     * @brief Constructor for the Thermostat class.
     * @param t_init Initial temperature of the system.
     * @param t_target Target temperature to maintain. Defaults to Tinit if not provided.
     * @param delta_t Maximum temperature change allowed per application.
     * @param nthermostat Number of steps after which the thermostat is applied.
     * @param seed Seed for the random number generator.
     */
  Thermostat(double t_init, std::optional<double> t_target, std::optional<double> delta_t, unsigned int seed = std::random_device{}())
      : t_init(t_init), seed(seed) {
    this->t_target = t_target.has_value() ? *t_target : this->t_init;
    this->delta_t = delta_t.has_value() ? *delta_t : std::numeric_limits<double>::infinity();

    SPDLOG_DEBUG("Thermostat initialized with Tinit={}, Ttarget={}, deltaT={},  seed={}", t_init, *t_target, *delta_t, seed);
  };

  Thermostat() = default;

  /**
     * @brief Applies the thermostat to adjust particle velocities.
     * @param particles Container of particles to apply the thermostat to.
     */
  void inline apply(container::ParticleContainer<DIMENSIONS> &particles) {
    double current_temperature = calculateCurrentTemperature(particles);
    SPDLOG_DEBUG("Current temperature calculated: {}", current_temperature);

    double temp_difference = t_target - current_temperature;
    if (std::abs(temp_difference) > delta_t) {
      temp_difference = std::signbit(temp_difference) ? -delta_t : delta_t;
    }

    double new_temperature = current_temperature + temp_difference;
    double scaling_factor = std::sqrt(new_temperature / current_temperature);
    SPDLOG_TRACE("Scaling velocities with factor: {}", scaling_factor);
    scaleVelocities(particles, scaling_factor);
  }

  /**
    * @brief Initializes particle velocities, optionally using Brownian motion.
    * @param particles Container of particles to initialize.
    * @param useBrownianMotion Whether to initialize using Brownian motion.
    * @param brownianMotion Specifies the Brownian motion constant.
    */
  void initializeVelocities(container::ParticleContainer<DIMENSIONS> &particles, bool useBrownianMotion, double brownianMotion) {
    if (useBrownianMotion) {
      SPDLOG_INFO("Initializing velocities with Brownian motion");

      particles.linear([this, brownianMotion](Particles<DIMENSIONS> &particles, size_t index) {
        double average_motion;
        if (brownianMotion == 0) {
          average_motion = std::sqrt(t_init / particles.mass[index]);
        } else {
          average_motion = brownianMotion;
        }
        std::array<double, DIMENSIONS> velocity = maxwellBoltzmannDistributedVelocity<DIMENSIONS>(average_motion, seed);
        for (size_t i = 0; i < DIMENSIONS; ++i) {
          particles.velocities[i][index] += velocity[i];
        }
      });
    }
  };

  /**
     * @brief Calculates the current temperature of the system.
     * @param particles Container of particles.
     * @return The current temperature.
     */
  auto calculateCurrentTemperature(container::ParticleContainer<DIMENSIONS> &particles) -> double {
    double kineticEnergy = calculateKineticEnergy(particles);
    double currentTemperature = (kineticEnergy) / (particles.size() * DIMENSIONS);
    SPDLOG_TRACE("Current temperature calculated: {}", currentTemperature);
    return currentTemperature;
  };

 private:
  /**
     * @brief Calculates the total kinetic energy of the system.
     * @param particles Container of particles.
     * @return The total kinetic energy.
     */
  auto calculateKineticEnergy(container::ParticleContainer<DIMENSIONS> &particles) -> double {
    double e_kin = 0;
    particles.linear([&](Particles<DIMENSIONS> &particles, size_t index) {
      e_kin += particles.mass[index] * ranges::fold_left(ranges::iota_view(0UL, DIMENSIONS), 0.0, [&](double acc, size_t i) {
                 return acc + particles.velocities[i][index] * particles.velocities[i][index];
               });
      SPDLOG_TRACE("Particle kinetic energy contribution: {}", e_kin);
    });
    SPDLOG_DEBUG("Total kinetic energy calculated: {}", e_kin);
    return e_kin;
  };

  /**
     * @brief Scales the velocities of all particles.
     * @param particles Container of particles.
     * @param scalingFactor The factor by which to scale the velocities.
     */
  void scaleVelocities(container::ParticleContainer<DIMENSIONS> &particles, double scalingFactor) {
    particles.linear([scalingFactor](Particles<DIMENSIONS> &particles, size_t index) {
      for (size_t i = 0; i < DIMENSIONS; ++i) {
        particles.velocities[i][index] *= scalingFactor;
      }
    });
  };

  double t_init;
  double t_target;
  double delta_t;
  unsigned int seed;
};
