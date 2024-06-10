//
// Created by Julius on 10.06.2024.
//
#ifndef THERMOSTAT_H
#define THERMOSTAT_H

#include "Particle.h"
#include "container/container.h"
#include <bits/random.h>
#include <cmath>
#include <limits>

class Thermostat {
 public:
  /**
     * @brief Constructor for the Thermostat class.
     * @param Tinit Initial temperature of the system.
     * @param Ttarget Target temperature to maintain. Defaults to Tinit if not provided.
     * @param deltaT Maximum temperature change allowed per application.
     * @param nthermostat Number of steps after which the thermostat is applied.
     * @param seed Seed for the random number generator.
     */
  Thermostat(double Tinit, double Ttarget = -1, double deltaT = std::numeric_limits<double>::infinity(), int nthermostat = 1, unsigned int seed = std::random_device{}());

  /**
     * @brief Applies the thermostat to adjust particle velocities.
     * @param particles Container of particles to apply the thermostat to.
     */
  void apply(container::particle_container &particles);

  /**
    * @brief Initializes particle velocities, optionally using Brownian motion.
    * @param particles Container of particles to initialize.
    * @param useBrownianMotion Whether to initialize using Brownian motion.
    * @param brownianMotion Specifies the Brownian motion constant.
    */
  void initializeVelocities(container::particle_container &particles, bool useBrownianMotion, double brownianMotion);

 private:
  /**
     * @brief Calculates the current temperature of the system.
     * @param particles Container of particles.
     * @return The current temperature.
     */
  double calculateCurrentTemperature(container::particle_container &particles);

  /**
     * @brief Calculates the total kinetic energy of the system.
     * @param particles Container of particles.
     * @return The total kinetic energy.
     */
  double calculateKineticEnergy(container::particle_container &particles);

  /**
     * @brief Scales the velocities of all particles.
     * @param particles Container of particles.
     * @param scalingFactor The factor by which to scale the velocities.
     */
  void scaleVelocities(container::particle_container &particles, double scalingFactor);

  double Tinit;
  double Ttarget;
  double deltaT;
  int nthermostat;
  int stepCount;
  unsigned int seed;
};

#endif// THERMOSTAT_H
