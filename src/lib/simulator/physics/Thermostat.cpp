//
// Created by Julius on 10.06.2024.
//
#include "Thermostat.h"
#include <cmath>
#include <limits>
#include "spdlog/spdlog.h"
#include "utils/MaxwellBoltzmannDistribution.h"

Thermostat::Thermostat(double Tinit, double Ttarget, double deltaT, int nthermostat, unsigned int seed)
    : Tinit(Tinit), Ttarget(Ttarget), deltaT(deltaT), nthermostat(nthermostat), stepCount(0),
      seed(seed) {
    if (Ttarget == -1) {
        this->Ttarget = Tinit;
    }
    SPDLOG_DEBUG("Thermostat initialized with Tinit={}, Ttarget={}, deltaT={}, nthermostat={}, seed={}", Tinit, Ttarget, deltaT, nthermostat, seed);
}

void Thermostat::apply(container::particle_container &particles) {
    if (stepCount % nthermostat == 0) {
        double currentTemperature = calculateCurrentTemperature(particles);
        SPDLOG_DEBUG("Current temperature calculated: {}", currentTemperature);

        double tempDifference = Ttarget - currentTemperature;
        if (std::abs(tempDifference) > deltaT) {
            tempDifference = (tempDifference > 0) ? deltaT : -deltaT;
        }

        double newTemperature = currentTemperature + tempDifference;
        double scalingFactor = std::sqrt(newTemperature / currentTemperature);
        SPDLOG_TRACE("Scaling velocities with factor: {}", scalingFactor);
        scaleVelocities(particles, scalingFactor);
    }
    stepCount++;
    SPDLOG_TRACE("Thermostat applied at step {}", stepCount);
}

void Thermostat::initializeVelocities(container::particle_container &particles, bool useBrownianMotion, double brownianMotion) {
    if (useBrownianMotion) {
        SPDLOG_INFO("Initializing velocities with Brownian motion");
        particles.linear([this, brownianMotion](Particle &particle) {
            std::array<double, 3> velocity = maxwellBoltzmannDistributedVelocity(brownianMotion, 3, seed);
            particle.velocity[0] = velocity[0];
            particle.velocity[1] = velocity[1];
            particle.velocity[2] = velocity[2];
            SPDLOG_TRACE("Particle initialized: position = ({}, {}, {}), velocity = ({}, {}, {}), mass = {}",
                         particle.position[0], particle.position[1], particle.position[2],
                         particle.velocity[0], particle.velocity[1], particle.velocity[2], particle.mass);
        });
    }
}

double Thermostat::calculateCurrentTemperature(container::particle_container &particles) {
    double kineticEnergy = calculateKineticEnergy(particles);
    double currentTemperature = (2.0 * kineticEnergy) / (particles.size() * 3);
    SPDLOG_TRACE("Current temperature calculated: {}", currentTemperature);
    return currentTemperature;
}

double Thermostat::calculateKineticEnergy(container::particle_container &particles) {
    double Ekin = 0;
    particles.linear([&Ekin](Particle &particle) {
        Ekin += 0.5 * particle.mass * (particle.velocity[0] * particle.velocity[0] +
                                       particle.velocity[1] * particle.velocity[1] +
                                       particle.velocity[2] * particle.velocity[2]);
        SPDLOG_TRACE("Particle kinetic energy contribution: {}", Ekin);
    });
    SPDLOG_DEBUG("Total kinetic energy calculated: {}", Ekin);
    return Ekin;
}

void Thermostat::scaleVelocities(container::particle_container &particles, double scalingFactor) {
    particles.linear([scalingFactor](Particle &particle) {
        particle.velocity[0] *= scalingFactor;
        particle.velocity[1] *= scalingFactor;
        particle.velocity[2] *= scalingFactor;
        SPDLOG_TRACE("Scaled particle velocity: ({}, {}, {})", particle.velocity[0], particle.velocity[1], particle.velocity[2]);
    });
}
