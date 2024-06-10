//
// Created by Julius on 10.06.2024.
//
#include "Particle.h"
#include "simulator/physics/Thermostat.h"
#include "container/container.h"
#include "utils/ArrayUtils.h"
#include <array>
#include <vector>
#include <gtest/gtest.h>
#include <memory>  // Required for std::construct_at

class ThermostatTest : public ::testing::Test {
protected:
    container::particle_container particles;
    Thermostat thermostat;
    const double tolerance = 0.0001;

    void SetUp() override {
        std::vector<Particle> particle_vector;

        // Initialize particles
        for (int i = 0; i < 10; ++i) {
            std::array<double, 3> coordinates = {static_cast<double>(i), static_cast<double>(i), static_cast<double>(i)};
            std::array<double, 3> velocity = {static_cast<double>(i + 1), static_cast<double>(i + 1), static_cast<double>(i + 1)};
            particle_vector.emplace_back(coordinates, velocity, 1.0, 0);
        }
        particles = container::particle_container(std::move(particle_vector));
    }

    void initializeThermostat(double Tinit, double Ttarget, double deltaT, int nthermostat, unsigned int seed) {
        thermostat = Thermostat(Tinit, Ttarget, deltaT, nthermostat, seed);
    }

    double calculateTemperature() {
        return thermostat.calculateCurrentTemperature(particles);
    }
};

/**
 * Test: test_heating
 *
 * Verifies that the thermostat correctly heats the particle system towards the target temperature.
 */
TEST_F(ThermostatTest, test_heating) {
    initializeThermostat(100.0, 200.0, 10.0, 1, 12345);
    thermostat.apply(particles);
    double newTemperature = calculateTemperature();

    EXPECT_GT(newTemperature, 100.0);
    EXPECT_LT(newTemperature, 110.0);
}

/**
 * Test: test_cooling
 *
 * Verifies that the thermostat correctly cools the particle system towards the target temperature.
 */
TEST_F(ThermostatTest, test_cooling) {
    initializeThermostat(200.0, 100.0, 10.0, 1, 12345);
    thermostat.apply(particles);
    double newTemperature = calculateTemperature();

    EXPECT_LT(newTemperature, 200.0);
    EXPECT_GT(newTemperature, 190.0);
}

/**
 * Test: test_holding_temperature
 *
 * Verifies that the thermostat holds the temperature constant when the current temperature is equal to the target temperature.
 */
TEST_F(ThermostatTest, test_holding_temperature) {
    initializeThermostat(150.0, 150.0, 10.0, 1, 12345);
    double initialTemperature = calculateTemperature();
    thermostat.apply(particles);
    double newTemperature = calculateTemperature();

    EXPECT_NEAR(initialTemperature, newTemperature, tolerance);
}