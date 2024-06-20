
#include "simulator/physics/Thermostat.h"
#include "Particle.h"
#include "container/container.h"
#include <array>
#include <gtest/gtest.h>
#include <memory>// Required for std::construct_at

/**
class ThermostatTest : public ::testing::Test {
 protected:
  container::ParticleContainer<2> particles = container::ParticleContainer<2>(Particles<2>());
  Thermostat<2> thermostat = Thermostat<2>();
  const double tolerance = 0.0001;

  // Constants for test usage
  const double t_init = 10.0;
  const double t_target = 20.0;
  const double delta_t = 10.0;
  const unsigned int seed = 12345;

  void SetUp() override {
    // Initialize particles
    for (int i = 0; i < 3; ++i) {
      std::array<double, 2> coordinates = {static_cast<double>(i), static_cast<double>(i)};
      std::array<double, 2> velocity = {static_cast<double>(i + 1), static_cast<double>(i + 1)};
      particles.insert(Particle<2>(coordinates, velocity, 1.0, 1));
    }
  }

  void initializeThermostat(double t_init, std::optional<double> t_target, std::optional<double> delta_t, unsigned int seed) {
    thermostat = Thermostat<2>(t_init, t_target, delta_t, seed);
    thermostat.initializeVelocities(particles, true, std::sqrt(t_init));
  }
};
**/

TEST_F(ThermostatTest, test_initialize) {
  initializeThermostat(t_init, t_target, delta_t, seed);
  double temp = thermostat.calculateCurrentTemperature(particles);

  EXPECT_EQ(temp, t_init);
}

/**
 * Test: test_heating
 *
 * Verifies that the thermostat correctly heats the particle system towards the target temperature.

TEST_F(ThermostatTest, test_heating) {

  initializeThermostat(t_init, t_target, delta_t, seed);
  double old_temp = thermostat.calculateCurrentTemperature(particles);
  thermostat.apply(particles);
  double newTemperature = thermostat.calculateCurrentTemperature(particles);
  EXPECT_GT(newTemperature, old_temp);
  EXPECT_GT(newTemperature, 10);
}
 */

/**
 * Test: test_cooling
 *
 * Verifies that the thermostat correctly cools the particle system towards the target temperature.

TEST_F(ThermostatTest, test_cooling) {
  initializeThermostat(t_target, t_init, delta_t, seed);
  double old_temp = thermostat.calculateCurrentTemperature(particles);
  thermostat.apply(particles);
  double newTemperature = thermostat.calculateCurrentTemperature(particles);

  EXPECT_LT(newTemperature, old_temp);
  EXPECT_LT(newTemperature, 20);
}
*/

/**
 * Test: test_holding_temperature
 *
 * Verifies that the thermostat holds the temperature constant when the current temperature is equal to the target temperature.

TEST_F(ThermostatTest, test_holding_temperature) {
  const double t_init = 100.0;
  const double t_target = 200.0;
  const double delta_t = 10.0;
  const unsigned int seed = 12345;

  initializeThermostat(t_init, t_target, delta_t, seed);
  double initialTemperature = thermostat.calculateCurrentTemperature(particles);
  thermostat.apply(particles);
  double newTemperature = thermostat.calculateCurrentTemperature(particles);

  EXPECT_NEAR(initialTemperature, newTemperature, tolerance);
}

*/