
#include "simulator/physics/Thermostat.h"
#include "Particle.h"
#include "container/container.h"
#include <array>
#include <gtest/gtest.h>
#include <memory>// Required for std::construct_at

/**
class ThermostatTest : public ::testing::Test {
 protected:
  container::ParticleContainer<3> particles;
  Thermostat<3> thermostat = Thermostat<3>();
  const double tolerance = 0.0001;

  // Constants for test usage
  const double t_init = 100.0;
  const double t_target = 200.0;
  const double delta_t = 10.0;
  const unsigned int seed = 12345;

  void SetUp() override {
    // Initialize particles
    for (int i = 0; i < 10; ++i) {
      std::array<double, 3> coordinates = {static_cast<double>(i), static_cast<double>(i), static_cast<double>(i)};
      std::array<double, 3> velocity = {static_cast<double>(i + 1), static_cast<double>(i + 1), static_cast<double>(i + 1)};
      particles.insert(Particle<3>(coordinates, velocity, 1.0, 1));
    }
  }

  void initializeThermostat(double t_init, std::optional<double> t_target, std::optional<double> delta_t, unsigned int seed) {
    thermostat = Thermostat<3>(t_init, t_target, delta_t, seed);
  }

  double calculateTemperature() {
    return thermostat.calculateCurrentTemperature(particles);
  }
};
**/

/**
 * Test: test_heating
 *
 * Verifies that the thermostat correctly heats the particle system towards the target temperature.

TEST_F(ThermostatTest, test_heating) {

  initializeThermostat(t_init, t_target, delta_t, seed);
  thermostat.apply(particles);
  double newTemperature = calculateTemperature();

  EXPECT_GT(newTemperature, 100.0);
  EXPECT_LT(newTemperature, 110.0);
}
 */

/**
 * Test: test_cooling
 *
 * Verifies that the thermostat correctly cools the particle system towards the target temperature.

TEST_F(ThermostatTest, test_cooling) {
  initializeThermostat(t_init, t_target, delta_t, seed);
  thermostat.apply(particles);
  double newTemperature = calculateTemperature();

  EXPECT_LT(newTemperature, 200.0);
  EXPECT_GT(newTemperature, 190.0);
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
  double initialTemperature = calculateTemperature();
  thermostat.apply(particles);
  double newTemperature = calculateTemperature();

  EXPECT_NEAR(initialTemperature, newTemperature, tolerance);
}

*/