// #include "simulator/physics/Thermostat.h"
// #include "Particle.h"
// #include <array>
// #include <gtest/gtest.h>
// #include <memory> // Required for std::construct_at
// #include "container/ParticleContainer.h"
//
// class ThermostatTest : public ::testing::Test {
//  protected:
//   std::unique_ptr<ParticleContainer<2>> particles;
//   Thermostat<2> thermostat;
//   const double tolerance = 0.0001;
//
//   // Constants for test usage
//   const double t_init = 10.0;
//   const double t_target = 200.0;
//   const double delta_t = 10.0;
//   const unsigned int seed = 12345;
//
//   void SetUp() override {
//     particles = std::make_unique<ParticleContainer<2>>();
//     thermostat = Thermostat<2>();
//
//     // Initialize particles
//     for (int i = 0; i < 100; ++i) {
//       std::array<double, 2> coordinates = {static_cast<double>(i), static_cast<double>(i)};
//       particles->insert(Particle<2>(coordinates, {0, 0}, 1.0, 1));
//     }
//   }
//
//   void initializeThermostat(double t_init, std::optional<double> t_target, std::optional<double> delta_t, unsigned int seed) {
//     thermostat = Thermostat<2>(t_init, t_target, delta_t, seed);
//     thermostat.initializeVelocities(*particles, true, 0.0);
//   }
// };
//
// TEST_F(ThermostatTest, test_initialize) {
//   initializeThermostat(t_init, t_target, delta_t, seed);
//   double temp = thermostat.calculateCurrentTemperature(*particles);
//   EXPECT_NEAR(temp, t_init, 2);
// }
//
// /**
//  * Test: test_heating
//  *
//  * Verifies that the thermostat correctly heats the particle system towards the target temperature.
// */
// TEST_F(ThermostatTest, test_heating) {
//   initializeThermostat(t_init, t_target, delta_t, seed);
//   double old_temp = thermostat.calculateCurrentTemperature(*particles);
//   for (int i = 0; i < 4; ++i) {
//     thermostat.apply(*particles);
//     double newTemperature = thermostat.calculateCurrentTemperature(*particles);
//     EXPECT_GT(newTemperature, old_temp);
//     EXPECT_GT(newTemperature, 10);
//     old_temp = newTemperature;
//   }
// }
//
// /**
//  * Test: test_cooling
//  *
//  * Verifies that the thermostat correctly cools the particle system towards the target temperature.
// */
// TEST_F(ThermostatTest, test_cooling) {
//   initializeThermostat(t_target, t_init, delta_t, seed);
//   double old_temp = thermostat.calculateCurrentTemperature(*particles);
//   for (int i = 0; i < 4; ++i) {
//     thermostat.apply(*particles);
//     double newTemperature = thermostat.calculateCurrentTemperature(*particles);
//     EXPECT_LT(newTemperature, old_temp);
//     EXPECT_LT(newTemperature, 200);
//     old_temp = newTemperature;
//   }
// }
//
// /**
//  * Test: test_holding_temperature
//  *
//  * Verifies that the thermostat holds the temperature constant when the current temperature is equal to the target temperature.
// */
// TEST_F(ThermostatTest, test_holding_temperature) {
//   const double t_init = 200.0;
//   const double t_target = 200.0;
//   const double delta_t = 10.0;
//   const unsigned int seed = 12345;
//
//   initializeThermostat(t_init, t_target, delta_t, seed);
//   // Let temp settle
//   for (int i = 0; i < 10; ++i) {
//     thermostat.apply(*particles);
//   }
//   double initialTemperature = thermostat.calculateCurrentTemperature(*particles);
//   EXPECT_NEAR(initialTemperature, 200, 1);
//
//   for (int i = 0; i < 4; ++i) {
//     thermostat.apply(*particles);
//     double newTemperature = thermostat.calculateCurrentTemperature(*particles);
//     EXPECT_NEAR(newTemperature, 200, 1);
//   }
// }
