// #include "simulator/physics/Gravity.h"
// #include "Particle.h"
// #include "utils/ArrayUtils.h"
// #include <array>
// #include <gtest/gtest.h>
//
// class PlanetaryCalculationTest : public ::testing::Test {
//  protected:
//   Particle<3> particleA, particleB, particleC, particleD;
//   simulator::physics::Gravity gravity;
//
//   void SetUp() override {
//     // Particle A
//     std::array<double, 3> coordinatesA = {0.0, 0.0, 0.0};
//     std::array<double, 3> velocityA = {2.0, 3.0, 3.0};
//     particleA = Particle<3>(coordinatesA, velocityA, 1.0, 0);
//
//     // Particle B
//     std::array<double, 3> coordinatesB = {1.0, 1.0, 2.0};
//     std::array<double, 3> velocityB = {3.0, 2.0, 3.0};
//     particleB = Particle<3>(coordinatesB, velocityB, 4.0, 0);
//
//     // Particle C
//     std::array<double, 3> coordinatesC = {2.0, 3.0, 4.0};
//     std::array<double, 3> velocityC = {5.0, 6.0, 3.0};
//     particleC = Particle<3>(coordinatesC, velocityC, 10.0, 0);
//
//     // Particle D
//     std::array<double, 3> coordinatesD = {219.4, 321.06, 0.45};
//     std::array<double, 3> velocityD = {523.0, 62.2, 303.67};
//     particleD = Particle<3>(coordinatesD, velocityD, 1005.34, 0);
//   }
// };
//
// /**
//  * Test: force_calculation_same_particles
//  *
//  * Verifies that the gravitational force calculation returns zero when the same particle is compared to itself.
//  *
//  * Ensures that the calculated force is {0.0, 0.0, 0.0} when using the same particle for both inputs.
//  */
// TEST_F(PlanetaryCalculationTest, force_calculation_same_particles) {
//   std::array<double, 3> result_force;
//   std::array<double, 3> correction = {0.0, 0.0, 0.0};
//   long type = 0;
//
//   gravity.calculateForce_3D(particleA.position[0], particleA.position[1], particleA.position[2], particleA.mass, type,
//                             particleA.position[0], particleA.position[1], particleA.position[2], particleA.mass, type,
//                             result_force[0], result_force[1], result_force[2], correction);
//
//   std::array<double, 3> zero_force = {0.0, 0.0, 0.0};
//   EXPECT_EQ(result_force, zero_force);
// }
//
// /**
//  * Test: force_calculation_simple_norm
//  *
//  * Verifies the gravitational force calculation between two particles using manual calculation for comparison.
//  *
//  * Ensures the calculated force vector is accurate to within seven decimal places for each axis (x, y, z).
//  */
// TEST_F(PlanetaryCalculationTest, force_calculation_simple_norm) {
//   std::array<double, 3> result_force;
//   std::array<double, 3> correction = {0.0, 0.0, 0.0};
//   long type = 0;
//
//   gravity.calculateForce_3D(particleA.position[0], particleA.position[1], particleA.position[2], particleA.mass, type,
//                             particleB.position[0], particleB.position[1], particleB.position[2], particleB.mass, type,
//                             result_force[0], result_force[1], result_force[2], correction);
//
//   auto x_difference = 1.0;
//   auto y_difference = 1.0;
//   auto z_difference = 2.0;
//
//   const auto position_difference = std::array<double, 3>({x_difference, y_difference, z_difference});
//
//   const auto norm = std::sqrt(x_difference * x_difference + y_difference * y_difference + z_difference * z_difference);
//   const auto temp = (particleA.mass * particleB.mass) / (norm * norm * norm);
//   const auto actual_force = temp * position_difference;
//
//   EXPECT_NEAR(result_force[0], actual_force[0], 1e-6);
//   EXPECT_NEAR(result_force[1], actual_force[1], 1e-6);
//   EXPECT_NEAR(result_force[2], actual_force[2], 1e-6);
// }
//
// /**
//  * Test: force_calculation_edge_norm
//  *
//  * Verifies the gravitational force calculation between two particles with a large positional difference.
//  *
//  * Ensures the calculated force vector is accurate to within seven decimal places for each axis (x, y, z) based on manual calculations.
//  */
// TEST_F(PlanetaryCalculationTest, force_calculation_edge_norm) {
//   std::array<double, 3> result_force;
//   std::array<double, 3> correction = {0.0, 0.0, 0.0};
//   long type = 0;
//
//   gravity.calculateForce_3D(particleA.position[0], particleA.position[1], particleA.position[2], particleA.mass, type,
//                             particleD.position[0], particleD.position[1], particleD.position[2], particleD.mass, type,
//                             result_force[0], result_force[1], result_force[2], correction);
//
//   auto x_difference = 219.4;
//   auto y_difference = 321.06;
//   auto z_difference = 0.45;
//
//   const auto position_difference = std::array<double, 3>({x_difference, y_difference, z_difference});
//
//   const auto norm = std::sqrt(x_difference * x_difference + y_difference * y_difference + z_difference * z_difference);
//   const auto temp = (particleA.mass * particleD.mass) / (norm * norm * norm);
//   const auto actual_force = temp * position_difference;
//
//   EXPECT_NEAR(result_force[0], actual_force[0], 1e-6);
//   EXPECT_NEAR(result_force[1], actual_force[1], 1e-6);
//   EXPECT_NEAR(result_force[2], actual_force[2], 1e-6);
// }
