#include "Particle.h"
#include "simulator/physics/Gravity.h"
#include "utils/ArrayUtils.h"
#include <array>
#include <gtest/gtest.h>

class PlanetaryCalculationTest : public ::testing::Test {
 protected:
  Particle<3> particleA, particleB, particleC, particleD, particleE, particleF;
  simulator::physics::Gravity gravity;

  PlanetaryCalculationTest()
      : particleA({{0.0, 0.0, 0.0}}, {{2.0, 3.0, 3.0}}, 1.0, 0),
        particleB({{1.0, 1.0, 2.0}}, {{3.0, 2.0, 3.0}}, 4.0, 0),
        particleC({{2.0, 3.0, 4.0}}, {{5.0, 6.0, 3.0}}, 10.0, 0),
        particleD({{219.4, 321.06, 0.45}}, {{523.0, 62.2, 303.67}}, 1005.34, 0),
        particleE({{0.1, 0.1, 0.1}}, {{1.0, 1.0, 1.0}}, 2.0, 0),
        particleF({{0.0, 0.0, 10.0}}, {{0.0, 0.0, -1.0}}, 5.0, 0) {}
};

/**
 * Test: force_calculation_simple_norm
 *
 * Verifies the gravitational force calculation between two particles using manual calculation for comparison.
 *
 * Ensures the calculated force vector is accurate to within seven decimal places for each axis (x, y, z).
 */
TEST_F(PlanetaryCalculationTest, force_calculation_simple_norm) {
  std::array<double, 3> result_force;
  std::array<double, 3> correction = {0.0, 0.0, 0.0};
  long type = 0;

  gravity.calculateForce_3D(particleA.position[0], particleA.position[1], particleA.position[2], particleA.mass, type,
                            particleB.position[0], particleB.position[1], particleB.position[2], particleB.mass, type,
                            result_force[0], result_force[1], result_force[2], correction);

  auto x_difference = 1.0;
  auto y_difference = 1.0;
  auto z_difference = 2.0;

  const auto position_difference = std::array<double, 3>({x_difference, y_difference, z_difference});

  const auto norm = std::sqrt(x_difference * x_difference + y_difference * y_difference + z_difference * z_difference);
  const auto temp = (particleA.mass * particleB.mass) / (norm * norm * norm);
  const auto actual_force = std::array<double, 3>({temp * x_difference, temp * y_difference, temp * z_difference});

  EXPECT_NEAR(result_force[0], actual_force[0], 1e-6);
  EXPECT_NEAR(result_force[1], actual_force[1], 1e-6);
  EXPECT_NEAR(result_force[2], actual_force[2], 1e-6);
}

/**
 * Test: force_calculation_edge_norm
 *
 * Verifies the gravitational force calculation between two particles with a large positional difference.
 *
 * Ensures the calculated force vector is accurate to within seven decimal places for each axis (x, y, z) based on manual calculations.
 */
TEST_F(PlanetaryCalculationTest, force_calculation_edge_norm) {
  std::array<double, 3> result_force;
  std::array<double, 3> correction = {0.0, 0.0, 0.0};
  long type = 0;

  gravity.calculateForce_3D(particleA.position[0], particleA.position[1], particleA.position[2], particleA.mass, type,
                            particleD.position[0], particleD.position[1], particleD.position[2], particleD.mass, type,
                            result_force[0], result_force[1], result_force[2], correction);

  auto x_difference = 219.4;
  auto y_difference = 321.06;
  auto z_difference = 0.45;

  const auto position_difference = std::array<double, 3>({x_difference, y_difference, z_difference});

  const auto norm = std::sqrt(x_difference * x_difference + y_difference * y_difference + z_difference * z_difference);
  const auto temp = (particleA.mass * particleD.mass) / (norm * norm * norm);
  const auto actual_force = std::array<double, 3>({temp * x_difference, temp * y_difference, temp * z_difference});

  EXPECT_NEAR(result_force[0], actual_force[0], 1e-6);
  EXPECT_NEAR(result_force[1], actual_force[1], 1e-6);
  EXPECT_NEAR(result_force[2], actual_force[2], 1e-6);
}

/**
 * Test: force_calculation_small_difference
 *
 * Verifies the gravitational force calculation between two particles with very small positional differences.
 *
 * Ensures the calculated force vector is accurate to within seven decimal places for each axis (x, y, z).
 */
TEST_F(PlanetaryCalculationTest, force_calculation_small_difference) {
  std::array<double, 3> result_force;
  std::array<double, 3> correction = {0.0, 0.0, 0.0};
  long type = 0;

  gravity.calculateForce_3D(particleA.position[0], particleA.position[1], particleA.position[2], particleA.mass, type,
                            particleE.position[0], particleE.position[1], particleE.position[2], particleE.mass, type,
                            result_force[0], result_force[1], result_force[2], correction);

  auto x_difference = 0.1;
  auto y_difference = 0.1;
  auto z_difference = 0.1;

  const auto position_difference = std::array<double, 3>({x_difference, y_difference, z_difference});

  const auto norm = std::sqrt(x_difference * x_difference + y_difference * y_difference + z_difference * z_difference);
  const auto temp = (particleA.mass * particleE.mass) / (norm * norm * norm);
  const auto actual_force = std::array<double, 3>({temp * x_difference, temp * y_difference, temp * z_difference});

  EXPECT_NEAR(result_force[0], actual_force[0], 1e-6);
  EXPECT_NEAR(result_force[1], actual_force[1], 1e-6);
  EXPECT_NEAR(result_force[2], actual_force[2], 1e-6);
}

/**
 * Test: force_calculation_along_one_axis
 *
 * Verifies the gravitational force calculation between two particles aligned along one axis (z-axis).
 *
 * Ensures the calculated force vector is accurate to within seven decimal places for each axis (x, y, z).
 */
TEST_F(PlanetaryCalculationTest, force_calculation_along_one_axis) {
  std::array<double, 3> result_force;
  std::array<double, 3> correction = {0.0, 0.0, 0.0};
  long type = 0;

  gravity.calculateForce_3D(particleA.position[0], particleA.position[1], particleA.position[2], particleA.mass, type,
                            particleF.position[0], particleF.position[1], particleF.position[2], particleF.mass, type,
                            result_force[0], result_force[1], result_force[2], correction);

  auto x_difference = 0.0;
  auto y_difference = 0.0;
  auto z_difference = 10.0;

  const auto position_difference = std::array<double, 3>({x_difference, y_difference, z_difference});

  const auto norm = std::sqrt(x_difference * x_difference + y_difference * y_difference + z_difference * z_difference);
  const auto temp = (particleA.mass * particleF.mass) / (norm * norm * norm);
  const auto actual_force = std::array<double, 3>({temp * x_difference, temp * y_difference, temp * z_difference});

  EXPECT_NEAR(result_force[0], actual_force[0], 1e-6);
  EXPECT_NEAR(result_force[1], actual_force[1], 1e-6);
  EXPECT_NEAR(result_force[2], actual_force[2], 1e-6);
}
