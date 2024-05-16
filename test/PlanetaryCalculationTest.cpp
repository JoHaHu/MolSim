#include "Particle.h"
#include "simulator/physics/Gravity.h"
#include "utils/ArrayUtils.h"
#include <array>
#include <cmath>
#include <gtest/gtest.h>

class PlanetaryCalculationTest : public ::testing::Test {
 public:
  Particle particleA, particleB, particleC, particleD;

  void SetUp() override {

    // Particle A
    std::array<double, 3> coordinatesA = {0.0, 0.0, 0.0};
    std::array<double, 3> velocityA = {2.0, 3.0, 3.0};

    // Particle B
    std::array<double, 3> coordinatesB = {1.0, 1.0, 2.0};
    std::array<double, 3> velocityB = {3.0, 2.0, 3.0};

    // Particle C
    std::array<double, 3> coordinatesC = {2.0, 3.0, 4.0};
    std::array<double, 3> velocityC = {5.0, 6.0, 3.0};

    // Particle D
    std::array<double, 3> coordinatesD = {219.4, 321.06, 0.45};
    std::array<double, 3> velocityD = {523.0, 62.2, 303.67};

    // construct particles with different parameter values
    particleA = Particle(coordinatesA, velocityA, 1.0, 0);
    particleB = Particle(coordinatesB, velocityB, 4.0, 0);
    particleC = Particle(coordinatesC, velocityC, 10.0, 0);
    particleD = Particle(coordinatesD, velocityD, 1005.34, 0);
  }
};

/**
 * Test: force_calculation_same_particles
 *
 * Verifies that the gravitational force calculation returns zero when the same particle is compared to itself.
 *
 * Ensures that the calculated force is {0.0, 0.0, 0.0} when using the same particle for both inputs.
 */
TEST_F(PlanetaryCalculationTest, force_calculation_same_particles) {
  SetUp();

  auto calculated_force = simulator::physics::Gravity::calculate_force(particleA, particleA);
  std::array<double, 3> zero_force = {0.0, 0.0, 0.0};

  EXPECT_EQ(calculated_force, zero_force);
}

/**
 * Test: force_calculation_simple_norm
 *
 * Verifies the gravitational force calculation between two particles using manual calculation for comparison.
 *
 * Ensures the calculated force vector is accurate to within seven decimal places for each axis (x, y, z).
 */
TEST_F(PlanetaryCalculationTest, force_calculation_simple_norm) {
  SetUp();

  // force between particleA and particleB
  auto calculated_force = simulator::physics::Gravity::calculate_force(particleA, particleB);

  // calculating force by hand and checking if the return is the expected result
  // for each axis (x, y and z are meant as axes here)

  auto x_difference = 1.0;
  auto y_difference = 1.0;
  auto z_difference = 2.0;

  const Container auto position_difference = std::array<double, 3>({x_difference, y_difference, z_difference});

  // L2-norm = sqrt(1² + 1² + 2²) = 2.449489743
  // multiplied mass of particles = 1.0 * 4.0 = 4.0
  // (L2-norm)³ = 14.69693846

  const Container auto actual_force = 4.0 / 14.69693846 * position_difference;

  // check if each axis force is exact enough (based on the manually computed numbers)
  // in this case 7 digits accuracy is enough
  EXPECT_TRUE(calculated_force[0] - actual_force[0] < 0.000001);
  EXPECT_TRUE(calculated_force[1] - actual_force[1] < 0.000001);
  EXPECT_TRUE(calculated_force[2] - actual_force[2] < 0.000001);
}

/**
 * Test: force_calculation_edge_norm
 *
 * Verifies the gravitational force calculation between two particles with a large positional difference.
 *
 * Ensures the calculated force vector is accurate to within seven decimal places for each axis (x, y, z) based on manual calculations.
 */
TEST_F(PlanetaryCalculationTest, force_calculation_edge_norm) {
  SetUp();

  auto calculated_force = simulator::physics::Gravity::calculate_force(particleA, particleD);

  // calculating value differences for each axis by hand (x, y and z are meant as axes here)

  // x_difference = 219.4
  // y_difference = 321.06
  // z_difference = 0.45

  auto x_difference = 219.4;
  auto y_difference = 321.06;
  auto z_difference = 0.45;

  const auto position_difference = std::array<double, 3>({x_difference, y_difference, z_difference});

  // computing the L2 / Euclidean norm manually
  // L2-norm = sqrt(219.4² + 321.06² + 0.45²) = 388.8651258499597
  // multiplied mass of particles = 1.0 * 1005.34 = 1005.34
  // (L2-norm)³ = 58802662.352711186

  const Container auto actual_force = 1005.34 / 58802662.352711186 * position_difference;

  // check if each axis force is exact enough (based on the manually computed numbers)
  // in this case 7 digits accuracy is enough
  EXPECT_TRUE(calculated_force[0] - actual_force[0] < 0.000001);
  EXPECT_TRUE(calculated_force[1] - actual_force[1] < 0.000001);
  EXPECT_TRUE(calculated_force[2] - actual_force[2] < 0.000001);
}
