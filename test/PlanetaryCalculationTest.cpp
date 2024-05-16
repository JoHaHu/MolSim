//
// Created by TimSc on 15.05.2024.
//

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-avoid-magic-numbers"

#include "lib/Particle.h"
#include "lib/simulator/physics/Gravity.h"
#include "lib/utils/ArrayUtils.h"
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

TEST_F(PlanetaryCalculationTest, force_calculation_same_particles) {
  SetUp();

  auto calculatedForce = simulator::physics::Gravity::calculate_force(particleA, particleA);
  std::array<double, 3> zeroForce = {0.0, 0.0, 0.0};

  EXPECT_EQ(calculatedForce, zeroForce);
}

TEST_F(PlanetaryCalculationTest, force_calculation_simple_norm) {
  SetUp();

  // force between particleA and particleB
  auto calculatedForce = simulator::physics::Gravity::calculate_force(particleA, particleB);

  // calculating force by hand and checking if the return is the expected result
  // for each axis (x, y and z are meant as axes here)

  auto xDifference = 1.0;
  auto yDifference = 1.0;
  auto zDifference = 2.0;

  const Container auto positionDifference = std::array<double, 3>({xDifference, yDifference, zDifference});

  // L2-norm = sqrt(1² + 1² + 2²) = 2.449489743
  // multiplied mass of particles = 1.0 * 4.0 = 4.0
  // (L2-norm)³ = 14.69693846

  const Container auto actualForce = 4.0 / 14.69693846 * positionDifference;

  // check if each axis force is exact enough (based on the manually computed numbers)
  // in this case 7 digits accuracy is enough
  EXPECT_TRUE(calculatedForce[0] - actualForce[0] < 0.000001);
  EXPECT_TRUE(calculatedForce[1] - actualForce[1] < 0.000001);
  EXPECT_TRUE(calculatedForce[2] - actualForce[2] < 0.000001);
}


TEST_F(PlanetaryCalculationTest, force_calculation_edge_norm) {
  SetUp();

  auto calculatedForce = simulator::physics::Gravity::calculate_force(particleA, particleD);

  // calculating value differences for each axis by hand (x, y and z are meant as axes here)

  // xDifference = 219.4
  // yDifference = 321.06
  // zDifference = 0.45

  auto xDifference = 219.4;
  auto yDifference = 321.06;
  auto zDifference = 0.45;

  const auto positionDifference = std::array<double, 3>({xDifference, yDifference, zDifference});

  // computing the L2 / Euclidean norm manually
  auto l2Norm = sqrt((pow(xDifference, 2) + pow(yDifference, 2) + pow(zDifference, 2)));

  // L2-norm = sqrt(219.4² + 321.06² + 0.45²) = 388.8651258499597
  // multiplied mass of particles = 1.0 * 1005.34 = 1005.34
  // (L2-norm)³ = 58802662.352711186

  const Container auto actualForce = 1005.34 / 58802662.352711186 * positionDifference;

  // check if each axis force is exact enough (based on the manually computed numbers)
  // in this case 7 digits accuracy is enough
  EXPECT_TRUE(calculatedForce[0] - actualForce[0] < 0.000001);
  EXPECT_TRUE(calculatedForce[1] - actualForce[1] < 0.000001);
  EXPECT_TRUE(calculatedForce[2] - actualForce[2] < 0.000001);
}



auto main(int argc, char **argv) -> int {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
#pragma clang diagnostic pop
