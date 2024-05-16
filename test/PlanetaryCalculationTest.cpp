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
  Particle particleA, particleB, particleC;

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

    // construct particles with different parameter values
    particleA = Particle(coordinatesA, velocityA, 1.0, 0);
    particleB = Particle(coordinatesB, velocityB, 4.0, 0);
    particleC = Particle(coordinatesC, velocityC, 10.0, 0);
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

  auto calculatedForce = simulator::physics::Gravity::calculate_force(particleA, particleB);

  // calculating value differences for each axis (x, y and z are meant as axes here
  // despite attribute x saving all coordinate values
  auto xDifference = particleB.position[0] - particleA.position[0];
  auto yDifference = particleB.position[1] - particleA.position[1];
  auto zDifference = particleB.position[2] - particleA.position[2];

  const auto positionDifference = std::array<double, 3>({xDifference, yDifference, zDifference});

  // computing the L2 / Euclidean norm manually
  auto l2Norm = sqrt((pow(xDifference, 2) + pow(yDifference, 2) + pow(zDifference, 2)));

  // TODO: fix multiply operation
  // const Container auto actualForce = (particleA.m * particleB.m) / pow(l2Norm, 3) * positionDifference;

  //EXPECT_EQ(calculatedForce, actualForce);
}

auto main(int argc, char **argv) -> int {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
#pragma clang diagnostic pop
