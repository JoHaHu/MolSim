#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-avoid-magic-numbers"

#include "lib/Particle.h"
#include "lib/simulator/physics/Gravity.h"
#include <gtest/gtest.h>
#include <array>

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

TEST_F(PlanetaryCalculationTest, ForceCalculation_SameParticles) {
  SetUp();

  auto force = Gravity::calculateF(particleA, particleA);
  std::array<double, 3> zeroForce = {0.0, 0.0, 0.0};

  EXPECT_EQ(force, zeroForce);
}



auto main(int argc, char **argv) -> int {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
#pragma clang diagnostic pop
