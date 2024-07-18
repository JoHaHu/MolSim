#include "simulator/physics/ProfileCalculator.h"
#include "Particle.h"
#include "container/ParticleContainer.h"
#include <array>
#include <fstream>
#include <gtest/gtest.h>
#include <sstream>
#include <vector>

class ProfileCalculatorTest : public ::testing::Test {
 protected:
  std::unique_ptr<ParticleContainer<3>> particles;
  ProfileCalculator<3> profileCalculator;

  ProfileCalculatorTest()
      : profileCalculator(10) {}

  void SetUp() override {
    particles = std::make_unique<ParticleContainer<3>>();

    // Initialize particles
    for (int i = 0; i < 100; ++i) {
      std::array<double, 3> coordinates = {static_cast<double>(i), 0.0, 0.0};
      std::array<double, 3> velocity = {static_cast<double>(i) / 10.0, 0.0, 0.0};
      particles->insert(Particle<3>(coordinates, velocity, 1.0, 1));
    }
  }
};

/**
 * @brief Test to verify the calculation of minimum and maximum x values.
 */
TEST_F(ProfileCalculatorTest, MinMaxXValues) {
  profileCalculator.updateProfiles(*particles);

  double xMin = profileCalculator.getXMin();
  double xMax = profileCalculator.getXMax();

  EXPECT_NEAR(xMin, 0.0, 1e-6);
  EXPECT_NEAR(xMax, 99.0, 1e-6);
}

/**
 * @brief Test to verify the calculation of bin size.
 */
TEST_F(ProfileCalculatorTest, BinSizeCalculation) {
  profileCalculator.updateProfiles(*particles);

  double binSize = profileCalculator.getBinSize();

  EXPECT_NEAR(binSize, 9.9, 1e-6);
}
