#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-avoid-magic-numbers"

#include "lib/Particle.h"
#include "lib/ParticleContainer.h"
#include <gtest/gtest.h>
#include <array>
#include <iostream>

class PlanetaryCalculationTest : public ::testing::Test {
 public:


  void SetUp() override {

  }
};

TEST_F(PlanetaryCalculationTest, BasicTest) {
  SetUp();

  EXPECT_NE(1, 2);
}



auto main(int argc, char **argv) -> int {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
#pragma clang diagnostic pop
