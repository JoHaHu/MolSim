#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-avoid-magic-numbers"

#include "lib/Particle.h"
#include "lib/ParticleContainer.h"
#include <gtest/gtest.h>

class ParticleTest : public ::testing::Test {
 private:
  ParticleContainer container;

 public:
  ParticleTest() : container(10) {}

  void SetUp() override {
    // Initialize container with enough particles to test pair iteration.
    std::vector<Particle> particles;

    // Particle A
    std::array<double, 3> coordinatesA{0.0, 0.0, 0.0};
    std::array<double, 3> velocityA{2.0, 3.0, 3.0};

    // Particle B
    std::array<double, 3> coordinatesB{1.0, 1.0, 2.0};
    std::array<double, 3> velocityB{3.0, 2.0, 3.0};

    // Particle C
    std::array<double, 3> coordinatesC{2.0, 3.0, 4.0};
    std::array<double, 3> velocityC{5.0, 6.0, 3.0};

    particles.emplace_back(coordinatesA, velocityA, 1.0);
    particles.emplace_back(coordinatesB, velocityB, 1.0);
    particles.emplace_back(coordinatesC, velocityC, 1.0);

    container = ParticleContainer(particles);
  }
};

TEST_F(ParticleTest, IncrementOperator_NotAtEnd) {
  SetUp();
  //auto it = container.begin();
  //auto original = it;

  //++it;

  //EXPECT_NE(it, original);
  //EXPECT_NE(it, container.end());
}

TEST(SimpleTest, BasicAssertions) {
  EXPECT_EQ(1, 1);
}

TEST(PairIteratorTest, BasicAssertions) {
  // Expect two strings not to be equal.
  EXPECT_STRNE("hello", "world");
  // Expect equality.
  EXPECT_EQ(7 * 6, 42);
}

auto main(int argc, char **argv) -> int {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
#pragma clang diagnostic pop