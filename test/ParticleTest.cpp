#include <gtest/gtest.h>
#include "../src/ParticleContainer.h"

class ParticleTest : public ::testing::Test {
 protected:
  ParticleContainer container;

  // custom ParticleContainer constructor with capacity 10 because there is no default otherwise
  ParticleTest() : container(10) {}

  void SetUp() override {
    // Initialize container with enough particles to test pair iteration.
    std::vector<Particle> particles;

    // Particle A
    std::array<double, 3> coordinatesA {0.0, 0.0, 0.0};
    std::array<double, 3> velocityA {2.0, 3.0, 3.0};

    // Particle B
    std::array<double, 3> coordinatesB {1.0, 1.0, 2.0};
    std::array<double, 3> velocityB {3.0, 2.0, 3.0};

    // Particle C
    std::array<double, 3> coordinatesC {2.0, 3.0, 4.0};
    std::array<double, 3> velocityC {5.0, 6.0, 3.0};

    particles.emplace_back(coordinatesA, velocityA, 1.0);
    particles.emplace_back(coordinatesA, velocityA, 1.0);
    particles.emplace_back(coordinatesA, velocityA, 1.0);

    container = ParticleContainer(std::move(particles));
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

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}