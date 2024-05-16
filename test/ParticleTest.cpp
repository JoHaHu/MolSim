#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-avoid-magic-numbers"

#include "lib/Particle.h"
#include "lib/ParticleContainer.h"
#include <array>
#include <gtest/gtest.h>
#include <iostream>

class ParticleTest : public ::testing::Test {
 public:
  ParticleContainer container, container2;
  std::vector<Particle> particles, particles2;

  // custom ParticleContainer constructor with capacity 10 because there is no default otherwise
  ParticleTest() : container(20), container2(20) {}

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
    Particle particleA = Particle(coordinatesA, velocityA, 1.0, 0);
    Particle particleB = Particle(coordinatesB, velocityB, 4.0, 0);
    Particle particleC = Particle(coordinatesC, velocityC, 10.0, 0);

    // add test particles to the particles vector
    particles.emplace_back(particleA);
    particles.emplace_back(particleB);
    particles.emplace_back(particleC);

    // add 15 particles to particles2 with different parameter values by using a loop
    for (int i = 0; i < 15; i++) {
      // Particle D
      std::array<double, 3> coordinatesD = {1.0 + i, 2.0 + i, 1.5 + i};
      std::array<double, 3> velocityD = {2.0, 3.0, 3.0};

      Particle particleD = Particle(coordinatesD, velocityD, 15.0 + i, 0);
      particles2.emplace_back(particleD);
    }

    // initialise containers
    container = ParticleContainer(particles);
    container2 = ParticleContainer(particles2);
  }
};

TEST_F(ParticleTest, IncrementOperator_NotAtEnd) {
  SetUp();
  auto it = container.begin();
  auto original = it;

  ++it;

  EXPECT_NE(it, original);
  EXPECT_NE(it, container.end());
}

TEST_F(ParticleTest, Container_Size) {
  SetUp();
  auto it = container.begin();
  int counter = 0;

  EXPECT_EQ(container.size(), 3);

  while (it != container.end()) {
    ++it;
    counter++;
  }
  EXPECT_EQ(counter, 3);
}

TEST_F(ParticleTest, IncrementOperator_ReachesEnd) {
  auto it = container.begin();

  while (it != container.end()) {
    ++it;
  }

  auto endIt = it;
  EXPECT_TRUE(it == endIt);
}

TEST_F(ParticleTest, IncrementOperator_BasicFunction) {
  auto pair = container.begin_pair();

  auto [pair1par1, pair1par2] = *pair;

  // Manually comparing the particles:

  // Particle A
  std::array<double, 3> coordinatesA = {0.0, 0.0, 0.0};
  std::array<double, 3> velocityA = {2.0, 3.0, 3.0};

  // Particle B
  std::array<double, 3> coordinatesB = {1.0, 1.0, 2.0};
  std::array<double, 3> velocityB = {3.0, 2.0, 3.0};

  // Particle C
  std::array<double, 3> coordinatesC = {2.0, 3.0, 4.0};
  std::array<double, 3> velocityC = {5.0, 6.0, 3.0};

  // first pair: 1 and 2
  EXPECT_EQ(pair1par1.position, coordinatesA);
  EXPECT_EQ(pair1par1.velocity, velocityA);

  EXPECT_EQ(pair1par2.position, coordinatesB);
  EXPECT_EQ(pair1par2.velocity, velocityB);

  // increment and read new pair
  pair++;
  auto [pair2par1, pair2par2] = *pair;

  // second pair: 1 and 3
  EXPECT_EQ(pair2par1.position, coordinatesA);
  EXPECT_EQ(pair2par1.velocity, velocityA);

  EXPECT_EQ(pair2par2.position, coordinatesC);
  EXPECT_EQ(pair2par2.velocity, velocityC);

  // increment and read new pair
  pair++;
  auto [pair3par1, pair3par2] = *pair;

  // third pair: 2 and 3
  EXPECT_EQ(pair3par1.position, coordinatesB);
  EXPECT_EQ(pair3par1.velocity, velocityB);

  EXPECT_EQ(pair3par2.position, coordinatesC);
  EXPECT_EQ(pair3par2.velocity, velocityC);

  /** TESTING WITH IDs (work in progress)
   * // TODO: initialise ID properly and compare IDs

  // check for pair 1 (particle 1 and 2) and increment after
  EXPECT_EQ(particle1.id, 1);
  EXPECT_EQ(particle2.id, 2);
  pair++;

  // check for pair 1 (particle 2 and 3) and increment after
  EXPECT_EQ(particle1.id, 2);
  EXPECT_EQ(particle2.id, 3);
  pair++;

  // check for pair 1 (particle 1 and 3) and increment after
  EXPECT_EQ(particle1.id, 1);
  EXPECT_EQ(particle2.id, 3);

  // check if pair reached end

  **/
}

TEST_F(ParticleTest, IncrementOperator_LargerVector) {
  auto pair1 = container.begin_pair();
  auto pair2 = container2.begin_pair();

  // number of pairs can be calculated with the handshake lemma
  auto amountOfPairs1 = 3;
  auto amountOfPairs2 = (15 * (15 - 1)) / 2;

  auto count1 = 0;
  auto count2 = 0;

  while (pair1 != container.end_pair()) {
    pair1++;
    count1++;
  }

  while (pair2 != container2.end_pair()) {
    pair2++;
    count2++;
  }

  // check if amount of pairs matches (according to number of particles)
  EXPECT_EQ(amountOfPairs1, count1);
  EXPECT_EQ(amountOfPairs2, count2);
}

TEST_F(ParticleTest, IncrementOperator_MultiIncrementStability) {
  auto it = container.begin();

  for (int i = 0; i < 8; ++i) {
    if (it != container.end()) {
      ++it;
    }
  }

  EXPECT_TRUE(it == container.end());
}

auto main(int argc, char **argv) -> int {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
#pragma clang diagnostic pop