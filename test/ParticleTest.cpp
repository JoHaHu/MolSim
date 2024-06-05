#include "Particle.h"
#include "container/container.h"
#include <array>
#include <gtest/gtest.h>

class ParticleTest : public ::testing::Test {
 public:
  std::vector<Particle> particles{}, particles2{};
  container::particle_container container, container2;

  // custom ParticleContainer constructor with capacity 10 because there is no default otherwise
  ParticleTest() : container(particles), container2(particles2) {}

  void SetUp() override {

    // Particle A
    std::array<double, 3> coordinates_A = {0.0, 0.0, 0.0};
    std::array<double, 3> velocity_A = {2.0, 3.0, 3.0};

    // Particle B
    std::array<double, 3> coordinates_B = {1.0, 1.0, 2.0};
    std::array<double, 3> velocity_B = {3.0, 2.0, 3.0};

    // Particle C
    std::array<double, 3> coordinates_C = {2.0, 3.0, 4.0};
    std::array<double, 3> velocity_C = {5.0, 6.0, 3.0};

    // construct particles with different parameter values
    Particle particleA = Particle(coordinates_A, velocity_A, 1.0, 0);
    Particle particleB = Particle(coordinates_B, velocity_B, 4.0, 0);
    Particle particleC = Particle(coordinates_C, velocity_C, 10.0, 0);

    // add test particles to the particles vector
    particles.push_back(particleA);
    particles.push_back(particleB);
    particles.push_back(particleC);

    // add 15 particles to particles2 with different parameter values by using a loop
    for (int i = 0; i < 15; i++) {
      // Particle D
      std::array<double, 3> coordinatesD = {1.0 + i, 2.0 + i, 1.5 + i};
      std::array<double, 3> velocityD = {2.0, 3.0, 3.0};

      Particle particleD = Particle(coordinatesD, velocityD, 15.0 + i, 0);
      particles2.push_back(particleD);
    }

    // initialise containers
    container = container::particle_container(particles);
    container2 = container::particle_container(particles2);
  }
};

/**
 * Test: iterator_container_size
 *
 * Verifies that the size of the container is correctly reported and that iterating through the container covers all elements.
 *
 * Ensures the container size is as expected and the iterator correctly counts all elements.
 */
TEST_F(ParticleTest, iterator_container_size) {
  EXPECT_EQ(container.size(), 3);
  auto count = 0;
  container.linear([&count](auto &p) {
    count++;
  });
  EXPECT_EQ(count, 3);
}

/**
 * Test: iterator_pair_building_large
 *
 * Verifies that the container correctly forms all possible pairs of particles and that these pairs can be iterated over.
 *
 * Checks the number of pairs using the handshake lemma.
 *
 * Ensures that the number of pairs in container1 and container2 matches the expected counts based on the number of particles.
 */
TEST_F(ParticleTest, iterator_pair_building_large) {

  auto inner_test = [](container::particle_container &p, long expected) {
    auto count = 0;
    p.pairwise([&count](auto t) { count++; });
    EXPECT_EQ(expected, count);
  };

  auto amount_of_pairs1 = 3;
  inner_test(container, amount_of_pairs1);

  auto amount_of_pairs2 = (15 * (15 - 1)) / 2;
  inner_test(container2, amount_of_pairs2);
}
