#include "Particle.h"
#include "container/ParticleContainer.h"
#include <array>
#include <gtest/gtest.h>

class ParticleTest : public ::testing::Test {
 public:
  std::vector<Particle<3>> particles{}, particles2{};
  std::unique_ptr<ParticleContainer<3>> container, container2;

  ParticleTest() = default;

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
    Particle<3> particleA(coordinates_A, velocity_A, 1.0, 0);
    Particle<3> particleB(coordinates_B, velocity_B, 4.0, 0);
    Particle<3> particleC(coordinates_C, velocity_C, 10.0, 0);

    // add test particles to the particles vector
    particles.push_back(particleA);
    particles.push_back(particleB);
    particles.push_back(particleC);

    // add 15 particles to particles2 with different parameter values by using a loop
    for (int i = 0; i < 15; ++i) {
      // Particle D
      std::array<double, 3> coordinatesD = {1.0 + i, 2.0 + i, 1.5 + i};
      std::array<double, 3> velocityD = {2.0, 3.0, 3.0};

      Particle<3> particleD(coordinatesD, velocityD, 15.0 + i, 0);
      particles2.push_back(particleD);
    }

    // initialise containers
    container = std::make_unique<ParticleContainer<3>>();
    container2 = std::make_unique<ParticleContainer<3>>();

    for (const auto &particle : particles) {
      container->insert(particle);
    }

    for (const auto &particle : particles2) {
      container2->insert(particle);
    }
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
  EXPECT_EQ(container->size(), 3);
  auto count = 0;
  container->linear([&count](Particles<3> &, size_t) {
    count++;
  });
  EXPECT_EQ(count, 3);
}
