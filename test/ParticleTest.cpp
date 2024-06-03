#include "Particle.h"
#include "container/container.h"
#include <array>
#include <gtest/gtest.h>

class ParticleTest : public ::testing::Test {
 public:
  container::particle_container container, container2;
  std::vector<Particle> particles, particles2;

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
 * Test: iterator_not_at_end
 *
 * Verifies that incrementing an iterator moves it from its original position and that it does not reach the end of the container.
 *
 * Ensures the iterator is correctly incremented and is not equal to the original position or the end of the container.
 */
TEST_F(ParticleTest, iterator_not_at_end) {
  auto range = container.linear();
  std::visit(
      [](auto &r) {
        auto it = r.begin();
        auto original = it;
        ++it;
        EXPECT_NE(it, original);
        EXPECT_NE(it, r.end());
      },
      range);
}

/**
 * Test: iterator_container_size
 *
 * Verifies that the size of the container is correctly reported and that iterating through the container covers all elements.
 *
 * Ensures the container size is as expected and the iterator correctly counts all elements.
 */
TEST_F(ParticleTest, iterator_container_size) {
  EXPECT_EQ(container.size(), 3);
  auto range = container.linear();
  std::visit(
      [](auto &r) {
        auto it = r.begin();
        int counter = 0;
        while (it != r.end()) {
          ++it;
          counter++;
        }
        EXPECT_EQ(counter, 3);
      },
      range);
}

/**
 * Test: iterator_reaching_end
 *
 * Verifies that an iterator can traverse the entire container and reach the end.
 *
 * Ensures that the iterator reaches the end of the container after iterating through all elements.
 */
TEST_F(ParticleTest, iterator_reaching_end) {
  auto range = container.linear();
  std::visit(
      [](auto &r) {
        auto it = r.begin();

        while (it != r.end()) {
          ++it;
        }

        auto endIt = it;
        EXPECT_TRUE(it == endIt);
      },
      range);
}

/**
 * Test: iterator_pair_building_simple
 *
 * Verifies that the container correctly forms pairs of particles and that these pairs can be iterated over.
 *
 * Manually compares the positions and velocities of particles in pairs to expected values.
 *
 * Ensures the first pair consists of Particle A and Particle B, the second pair consists of Particle A and Particle C,
 * and the third pair consists of Particle B and Particle C.
 */
TEST_F(ParticleTest, iterator_pair_building_simple) {
  auto range = container.pairwise();
  std::visit(
      [](auto &r) {
        auto pair = r.begin();

        auto [pair1_par1, pair1_par2] = *pair;

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
        EXPECT_EQ(pair1_par1.position, coordinatesA);
        EXPECT_EQ(pair1_par1.velocity, velocityA);

        EXPECT_EQ(pair1_par2.position, coordinatesB);
        EXPECT_EQ(pair1_par2.velocity, velocityB);

        // increment and read new pair
        pair++;
        auto [pair2_par1, pair2_par2] = *pair;

        // second pair: 1 and 3
        EXPECT_EQ(pair2_par1.position, coordinatesA);
        EXPECT_EQ(pair2_par1.velocity, velocityA);

        EXPECT_EQ(pair2_par2.position, coordinatesC);
        EXPECT_EQ(pair2_par2.velocity, velocityC);

        // increment and read new pair
        pair++;
        auto [pair3_par1, pair3_par2] = *pair;

        // third pair: 2 and 3
        EXPECT_EQ(pair3_par1.position, coordinatesB);
        EXPECT_EQ(pair3_par1.velocity, velocityB);

        EXPECT_EQ(pair3_par2.position, coordinatesC);
        EXPECT_EQ(pair3_par2.velocity, velocityC);
      },
      range);
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
    std::visit(
        [expected](auto &r) {
          auto pair = r.begin();
          auto count = 0;

          while (pair != r.pair()) {
            pair++;
            count++;
          }
          EXPECT_EQ(expected, count);
        },
        p.pairwise());
  };

  auto amount_of_pairs1 = 3;
  inner_test(container, amount_of_pairs1);

  auto amount_of_pairs2 = (15 * (15 - 1)) / 2;
  inner_test(container2, amount_of_pairs2);
}

/**
 * Test: multi_increment_stability
 *
 * Verifies that an iterator remains stable and correctly reaches the end of the container after multiple increments.
 *
 * Ensures that incrementing the iterator eight times results in it being equal to the end of the container.
 */
TEST_F(ParticleTest, multi_increment_stability) {

  auto range = container.linear();
  std::visit(
      [](auto &r) {
        auto it = r.begin();

        for (int i = 0; i < 8; ++i) {
          if (it != r.end()) {
            ++it;
          }
        }

        EXPECT_TRUE(it == r.end());
      },
      range);
}