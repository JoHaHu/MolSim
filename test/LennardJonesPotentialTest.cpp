#include "simulator/physics/LennardJones.h"
#include <cmath>
#include <gtest/gtest.h>

class LennardJonesPotentialTest : public ::testing::Test {
 public:
  std::unique_ptr<simulator::physics::LennardJonesForce> lennardJonesForce;

  void SetUp() override {
    const std::vector<double> epsilons = {5.0};
    const std::vector<double> sigmas = {1.0};
    const double cutoff = 3.0;
    lennardJonesForce = std::make_unique<simulator::physics::LennardJonesForce>(cutoff, epsilons, sigmas);
  }
};

/**
 * Test: lennard_jones_forces_simple
 *
 * Verifies the correctness of the Lennard-Jones force calculation between two particles.
 *
 * Verifies that each component of the calculated force vector is accurate within 0.000001.
 */
TEST_F(LennardJonesPotentialTest, lennard_jones_forces_simple) {
  // Particle i
  const std::array<double, 3> x_i = {1.0, 0.0, 0.0};
  const std::array<double, 3> velocity_i = {0.0, 0.0, 0.0};

  // Particle j
  std::array<double, 3> x_j = {0.0, 1.0, 1.0};
  std::array<double, 3> velocity_j = {0.0, 0.0, 0.0};

  Particle<3> particle_i(x_i, velocity_i, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 0, 0);
  Particle<3> particle_j(x_j, velocity_j, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 0, 0);

  // Manual computation of forces
  const std::array<double, 3> actual_forces = {-1.37174211, 1.37174211, 1.37174211};

  // Computing forces with method
  std::array<double, 3> calculated_forces{};
  std::array<double, 3> correction = {0.0, 0.0, 0.0};
  long type_i = particle_i.type;
  long type_j = particle_j.type;
  lennardJonesForce->calculateForce_3D(
      particle_i.position[0], particle_i.position[1], particle_i.position[2], particle_i.mass, type_i,
      particle_j.position[0], particle_j.position[1], particle_j.position[2], particle_j.mass, type_j,
      calculated_forces[0], calculated_forces[1], calculated_forces[2], correction);

  // Check if each result of the force vector is accurate enough
  EXPECT_NEAR(calculated_forces[0], actual_forces[0], 0.000001);
  EXPECT_NEAR(calculated_forces[1], actual_forces[1], 0.000001);
  EXPECT_NEAR(calculated_forces[2], actual_forces[2], 0.000001);
}

/**
* Test: lennard_jones_small_forces
*
* Verifies the correctness of the Lennard-Jones force calculation between two particles
* with small force magnitudes -> cutoff used.
*/
TEST_F(LennardJonesPotentialTest, lennard_jones_small_forces) {
  // Particle i
  const std::array<double, 3> x_i = {2.1, 3.5, 1.4};
  const std::array<double, 3> velocity_i = {0.0, 0.0, 0.0};

  // Particle j
  std::array<double, 3> x_j = {0.004, 0.0584, 0.0372};
  std::array<double, 3> velocity_j = {0.0, 0.0, 0.0};

  Particle<3> particle_i(x_i, velocity_i, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 0, 0);
  Particle<3> particle_j(x_j, velocity_j, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 0, 0);

  // Using cutoff so the forces are too small
  std::array<double, 3> actual_forces = {0.0, 0.0, 0.0};

  // Computing forces with method
  std::array<double, 3> calculated_forces{};
  std::array<double, 3> correction = {0.0, 0.0, 0.0};
  long type_i = particle_i.type;
  long type_j = particle_j.type;
  lennardJonesForce->calculateForce_3D(
      particle_i.position[0], particle_i.position[1], particle_i.position[2], particle_i.mass, type_i,
      particle_j.position[0], particle_j.position[1], particle_j.position[2], particle_j.mass, type_j,
      calculated_forces[0], calculated_forces[1], calculated_forces[2], correction);

  // Check if each result of the force vector is accurate enough
  EXPECT_NEAR(calculated_forces[0], actual_forces[0], 0.000001);
  EXPECT_NEAR(calculated_forces[1], actual_forces[1], 0.000001);
  EXPECT_NEAR(calculated_forces[2], actual_forces[2], 0.000001);
}

/**
 * Test: lennard_jones_same_particle_not_zero
 *
 * Verifies the behavior of the Lennard-Jones force calculation when the same particle is compared.
 *
 * Ensures that the calculated forces result in NaN values, as the displacement should be zero.
 */
TEST_F(LennardJonesPotentialTest, lennard_jones_same_particle) {
  // Particle i
  std::array<double, 3> x_i = {1.0, 1.0, 0.0};
  std::array<double, 3> velocity_i = {0.0, 0.0, 0.0};

  // Particle j
  std::array<double, 3> x_j = {1.0, 1.0, 0.0};
  std::array<double, 3> velocity_j = {0.0, 0.0, 0.0};

  Particle<3> particle_i(x_i, velocity_i, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 0, 0);
  Particle<3> particle_j(x_j, velocity_j, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0, 0, 0);

  std::array<double, 3> actual_forces = {NAN, NAN, NAN};

  // Computing forces with method
  std::array<double, 3> calculated_forces{};
  std::array<double, 3> correction = {0.0, 0.0, 0.0};
  long type_i = particle_i.type;
  long type_j = particle_j.type;
  lennardJonesForce->calculateForce_3D(
      particle_i.position[0], particle_i.position[1], particle_i.position[2], particle_i.mass, type_i,
      particle_j.position[0], particle_j.position[1], particle_j.position[2], particle_j.mass, type_j,
      calculated_forces[0], calculated_forces[1], calculated_forces[2], correction);

  // Check if each result of the force vector is accurate enough
  EXPECT_TRUE(std::isnan(calculated_forces[0]));
  EXPECT_TRUE(std::isnan(calculated_forces[1]));
  EXPECT_TRUE(std::isnan(calculated_forces[2]));
}
