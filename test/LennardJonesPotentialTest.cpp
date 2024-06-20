#include "simulator/physics/LennardJones.h"
#include <cmath>
#include <gtest/gtest.h>

class LennardJonesPotentialTest : public ::testing::Test {
 public:
  void SetUp() override {
    const std::vector<double> epsilons = {5.0};
    const std::vector<double> sigmas = {1.0};
    const double cutoff = 3.0;
    simulator::physics::lennard_jones::initialize_constants(epsilons, sigmas, cutoff);
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


  // velocity and mass do not matter in this calculation, therefore set to zero an 1.0

  // Particle i
  const std::array<double_v, 3> x_i = {1.0, 0.0, 0.0};
  const std::array<double_v, 3> velocity_i = {0.0, 0.0, 0.0};

  // Particle j
  std::array<double_v, 3> x_j = {0.0, 1.0, 1.0};
  std::array<double_v, 3> velocity_j = {0.0, 0.0, 0.0};

  VectorizedParticle particle_i = VectorizedParticle(
      x_i,
      velocity_i,
      {double_v(0.0), double_v(0.0), double_v(0.0)},
      {double_v(0.0), double_v(0.0), double_v(0.0)},
      double_v(0.0), long_v(0), double_mask(true));

  VectorizedParticle particle_j = VectorizedParticle(
      x_j,
      velocity_j,
      {0, 0, 0},
      {0, 0, 0},
      double_v(1.0), long_v(0), double_mask(true));

  // r_ij = displacement of vectors = || x_i - x_j ||
  // manual computation of forces
  //    x_i = [1.0, 0.0, 0.0]
  //    x_j = [0.0, 1.0, 1.0]
  // => r_ij = [1.0, -1.0, -1.0]

  // computation of Lennard-Jones force according to the formulas on the worksheet

  // Lennard-Jones potential U(x_i, x_j) = -0.7133058984910836
  // Lennard-Jones force F_ij = [-1.37174211, 1.37174211, 1.37174211]

  const std::array<double_v, 3> actual_forces = {-1.37174211, 1.37174211, 1.37174211};

  // computing forces with method
  std::array<double_v, 3> calculated_forces{};
  std::array<double_v, 3> correction = {double_v(0), double_v(0), double_v(0)};
  simulator::physics::lennard_jones::calculate_force_vectorized(
      particle_i,
      particle_j,
      double_mask(true),
      calculated_forces,
      correction);

  // check if each result of the force vector is accurate enough
  EXPECT_TRUE(stdx::all_of(calculated_forces[0] - actual_forces[0] < 0.000001));
  EXPECT_TRUE(stdx::all_of(calculated_forces[1] - actual_forces[1] < 0.000001));
  EXPECT_TRUE(stdx::all_of(calculated_forces[2] - actual_forces[2] < 0.000001));
}

/**
* Test: lennard_jones_small_forces
*
* Verifies the correctness of the Lennard-Jones force calculation between two particles
* with small force magnitudes -> cutoff used.
*/
TEST_F(LennardJonesPotentialTest, lennard_jones_small_forces) {

  // to change constants later, when type of atom/molecule is adjustable
  // double epsilon = 4.30483;
  // double sigma = 2.121;

  // velocity and mass do not matter in this calculation, therefore set to zero an 1.0


  // Particle i
  const std::array<double_v, 3> x_i = {2.1, 3.5, 1.4};
  const std::array<double_v, 3> velocity_i = {0.0, 0.0, 0.0};

  // Particle j
  std::array<double_v, 3> x_j = {0.004, 0.0584, 0.0372};
  std::array<double_v, 3> velocity_j = {0.0, 0.0, 0.0};

  VectorizedParticle particle_i = VectorizedParticle(
      x_i,
      velocity_i,
      {double_v(0.0), double_v(0.0), double_v(0.0)},
      {double_v(0.0), double_v(0.0), double_v(0.0)},
      double_v(0.0), long_v(0), double_mask(true));

  VectorizedParticle particle_j = VectorizedParticle(
      x_j,
      velocity_j,
      {0, 0, 0},
      {0, 0, 0},
      double_v(1.0), long_v(0), double_mask(true));

  // r_ij = displacement of vectors = || x_i - x_j ||
  // manual computation of forces
  //    x_i = [2.1, 3.5, 1.4]
  //    x_j = [0.004, 0.0584, 0.0372]
  // => r_ij = [2.096  3.4416 1.3628]

  // computation of Lennard-Jones force according to the formulas on the worksheet

  // Lennard-Jones potential U(x_i, x_j) = -0.0033750273618335614
  // Lennard-Jones force F_ij = [-0.00234524 -0.00385084 -0.00152485]

  std::array<double_v, 3> actual_forces = {0.0, 0.0, 0.0}; // using cutoff so the forces are too small

  // computing forces with method
  std::array<double_v, 3> calculated_forces{};
  std::array<double_v, 3> correction = {double_v(0), double_v(0), double_v(0)};
  simulator::physics::lennard_jones::calculate_force_vectorized(
      particle_i,
      particle_j,
      double_mask(true),
      calculated_forces,
      correction);

  // check if each result of the force vector is accurate enough
  EXPECT_TRUE(stdx::all_of(calculated_forces[0] - actual_forces[0] < 0.000001));
  EXPECT_TRUE(stdx::all_of(calculated_forces[1] - actual_forces[1] < 0.000001));
  EXPECT_TRUE(stdx::all_of(calculated_forces[2] - actual_forces[2] < 0.000001));
}


/**
 * Test: lennard_jones_same_particle_not_zero
 *
 * Verifies the behavior of the Lennard-Jones force calculation when the same particle is compared.
 *
 * Ensures that the calculated forces result in NaN values, as the displacement should be zero.
 */
TEST_F(LennardJonesPotentialTest, lennard_jones_same_particle_not_zero) {

  // velocity and mass do not matter in this calculation, therefore set to zero an 1.0

  // Particle i
  std::array<double_v, 3> x_i = {1.0, 1.0, 0.0};
  std::array<double_v, 3> velocity_i = {0.0, 0.0, 0.0};

  // Particle j
  std::array<double_v, 3> x_j = {1.0, 1.0, 0.0};
  std::array<double_v, 3> velocity_j = {0.0, 0.0, 0.0};

  VectorizedParticle particle_i = VectorizedParticle(
      x_i,
      velocity_i,
      {double_v(0.0), double_v(0.0), double_v(0.0)},
      {double_v(0.0), double_v(0.0), double_v(0.0)},
      double_v(0.0), long_v(0), double_mask(true));

  VectorizedParticle particle_j = VectorizedParticle(
      x_j,
      velocity_j,
      {0, 0, 0},
      {0, 0, 0},
      double_v(1.0), long_v(0), double_mask(true));

  std::array<double_v, 3> actual_forces = {NAN, NAN, NAN}; // using cutoff so the forces are too small

  // computing forces with method
  std::array<double_v, 3> calculated_forces{};
  std::array<double_v, 3> correction = {double_v(0), double_v(0), double_v(0)};
  simulator::physics::lennard_jones::calculate_force_vectorized(
      particle_i,
      particle_j,
      double_mask(true),
      calculated_forces,
      correction);

  // check if each result of the force vector is accurate enough

  // TODO: add according EXPECT or ASSERT

}
