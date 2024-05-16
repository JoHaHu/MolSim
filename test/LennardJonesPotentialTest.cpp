//
// Created by TimSc on 16.05.2024.
//

#include <gtest/gtest.h>

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-avoid-magic-numbers"



class LennardJonesPotentialTest : public ::testing::Test {
 public:


  // no additional setup necessary here

};

TEST_F(LennardJonesPotentialTest, lennard_jones_forces_simple) {

  // initialise input parameters with simple values
  double epsilon = 1.0;
  double sigma = 1.0;

  std::array<double, 3> x_i = {1.0, 0.0, 0.0};
  std::array<double, 3> x_j = {1.0, 0.0, 0.0};

  // r_ij = displacement of vectors = || x_i - x_j ||
  // manual computation of forces
  //    x_i = [1.0, 0.0, 0.0]
  //    x_j = [1.0, 0.0, 0.0]
  // => r_ij = [1.0, 1.0, 0.1]

  // computation of Lennard-Jones force according to the formulas on the worksheet

  // Lennard-Jones potential U(x_i, x_j) = -0.4375
  // Lennard-Jones force F_ij = [-1.125, 1.125, -0.0]

  std::array<double,3> actual_forces = {-1.125, 1.125, -0.0};

  // computing forces with method
  // TODO: insert method name and replace this hard-coded array
  std::array<double,3> calculated_forces = {-1.125, 1.125, -0.0};

  EXPECT_EQ(actual_forces, calculated_forces);
}

TEST_F(LennardJonesPotentialTest, lennard_jones_small_forces) {

  // set input parameters with arbitrary values (no real life relation yet) and calculating forces
  double epsilon = 4.30483;
  double sigma = 0.397493;

  std::array<double, 3> x_i = {2.1, 3.5, 1.4};
  std::array<double, 3> x_j = {0.004, 0.0584, 0.0372};

  // r_ij = displacement of vectors = || x_i - x_j ||
  // manual computation of forces
  //    x_i = [2.1, 3.5, 1.4]
  //    x_j = [0.004, 0.0584, 0.0372]
  // => r_ij = [2.096  3.4416 1.3628]

  // computation of Lennard-Jones force according to the formulas on the worksheet

  // Lennard-Jones potential U(x_i, x_j) = -1.1463394182511578e-05
  // Lennard-Jones force F_ij = [-7.96701562e-06, -1.30817180e-05, -5.18008058e-06]

  std::array<double,3> actual_forces = {-7.96701562e-06, -1.30817180e-05, -5.18008058e-06};

  // computing forces with method
  // TODO: insert method name and replace this hard-coded array
  std::array<double,3> calculated_forces = {-7.96701562e-06, -1.30817180e-05, -5.18008058e-06};

  EXPECT_EQ(actual_forces, calculated_forces);
}

TEST_F(LennardJonesPotentialTest, lennard_jones_large_forces) {

  // set input parameters with arbitrary values (no real life relation yet) and calculating forces
  double epsilon = 0.2453234;
  double sigma = 4.482;

  std::array<double, 3> x_i = {0.001, 0.0532, 0.0361};
  std::array<double, 3> x_j = {0.004, 0.0584, 0.0372};

  // r_ij = displacement of vectors = || x_i - x_j ||
  // manual computation of forces
  //    x_i = [0.001, 0.0532, 0.0361]
  //    x_j = [0.004, 0.0584, 0.0372]
  // => r_ij = [-0.003, -0.0052, -0.0011]

  // computation of Lennard-Jones force according to the formulas on the worksheet

  // Lennard-Jones potential U(x_i, x_j) = 2.4138124905450815e+34
  // Lennard-Jones force F_ij = [-2.33281207e+37,-4.04354092e+37, -8.55364426e+36]

  std::array<double,3> actual_forces = {-2.33281207e+37, -4.04354092e+37, -8.55364426e+36};

  // computing forces with method
  // TODO: insert method name and replace this hard-coded array
  std::array<double,3> calculated_forces = {-2.33281207e+37, -4.04354092e+37, -8.55364426e+36};

  EXPECT_EQ(actual_forces, calculated_forces);
}


TEST_F(LennardJonesPotentialTest, lennard_jones_forces_argon) {

  // set input parameters with the values for argon
  double epsilon = 0.238;
  double sigma = 3.4;

  std::array<double, 3> x_i = {1.1, 2.0, 0.9};
  std::array<double, 3> x_j = {2.0, 1.1, 2.2};

  // r_ij = displacement of vectors = || x_i - x_j ||
  // manual computation of forces
  //    x_i = [1.1, 2.0, 0.9]
  //    x_j = [2.0, 1.1, 2.2]
  // => r_ij = [-0.9, 0.9, -1.3]

  // computation of Lennard-Jones force according to the formulas on the worksheet

  // Lennard-Jones potential U(x_i, x_j) = 1686.938457711386
  // Lennard-Jones force F_ij = [-5570.3695767, 5570.3695767, -8046.08938856]

  std::array<double,3> actual_forces = {-5570.3695767, 5570.3695767, -8046.08938856};

  // computing forces with method
  // TODO: insert method name and replace this hard-coded array
  std::array<double,3> calculated_forces = {-5570.3695767, 5570.3695767, -8046.08938856};

  EXPECT_EQ(actual_forces, calculated_forces);
}

TEST_F(LennardJonesPotentialTest, lennard_jones_forces_xenon) {

  // set input parameters with the values for xenon
  double epsilon = 0.399;
  double sigma = 4.1;

  std::array<double, 3> x_i = {1.1, 2.0, 0.9};
  std::array<double, 3> x_j = {2.0, 1.1, 2.2};

  // r_ij = displacement of vectors = || x_i - x_j ||
  // manual computation of forces
  //    x_i = [1.1, 2.0, 0.9]
  //    x_j = [2.0, 1.1, 2.2]
  // => r_ij = [-0.9, 0.9, -1.3]

  // computation of Lennard-Jones force according to the formulas on the worksheet

  // Lennard-Jones potential U(x_i, x_j) = 27173.32984787372
  // Lennard-Jones force F_ij = [-89003.27408715, 89003.27408715, -128560.28479256]

  std::array<double,3> actual_forces = {-89003.27408715, 89003.27408715, -128560.28479256};

  // computing forces with method
  // TODO: insert method name and replace this hard-coded array
  std::array<double,3> calculated_forces = {-89003.27408715, 89003.27408715, -128560.28479256};

  EXPECT_EQ(actual_forces, calculated_forces);
}


auto main(int argc, char **argv) -> int {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
#pragma clang diagnostic pop