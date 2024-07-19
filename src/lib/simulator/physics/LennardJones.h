#pragma once

#include "Force.h"
#include "Particle.h"
#include "experimental/simd"
#include "immintrin.h"
#include "spdlog/spdlog.h"
#include "utils/ArrayUtils.h"

namespace simulator::physics {

/**
 * @brief Implements Lennard-Jones force model for simulations.
 */
class LennardJonesForce final : public Force {

  double cutoff = 3.0;

  static const long NUM_TYPES = 8;

  std::array<double, NUM_TYPES> membrane_cutoff;
  // For now supports up to 8 types, can be increased

  alignas(64) std::array<double, NUM_TYPES * NUM_TYPES> epsilons_24 = std::array<double, NUM_TYPES * NUM_TYPES>();
  alignas(64) std::array<double, NUM_TYPES * NUM_TYPES> epsilons_48 = std::array<double, NUM_TYPES * NUM_TYPES>();

 public:
  /**
   * @brief Constructor for LennardJonesForce.
   * @param cutoff Cutoff distance for the force calculation.
   * @param epsilons Vector of epsilon values for different particle types.
   * @param sigmas Vector of sigma values for different particle types.
   */
  explicit LennardJonesForce(double cutoff, std::vector<double> epsilons, std::vector<double> sigmas) {
    initialize_constants(epsilons, sigmas, cutoff);
  }

  /**
   * @brief Calculates 2D Lennard-Jones force between particles.
   * @param x1 X-coordinate of the first particle.
   * @param y1 Y-coordinate of the first particle.
   * @param mass1 Mass of the first particle.
   * @param type1 Type of the first particle.
   * @param x2 X-coordinate of the second particle.
   * @param y2 Y-coordinate of the second particle.
   * @param mass2 Mass of the second particle.
   * @param type2 Type of the second particle.
   * @param result_x Resultant force in the x-direction.
   * @param result_y Resultant force in the y-direction.
   * @param correction Correction factors.
   */
#pragma omp declare simd simdlen(4) uniform(this, x1, y1, mass1, type1, correction, membrane1) linear(ref(x2, y2, mass2, type2, membrane2))
#pragma omp declare simd simdlen(8) uniform(this, x1, y1, mass1, type1, correction, membrane1) linear(ref(x2, y2, mass2, type2, membrane2))
  inline void calculateForce_2D(
      double const &x1,
      double const &y1,
      double const &mass1,
      long const &type1,
      uint8_t const &membrane1,
      double &x2,
      double &y2,
      double &mass2,
      long &type2,
      uint8_t &membrane2,
      double &result_x,
      double &result_y,
      std::array<double, 2> &correction) override {

    const auto diff = (std::array<double, 2>({x2, y2}) + correction) - std::array<double, 2>({x1, y1});
    const auto norm_2 = ArrayUtils::L2NormSquared(diff);
    const auto norm_6 = norm_2 * norm_2 * norm_2;
    auto epsilon24 = epsilons_24[8 * type1 + type2];
    auto epsilon48 = epsilons_48[8 * type1 + type2];
    const auto temp = (norm_6 * epsilon24 - epsilon48) / (norm_6 * norm_6 * norm_2);
    auto force = temp * diff;
    if (norm_2 <= cutoff && ((membrane1 > 0 && membrane2 > 0 && norm_2 <= membrane_cutoff[type1]) || membrane1 <= 0 || membrane2 <= 0)) {
      result_x = force[0];
      result_y = force[1];
    }
  }

  /**
   * @brief Calculates 3D Lennard-Jones force between particles.
   * @param x1 X-coordinate of the first particle.
   * @param y1 Y-coordinate of the first particle.
   * @param z1 Z-coordinate of the first particle.
   * @param mass1 Mass of the first particle.
   * @param type1 Type of the first particle.
   * @param x2 X-coordinate of the second particle.
   * @param y2 Y-coordinate of the second particle.
   * @param z2 Z-coordinate of the second particle.
   * @param mass2 Mass of the second particle.
   * @param type2 Type of the second particle.
   * @param result_x Resultant force in the x-direction.
   * @param result_y Resultant force in the y-direction.
   * @param result_z Resultant force in the z-direction.
   * @param correction Correction factors.
   */
#pragma omp declare simd simdlen(4) uniform(this, x1, y1, mass1, type1, correction, membrane1) linear(ref(x2, y2, mass2, type2, membrane2))
#pragma omp declare simd simdlen(8) uniform(this, x1, y1, mass1, type1, correction, membrane1) linear(ref(x2, y2, mass2, type2, membrane2))
  inline void calculateForce_3D(
      double const &x1,
      double const &y1,
      double const &z1,
      double const &mass1,
      long const &type1,
      uint8_t const &membrane1,
      double &x2,
      double &y2,
      double &z2,
      double &mass2,
      long &type2,
      uint8_t &membrane2,
      double &result_x,
      double &result_y,
      double &result_z,
      std::array<double, 3> &correction) override {

    const auto diff = (std::array<double, 3>({x2, y2, z2}) + correction) - std::array<double, 3>({x1, y1, z1});
    const auto norm_2 = ArrayUtils::L2NormSquared(diff);
    const auto norm_6 = norm_2 * norm_2 * norm_2;
    auto epsilon24 = epsilons_24[8 * type1 + type2];
    auto epsilon48 = epsilons_48[8 * type1 + type2];
    const auto temp = (norm_6 * epsilon24 - epsilon48) / (norm_6 * norm_6 * norm_2);
    auto force = temp * diff;

    if (norm_2 <= cutoff && ((membrane1 > 0 && membrane2 > 0 && norm_2 <= membrane_cutoff[type1]) || membrane1 <= 0 || membrane2 <= 0)) {
      result_x = force[0];
      result_y = force[1];
      result_z = force[2];
    }
  }

  /**
   * @brief Calculates Lennard-Jones force for a boundary particle.
   * @param diff Difference in position.
   * @param type Type of the particle.
   * @return Calculated boundary force.
   */
  double inline calculate_boundary_force(double diff, int type) override {
    const auto norm_2 = diff * diff;
    const auto norm_6 = norm_2 * norm_2 * norm_2;
    const auto temp = (norm_6 * epsilons_24[8 * type + type] - epsilons_48[8 * type + type]) / (norm_6 * norm_6 * norm_2);
    return temp * diff;
  }

 private:
  /**
   * @brief Initializes the Lennard-Jones constants.
   * @param epsilons Vector of epsilon values for different particle types.
   * @param sigmas Vector of sigma values for different particle types.
   * @param c Cutoff distance for the force calculation.
   */
  auto initialize_constants(std::vector<double> epsilons, std::vector<double> sigmas, double c) -> void {

    cutoff = c * c;
    for (int i = 0; i < sigmas.size(); ++i) {
      membrane_cutoff[i] = pow(2, 1 / 6) * sigmas[i] * pow(2, 1 / 6) * sigmas[i];
    }

    for (size_t i = 0; i < epsilons.size(); ++i) {
      for (size_t j = 0; j < epsilons.size(); ++j) {
        auto sigma = (sigmas[i] + sigmas[j]) / 2;
        auto sigma3 = sigma * sigma * sigma;
        auto sigma6 = sigma3 * sigma3;
        auto epsilon = std::sqrt(epsilons[i] * epsilons[j]);
        auto epsilon24 = 24 * epsilon * sigma6;
        auto epsilon48 = 2 * epsilon24 * sigma6;
        epsilons_24[i + 8 * j] = epsilon24;
        epsilons_48[i + 8 * j] = epsilon48;
      }
    }
  }
};

}// namespace simulator::physics
