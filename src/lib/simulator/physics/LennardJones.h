#pragma once

#include "Particle.h"
#include "experimental/simd"
#include "spdlog/spdlog.h"
#include "utils/ArrayUtils.h"

#include "Force.h"
#include "immintrin.h"
namespace simulator::physics {

/**
 * Implements physics for simulations using lennard-jones as force Model
 * */

class LennardJonesForce final : public Force {

  double cutoff = 3.0;

  static const long NUM_TYPES = 8;
  // For now supports up to 8 types, can be increased

  alignas(64) std::array<double, NUM_TYPES *NUM_TYPES> epsilons_24 = std::array<double, NUM_TYPES * NUM_TYPES>();
  alignas(64) std::array<double, NUM_TYPES *NUM_TYPES> epsilons_48 = std::array<double, NUM_TYPES * NUM_TYPES>();

 public:
  explicit LennardJonesForce(double cutoff, std::vector<double> epsilons, std::vector<double> sigmas) {
    initialize_constants(epsilons, sigmas, cutoff);
  }

#pragma omp declare simd simdlen(4) uniform(this, x1, y1, mass1, type1, correction) linear(ref(x2, y2, mass2, type2))
#pragma omp declare simd simdlen(8) uniform(this, x1, y1, mass1, type1, correction) linear(ref(x2, y2, mass2, type2))
  inline void calculateForce_2D(
      double const &x1,
      double const &y1,
      double const &mass1,
      long const &type1,
      double &x2,
      double &y2,
      double &mass2,
      long &type2,
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
    if (norm_2 <= cutoff) {
      result_x = force[0];
      result_y = force[1];
    }
  }

#pragma omp declare simd simdlen(4) uniform(this, x1, y1, z1, mass1, type1, correction) linear(ref(x2, y2, z2, mass2, type2))
#pragma omp declare simd simdlen(8) uniform(this, x1, y1, z1, mass1, type1, correction) linear(ref(x2, y2, z2, mass2, type2))
  inline void calculateForce_3D(
      double const &x1,
      double const &y1,
      double const &z1,
      double const &mass1,
      long const &type1,
      double &x2,
      double &y2,
      double &z2,
      double &mass2,
      long &type2,
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

    if (norm_2 <= cutoff) {
      result_x = force[0];
      result_y = force[1];
      result_z = force[2];
    }
  }

  /**
 * simplified Lennard-Jones force for boundary particle. since the vector is perpendicular to the boundary plane, only one position component is relevant to calculate the norm
 * */
  double inline calculate_boundary_force(double diff, int type) override {
    const auto norm_2 = diff * diff;
    const auto norm_6 = norm_2 * norm_2 * norm_2;
    const auto temp = (norm_6 * epsilons_24[8 * type + type] - epsilons_48[8 * type + type]) / (norm_6 * norm_6 * norm_2);
    return temp * diff;
  }

 private:
  auto initialize_constants(std::vector<double> epsilons, std::vector<double> sigmas, double c) -> void {

    cutoff = c * c;

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
