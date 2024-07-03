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
template<const size_t DIMENSIONS>
class LennardJonesForce final : public Force<DIMENSIONS> {

  double cutoff = 3.0;

  static const long NUM_TYPES = 8;
  // For now supports up to 8 types, can be increased

  alignas(64) std::array<double, NUM_TYPES * NUM_TYPES> epsilons_24;
  alignas(64) std::array<double, NUM_TYPES * NUM_TYPES> epsilons_48;

 public:
  explicit LennardJonesForce(double cutoff, std::vector<double> epsilons, std::vector<double> sigmas) {
    initialize_constants(epsilons, sigmas, cutoff);
  }

#pragma omp declare simd inbranch simdlen(4) uniform(this, position1, mass1, type1) linear(ref(mass2, type2))
  inline void calculateForce(
      std::array<double, DIMENSIONS> position1,
      double mass1,
      long type1,
      std::array<double, DIMENSIONS> position2,
      double &mass2,
      long &type2,
      std::array<double, DIMENSIONS> &force,
      std::array<double, DIMENSIONS> &correction) override {

    const auto diff = (position2 + correction) - position1;
    const auto norm_2 = ArrayUtils::L2NormSquared(diff);
    const auto norm_6 = norm_2 * norm_2 * norm_2;
    auto epsilon24 = epsilons_24[8 * type1 + type2];
    auto epsilon48 = epsilons_48[8 * type1 + type2];
    const auto temp = (norm_6 * epsilon24 - epsilon48) / (norm_6 * norm_6 * norm_2);
    if (norm_2 > cutoff) {
      return;
    }
    force = temp * diff;
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
  auto initialize_constants(std::vector<double> epsilons, std::vector<double> sigmas, double c) {

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
