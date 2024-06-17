#pragma once

#include "Particle.h"
#include "experimental/simd"
#include "spdlog/spdlog.h"
#include "utils/ArrayUtils.h"

#include "immintrin.h"
namespace simulator::physics::lennard_jones {
/**
 * Implements physics for simulations using lennard-jones as force Model
 * */

inline static double cutoff = 3.0;

static const int NUM_TYPES = 8;
// For now supports up to 8 types, can be increased
alignas(64) inline static std::array<double, NUM_TYPES * NUM_TYPES> epsilons_24;
alignas(64) inline static std::array<double, NUM_TYPES * NUM_TYPES> epsilons_48;

auto static initialize_constants(std::vector<double> epsilons, std::vector<double> sigmas, double c) {

  cutoff = c * c;

  for (int i = 0; i < epsilons.size(); ++i) {
    for (int j = 0; j < epsilons.size(); ++j) {
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
template<const size_t DIMENSIONS>
__attribute__((__always_inline__)) auto static calculate_force_vectorized(const VectorizedParticle<DIMENSIONS> &p1, const VectorizedParticle<DIMENSIONS> &p2, double_mask mask, std::array<double_v, DIMENSIONS> &force, std::array<double_v, DIMENSIONS> &corrections) {
  SPDLOG_TRACE("Entering LennardJones calculate_force_vectorized");

  const auto diff = (p2.position + corrections) - p1.position;

  const auto norm_2 = ArrayUtils::L2NormSquared(diff);

  const auto norm_6 = norm_2 * norm_2 * norm_2;

  long_v vindex = long_v(0);

  stdx::where(vindex < NUM_TYPES, vindex) = 8 * p1.type + p2.type;

  double_v epsilon24;
  double_v epsilon48;


  // Faster than vgatherpd
  for (int i = 0; i < double_v::size(); ++i) {
    epsilon24[i] = epsilons_24[vindex[i]];
    epsilon48[i] = epsilons_48[vindex[i]];
  }

  const auto temp = (norm_6 * epsilon24 - epsilon48) / (norm_6 * norm_6 * norm_2);
  force = temp * diff;

  if (stdx::any_of(norm_2 == 0 && mask)) [[unlikely]] {
    SPDLOG_WARN("zero reflecting_distance between particle");
  }

  auto norm_mask = norm_2 > cutoff;

  for (int i = 0; i < DIMENSIONS; ++i) {
    stdx::where(norm_mask, force[i]) = 0;
  }

  SPDLOG_TRACE("Exiting LennardJones calculate_force_vectorized");
}

/**
 * simplified Lennard-Jones force for boundary particle. since the vector is perpendicular to the boundary plane, only one position component is relevant to calculate the norm
 * */
__attribute__((__always_inline__)) auto static calculate_force(double diff, int type) -> double {
  const auto norm_2 = diff * diff;
  const auto norm_6 = norm_2 * norm_2 * norm_2;
  const auto temp = (norm_6 * epsilons_24[8 * type + type] - epsilons_48[8 * type + type]) / (norm_6 * norm_6 * norm_2);
  return temp * diff;
}

}// namespace simulator::physics::lennard_jones
