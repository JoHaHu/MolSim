#pragma once

#include "Particle.h"
#include "experimental/simd"
#include "spdlog/spdlog.h"
#include "utils/ArrayUtils.h"

namespace simulator::physics::lennard_jones {
/**
 * Implements physics for simulations using lennard-jones as force Model
 * */

inline static double epsilon = 5.0;
inline static double sigma = 1.0;
inline static double cutoff = 3.0;

static double epsilon24;
static double epsilon48;

auto static initialize_constants(double e, double s, double c) {
  epsilon = e;
  sigma = s;
  cutoff = c * c;

  auto sigma3 = sigma * sigma * sigma;
  auto sigma6 = sigma3 * sigma3;
  epsilon24 = 24 * epsilon * sigma6;
  epsilon48 = 2 * epsilon24 * sigma6;
}
template<const size_t DIMENSIONS>
auto static calculate_force_vectorized(const VectorizedParticle<DIMENSIONS> &p1, const VectorizedParticle<DIMENSIONS> &p2, double_mask mask, std::array<double_v, DIMENSIONS> &force, std::array<double_v, DIMENSIONS> &corrections) {
  SPDLOG_TRACE("Entering LennardJones calculate_force_vectorized");

  const auto diff = (p2.position + corrections) - p1.position;

  const auto norm_2 = ArrayUtils::L2NormSquared(diff);

  const auto norm_6 = norm_2 * norm_2 * norm_2;

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
 * simplified jennard jones force for boundary particle. since the vector is perpendicular to the boundary plane, only one position component is relevant to calculate the norm
 * */
auto static calculate_force(double diff) -> double {
  const auto norm_2 = diff * diff;
  const auto norm_6 = norm_2 * norm_2 * norm_2;
  const auto temp = (norm_6 * epsilon24 - epsilon48) / (norm_6 * norm_6 * norm_2);
  return temp * diff;
}

}// namespace simulator::physics::lennard_jones
