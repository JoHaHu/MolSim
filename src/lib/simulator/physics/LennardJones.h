#pragma once

#include "Particle.h"
#include "experimental/simd"
#include "utils/ArrayUtils.h"

namespace simulator::physics::lennard_jones {
/**
 * Implements physics for simulations using lennard-jones as force Model
 * */

inline static double epsilon;
inline static double sigma;
inline static double cutoff;


auto static calculate_force(const VectorizedParticle& p1, const VectorizedParticle& p2) -> std::array<double_v, 3> {
  SPDLOG_TRACE("Entering LennardJones calculate_force");

  const auto diff = p2.position - p1.position;

  const auto norm = ArrayUtils::L2Norm(diff);

  // TODO precalculate sigma^3 , 24 * epsilon

  const auto sigma_over_norm_3 = (sigma / norm) * (sigma / norm) * (sigma / norm);
  const auto sigma_over_norm_6 = sigma_over_norm_3 * sigma_over_norm_3;

  const auto temp = 24 * epsilon / (norm * norm) * (1 - 2 * sigma_over_norm_6) * sigma_over_norm_6;
  auto force = temp * diff;

  SPDLOG_TRACE("Exiting LennardJones calculate_force");

  auto norm_mask = norm > cutoff;
  stdx::where(norm_mask, force[0]) = 0;
  stdx::where(norm_mask, force[1]) = 0;
  stdx::where(norm_mask, force[2]) = 0;

  return force;
}

}// namespace simulator::physics::lennard_jones
