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

static double epsilon24;
static double epsilon48;

auto static initialize_constants(double e, double s, double c) {
  epsilon = e;
  sigma = s;
  cutoff = c;

  auto sigma3 = sigma * sigma * sigma;
  auto sigma6 = sigma3 * sigma3;
  epsilon24 = 24 * epsilon * sigma6;
  epsilon48 = 48 * epsilon * sigma6 * sigma6;
}

auto static calculate_force(const VectorizedParticle &p1, const VectorizedParticle &p2) -> std::array<double_v, 3> {
  SPDLOG_TRACE("Entering LennardJones calculate_force");

  const auto diff = p2.position - p1.position;

  const auto norm = ArrayUtils::L2Norm(diff);

  const auto norm_2 = norm * norm;
  const auto norm_6 = norm_2 * norm_2 * norm_2;

  const auto temp = (norm_6 * epsilon24 - epsilon48) / (norm_6 * norm_6 * norm_2);
  auto force = temp * diff;

  SPDLOG_TRACE("Exiting LennardJones calculate_force");

  auto norm_mask = norm > cutoff;
  stdx::where(norm_mask, force[0]) = 0;
  stdx::where(norm_mask, force[1]) = 0;
  stdx::where(norm_mask, force[2]) = 0;

  return force;
}

}// namespace simulator::physics::lennard_jones
