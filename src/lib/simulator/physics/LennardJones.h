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
  epsilon48 = 48 * epsilon * sigma6 * sigma6;
}

auto static calculate_force_vectorized(const VectorizedParticle &p1, const VectorizedParticle &p2, double_mask mask, std::array<double_v, 3> &force) {
  SPDLOG_TRACE("Entering LennardJones calculate_force_vectorized");

  const auto diff = p2.position - p1.position;

  const auto norm_2 = ArrayUtils::L2NormSquared(diff);

  const auto norm_6 = norm_2 * norm_2 * norm_2;

  const auto temp = (norm_6 * epsilon24 - epsilon48) / (norm_6 * norm_6 * norm_2);
  force = temp * diff;

  if (stdx::any_of(norm_2 == 0 && mask)) [[unlikely]] {
    SPDLOG_WARN("zero reflecting_distance between particle");
  }
  auto norm_mask = norm_2 > cutoff;

  stdx::where(norm_mask, force[0]) = 0;
  stdx::where(norm_mask, force[1]) = 0;
  stdx::where(norm_mask, force[2]) = 0;

  SPDLOG_TRACE("Exiting LennardJones calculate_force_vectorized");
}

auto static calculate_force(Particles &p, size_t index, std::array<double, 3>& position2) {
  SPDLOG_TRACE("Entering LennardJones calculate_force_vectorized");

  const auto diff = position2 - std::array<double, 3>({p.position_x[index], p.position_y[index], p.position_z[index]});

  const auto norm_2 = ArrayUtils::L2NormSquared(diff);

  const auto norm_6 = norm_2 * norm_2 * norm_2;

  const auto temp = (norm_6 * epsilon24 - epsilon48) / (norm_6 * norm_6 * norm_2);
  auto force = temp * diff;

  if (norm_2 == 0) [[unlikely]] {
    SPDLOG_WARN("zero distance between particle");
  }

  if (norm_2 > cutoff) {
    force[0] = 0;
    force[1] = 0;
    force[2] = 0;
  }

  p.force_x[index] += force[0];
  p.force_y[index] += force[1];
  p.force_z[index] += force[2];

  SPDLOG_TRACE("Exiting LennardJones calculate_force_vectorized");
}

}// namespace simulator::physics::lennard_jones
