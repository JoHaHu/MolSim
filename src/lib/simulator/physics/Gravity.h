#pragma once

#include "Particle.h"
#include "experimental/simd"
#include "spdlog/spdlog.h"
#include "utils/ArrayUtils.h"

namespace simulator::physics::gravity {

/**
 * Implements physics for simulations using gravity as force Model
 *
 * */
template<const size_t DIMENSIONS>
__attribute__((__always_inline__))
auto static calculate_force_vectorized(VectorizedParticle<DIMENSIONS> &p1,
                                       VectorizedParticle<DIMENSIONS> &p2,
                                       double_mask mask,
                                       std::array<double_v, DIMENSIONS> &force,
                                       std::array<double_v, DIMENSIONS> &corrections) {
  SPDLOG_TRACE("Entering Gravity calculate_force_vectorized");

  const auto diff = (p2.position - corrections) - p1.position;

  const auto norm = ArrayUtils::L2Norm(diff);

  //  if (norm == 0) {
  //    SPDLOG_DEBUG("Zero distance between particles encountered");
  //  }

  const auto temp = (p1.mass * p2.mass) / (norm * norm * norm);
  force = temp * diff;
  SPDLOG_TRACE("Exiting Gravity calculate_force_vectorized");

  auto norm_mask = norm == 0;
  for (int i = 0; i < DIMENSIONS; ++i) {
    stdx::where(norm_mask, force[i]) = 0;
  }
};

}// namespace simulator::physics::gravity
