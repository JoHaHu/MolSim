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
auto static calculate_force(VectorizedParticle& p1,
                            VectorizedParticle& p2) -> std::array<double_v, 3> {
  SPDLOG_TRACE("Entering Gravity calculate_force");

  const auto diff = p2.position - p1.position;

  const auto norm = ArrayUtils::L2Norm(diff);

  //  if (norm == 0) {
  //    SPDLOG_DEBUG("Zero distance between particles encountered");
  //  }

  const auto temp = (p1.mass * p2.mass) / (norm * norm * norm);
  auto force = temp * diff;
  SPDLOG_TRACE("Exiting Gravity calculate_force");
  auto norm_mask = norm == 0;

  stdx::where(norm_mask, force[0]) = 0;
  stdx::where(norm_mask, force[1]) = 0;
  stdx::where(norm_mask, force[2]) = 0;

  return force;
};

}// namespace simulator::physics::gravity
