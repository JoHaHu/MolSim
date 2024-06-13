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
auto static calculate_force_vectorized(VectorizedParticle &p1,
                                       VectorizedParticle &p2,
                                       double_mask mask,
                                       std::array<double_v, 3> &force) {
  SPDLOG_TRACE("Entering Gravity calculate_force_vectorized");

  const auto diff = p2.position - p1.position;

  const auto norm = ArrayUtils::L2Norm(diff);

  //  if (norm == 0) {
  //    SPDLOG_DEBUG("Zero distance between particles encountered");
  //  }

  const auto temp = (p1.mass * p2.mass) / (norm * norm * norm);
  force = temp * diff;
  SPDLOG_TRACE("Exiting Gravity calculate_force_vectorized");

  auto norm_mask = norm == 0;
  stdx::where(norm_mask, force[0]) = 0;
  stdx::where(norm_mask, force[1]) = 0;
  stdx::where(norm_mask, force[2]) = 0;
};
// TODO boundary condition gravity
auto static calculate_force(Particles &p, size_t index, std::array<double, 3> &position2) {
  SPDLOG_TRACE("Entering Gravity calculate_force_vectorized");

  const auto diff = position2 - std::array<double, 3>({p.position_x[index], p.position_y[index], p.position_z[index]});

  const auto norm = ArrayUtils::L2Norm(diff);

  const auto temp = (p.mass[index] * p.mass[index]) / (norm * norm * norm);
  auto force = temp * diff;
  SPDLOG_TRACE("Exiting Gravity calculate_force_vectorized");

  if (norm == 0) {
    force[0] = 0;
    force[1] = 0;
    force[2] = 0;
  }
  p.force_x[index] = force[0];
  p.force_y[index] = force[1];
  p.force_z[index] = force[2];
};

}// namespace simulator::physics::gravity
