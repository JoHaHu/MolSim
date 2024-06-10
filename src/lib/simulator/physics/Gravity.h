#pragma once

#include "Particle.h"
#include "experimental/simd"
#include "spdlog/spdlog.h"
#include "utils/ArrayUtils.h"

namespace stdx {
using namespace std::experimental;
using namespace std::experimental::__proposed;
}
namespace simulator::physics::gravity {

using double_v = stdx::native_simd<double>;
using int_v = stdx::native_simd<int>;

/**
 * Implements physics for simulations using gravity as force Model
 *
 * */
auto static calculate_force(const double_v &position_x, const double_v &position_y, const double_v &position_z, const double_v &mass, const int_v &type,
                            const double_v &position_x2, const double_v &position_y2, const double_v &position_z2, const double_v &mass2, const int_v &type2) -> std::array<double_v, 3> {
  SPDLOG_TRACE("Entering Gravity calculate_force");

  const auto x_diff = position_x2 - position_x;
  const auto y_diff = position_y2 - position_y;
  const auto z_diff = position_z2 - position_z;

  const auto diff = {x_diff, y_diff, z_diff};

  const auto norm = stdx::sqrt(std::accumulate(std::cbegin(diff), std::cend(diff), double_v(0.0), [](auto a, auto b) { return a + b * b; }));

  //  if (norm == 0) {
  //    SPDLOG_DEBUG("Zero distance between particles encountered");
  //  }

  const auto temp = (mass * mass2) / (norm * norm * norm);
  auto force_x = temp * x_diff;
  auto force_y = temp * y_diff;
  auto force_z = temp * z_diff;
  SPDLOG_TRACE("Exiting Gravity calculate_force");
  auto norm_mask = norm == 0;

  stdx::where(norm_mask, force_x) = 0;
  stdx::where(norm_mask, force_y) = 0;
  stdx::where(norm_mask, force_z) = 0;

  return {force_x, force_y, force_z};
};

}// namespace simulator::physics::gravity
