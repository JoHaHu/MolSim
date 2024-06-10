#pragma once

#include "Particle.h"
#include "experimental/simd"
#include "utils/ArrayUtils.h"

namespace stdx {
using namespace std::experimental;
using namespace std::experimental::__proposed;
}// namespace stdx

namespace simulator::physics::lennard_jones {
/**
 * Implements physics for simulations using lennard-jones as force Model
 * */

using double_v = stdx::native_simd<double>;
using double_mask = stdx::native_simd_mask<double>;
using int_v = stdx::native_simd<int>;

inline static double epsilon;
inline static double sigma;
inline static double cutoff;

auto static calculate_force(const double_v &position_x, const double_v &position_y, const double_v &position_z, const double_v &mass, const int_v &type,
                            const double_v &position_x2, const double_v &position_y2, const double_v &position_z2, const double_v &mass2, const int_v &type2) -> std::array<double_v, 3> {
  SPDLOG_TRACE("Entering LennardJones calculate_force");

  const auto x_diff = position_x2 - position_x;
  const auto y_diff = position_y2 - position_y;
  const auto z_diff = position_z2 - position_z;

  const auto diff = {x_diff, y_diff, z_diff};

  const auto norm = stdx::sqrt(std::accumulate(std::cbegin(diff), std::cend(diff), double_v(0.0), [](auto a, auto b) { return a + b * b; }));

  if (stdx::any_of(norm == 0)) {
    SPDLOG_DEBUG("Zero distance between particles encountered. This should not be possible.");
  }

  const auto sigma_over_norm_3 = (sigma / norm) * (sigma / norm) * (sigma / norm);
  const auto sigma_over_norm_6 = sigma_over_norm_3 * sigma_over_norm_3;

  const auto temp = 24 * epsilon / (norm * norm) * (1 - 2 * sigma_over_norm_6) * sigma_over_norm_6;
  auto force_x = temp * x_diff;
  auto force_y = temp * y_diff;
  auto force_z = temp * z_diff;
  SPDLOG_TRACE("Exiting LennardJones calculate_force");

  auto norm_mask = norm > cutoff;
  stdx::where(norm_mask, force_x) = 0;
  stdx::where(norm_mask, force_y) = 0;
  stdx::where(norm_mask, force_z) = 0;

  return {force_x, force_y, force_z};
}

}// namespace simulator::physics::lennard_jones
