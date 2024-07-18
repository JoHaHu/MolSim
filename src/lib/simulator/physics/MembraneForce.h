#pragma once

#include "Particle.h"
#include "experimental/simd"
#include "spdlog/spdlog.h"
#include "utils/ArrayUtils.h"

#include "Force.h"
#include "immintrin.h"
namespace simulator::physics {

/**
 * Implements physics for simulations using lennard-jones as force Model
 * */

class MembraneForce {
 private:
  double k;
  double rzero;

 public:
  MembraneForce(double k, double rzero) : k(k), rzero(rzero) {}

  void apply_force2D(
      double const &x1,
      double const &y1,
      double &x2,
      double &y2,
      double &result_x,
      double &result_y,
      bool diagonal,
      std::array<double, 2> &correction) {

    const auto diff = (std::array<double, 2>({x2, y2}) + correction) - std::array<double, 2>({x1, y1});
    const auto norm = ArrayUtils::L2Norm(diff);

    std::array<double, 2> result = {0, 0};
    if (diagonal) {
      result = k * (1 - std::numbers::sqrt2 * rzero / norm) * diff;
    } else {
      result = k * (1 - rzero / norm) * diff;
    }
    result_x = result[0];
    result_y = result[1];
  }

  void apply_force3D(
      double const &x1,
      double const &y1,
      double const &z1,
      double &x2,
      double &y2,
      double &z2,
      double &result_x,
      double &result_y,
      double &result_z,
      bool diagonal,
      std::array<double, 3> &correction) {

    const auto diff = (std::array<double, 3>({x2, y2, z2}) + correction) - std::array<double, 3>({x1, y1, z1});
    const auto norm = ArrayUtils::L2Norm(diff);

    std::array<double, 3> result = {0, 0};
    if (diagonal) {
      result = k * (1 - std::numbers::sqrt2 * rzero / norm) * diff;
    } else {
      result = k * (1 - rzero / norm) * diff;
    }
    result_x = result[0];
    result_y = result[1];
    result_z = result[2];
  }
};

}// namespace simulator::physics
