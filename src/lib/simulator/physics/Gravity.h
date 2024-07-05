#pragma once

#include "Force.h"
#include "Particle.h"
#include "experimental/simd"
#include "spdlog/spdlog.h"
#include "utils/ArrayUtils.h"

namespace simulator::physics {

class Gravity final : public Force {

 public:
/**
 * Implements physics for simulations using gravity as force Model
 *
 * */
#pragma omp declare simd inbranch simdlen(8) uniform(this, x1, y1, mass1, type1) linear(ref(x2, y2, mass2, type2))
  inline void calculateForce_2D(
      double x1,
      double y1,
      double mass1,
      long type1,
      double &x2,
      double &y2,
      double &mass2,
      long &type2,
      std::array<double, 2> &force,
      std::array<double, 2> &correction) override {

    SPDLOG_TRACE("Entering Gravity calculate_force_vectorized");

    const auto diff = (std::array<double, 2>({x2, y2}) + correction) - std::array<double, 2>({x1, y1});

    const auto norm = ArrayUtils::L2Norm(diff);

    const auto temp = (mass1 * mass2) / (norm * norm * norm);
    force = temp * diff;
    SPDLOG_TRACE("Exiting Gravity calculate_force_vectorized");
  };

#pragma omp declare simd inbranch simdlen(8) uniform(this, x1, y1, z1, mass1, type1) linear(ref(x2, y2, z2, mass2, type2))
  inline void calculateForce_3D(
      double x1,
      double y1,
      double z1,
      double mass1,
      long type1,
      double &x2,
      double &y2,
      double &z2,
      double &mass2,
      long &type2,
      std::array<double, 3> &force,
      std::array<double, 3> &correction) override {

    SPDLOG_TRACE("Entering Gravity calculate_force_vectorized");
    const auto diff = (std::array<double, 3>({x2, y2, z2}) + correction) - std::array<double, 3>({x1, y1, z1});

    const auto norm = ArrayUtils::L2Norm(diff);

    const auto temp = (mass1 * mass2) / (norm * norm * norm);
    force = temp * diff;
    SPDLOG_TRACE("Exiting Gravity calculate_force_vectorized");
  };

  double calculate_boundary_force(double diff, int type) override {
    return 0;
  }
};

/**
 * use for apply constant gravity
 * */
auto static constant_gravity(double mass, double constant) {
  return mass * constant;
}

}// namespace simulator::physics
