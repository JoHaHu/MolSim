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
#pragma omp declare simd simdlen(4) uniform(this, x1, y1, mass1, type1, correction, membrane1) linear(ref(x2, y2, mass2, type2, membrane2))
#pragma omp declare simd simdlen(8) uniform(this, x1, y1, mass1, type1, correction, membrane1) linear(ref(x2, y2, mass2, type2, membrane2))
  inline void calculateForce_2D(
      double const &x1,
      double const &y1,
      double const &mass1,
      long const &type1,
      uint8_t const &membrane1,
      double &x2,
      double &y2,
      double &mass2,
      long &type2,
      uint8_t &membrane2,
      double &result_x,
      double &result_y,
      std::array<double, 2> &correction) override {

    SPDLOG_TRACE("Entering Gravity calculate_force_vectorized");

    const auto diff = (std::array<double, 2>({x2, y2}) + correction) - std::array<double, 2>({x1, y1});

    const auto norm = ArrayUtils::L2Norm(diff);

    const auto temp = (mass1 * mass2) / (norm * norm * norm);
    auto force = temp * diff;
    result_x = force[0];
    result_y = force[1];
    SPDLOG_TRACE("Exiting Gravity calculate_force_vectorized");
  };

#pragma omp declare simd simdlen(4) uniform(this, x1, y1, mass1, type1, correction, membrane1) linear(ref(x2, y2, mass2, type2, membrane2))
#pragma omp declare simd simdlen(8) uniform(this, x1, y1, mass1, type1, correction, membrane1) linear(ref(x2, y2, mass2, type2, membrane2))
  inline void calculateForce_3D(
      double const &x1,
      double const &y1,
      double const &z1,
      double const &mass1,
      long const &type1,
      uint8_t const &membrane1,
      double &x2,
      double &y2,
      double &z2,
      double &mass2,
      long &type2,
      uint8_t &membrane2,
      double &result_x,
      double &result_y,
      double &result_z,
      std::array<double, 3> &correction) override {

    SPDLOG_TRACE("Entering Gravity calculate_force_vectorized");
    const auto diff = (std::array<double, 3>({x2, y2, z2}) + correction) - std::array<double, 3>({x1, y1, z1});

    const auto norm = ArrayUtils::L2Norm(diff);

    const auto temp = (mass1 * mass2) / (norm * norm * norm);
    auto force = temp * diff;
    result_x = force[0];
    result_y = force[1];
    result_z = force[2];
    SPDLOG_TRACE("Exiting Gravity calculate_force_vectorized");
  };

  double calculate_boundary_force(double diff, int type) override {
    return 0;
  }
};

/**
 * use for applying constant gravity
 * */
auto static constant_gravity(double mass, double constant) {
  return mass * constant;
}

}// namespace simulator::physics
