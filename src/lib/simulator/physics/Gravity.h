#pragma once

#include "Force.h"
#include "Particle.h"
#include "experimental/simd"
#include "spdlog/spdlog.h"
#include "utils/ArrayUtils.h"

namespace simulator::physics {

template<const size_t DIMENSIONS>
class Gravity final : public Force<DIMENSIONS> {

 public:
/**
 * Implements physics for simulations using gravity as force Model
 *
 * */
#pragma omp declare simd inbranch uniform(this, position1, mass1, type1) linear(ref(position2, mass2, type2))
  void inline calculateForce(
      std::array<double, DIMENSIONS> position1,
      double mass1,
      long type1,
      std::array<double, DIMENSIONS> position2,
      double &mass2,
      long &type2,
      std::array<double, DIMENSIONS> &force,
      std::array<double, DIMENSIONS> &correction) override {

    SPDLOG_TRACE("Entering Gravity calculate_force_vectorized");

    const auto diff = (position2 - correction) - position2;

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
