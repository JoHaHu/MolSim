#pragma once

#include "Particle.h"

namespace simulator::physics {
/**
 * Implements physics for simulations using lennard-jones as force Model
 * */
class LennardJones {
 public:
  LennardJones() = default;
  explicit LennardJones(double cutoff, double sigma = 1.0, double epsilon = 5.0);
  auto calculate_force(const Particle &particle1, const Particle &particle2) -> std::array<double, 3>;

 private:
  double epsilon;
  double sigma;
  double cutoff = 3.0;
};
}// namespace simulator::physics
