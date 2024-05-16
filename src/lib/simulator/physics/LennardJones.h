#pragma once

#include "Particle.h"

namespace simulator::physics {
class LennardJones {
 public:
  LennardJones() = delete;
  auto static calculate_force(const Particle &particle1, const Particle &particle2) -> std::array<double, 3>;
};
}// namespace simulator::physics
