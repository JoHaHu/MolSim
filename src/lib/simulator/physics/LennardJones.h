#pragma once

#include "lib/Particle.h"

namespace simulator::physics {
class LennardJones {
 public:
  LennardJones() = delete;
  auto static calculateF(Particle const &particle1, Particle const &particle2) -> std::array<double, 3>;
};
}// namespace simulator::physics
