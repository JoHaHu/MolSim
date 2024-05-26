#pragma once

#include "Particle.h"

namespace simulator::physics {
/**
 * Implements physics for simulations using lennard-jones as force Model
 * */
class LennardJones {
 public:
  auto calculate_force(const Particle &particle1, const Particle &particle2) -> std::array<double, 3>;

};
}// namespace simulator::physics
