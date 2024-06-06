#pragma once

#include "Particle.h"
namespace simulator::physics {

/**
 * Implements physics for simulations using gravity as force Model
 *
 * */
class Gravity {
 public:
  auto calculate_force(const Particle &particle1, const Particle &particle2) -> std::array<double, 3>;
};

}// namespace simulator::physics
