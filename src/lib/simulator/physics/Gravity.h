#pragma once

#include "lib/Particle.h"
class Gravity {

 public:
  auto static calculateF(Particle const &p1, Particle const &p2) -> std::array<double, 3>;
};