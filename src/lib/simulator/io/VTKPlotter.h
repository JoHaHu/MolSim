#pragma once

#include "lib/ParticleContainer.h"
class VTKPlotter {
 public:
  auto static plotParticles(ParticleContainer &pc, int iteration) -> void;
};
