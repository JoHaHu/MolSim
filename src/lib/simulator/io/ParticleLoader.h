#pragma once

#include "lib/ParticleContainer.h"

namespace simulator::io {
class ParticleLoader {

 public:
  ParticleLoader() = default;
  virtual ~ParticleLoader() = default;

  virtual auto load_particles() -> ParticleContainer = 0;
};

}// namespace simulator::io
