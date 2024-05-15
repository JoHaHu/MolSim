#pragma once

#include "ParticleLoader.h"
#include "lib/Particle.h"
#include "lib/ParticleContainer.h"
#include "lib/config/config.h"

#include <vector>

namespace simulator::io {

class FileReader : public ParticleLoader {
 public:
  explicit FileReader(const std::shared_ptr<config::Config> &config);

 private:
  std::shared_ptr<config::Config> config;

 public:
  auto load_particles() -> ParticleContainer override;
};
}// namespace simulator::io