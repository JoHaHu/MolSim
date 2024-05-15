#pragma once

#include "lib/Particle.h"
#include "lib/ParticleContainer.h"

#include <vector>

namespace simulator::io {

class FileReader {

 public:
  static void read_file(std::vector<Particle> &particles, std::string filename);
};
}// namespace simulator::io