#pragma once

#include "Particle.h"
#include "ParticleContainer.h"

#include <vector>

class FileReader {

 public:
  static void read_file(std::vector<Particle> &particles, std::string filename);
};
