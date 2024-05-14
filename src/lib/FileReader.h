/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include "Particle.h"
#include "ParticleContainer.h"

#include <vector>

class FileReader {

 public:
  static void readFile(std::vector<Particle> &particles, char *filename);
};
