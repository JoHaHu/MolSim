#pragma once

#include "lib/ParticleContainer.h"

namespace simulator::io {

class Plotter {
 public:
  Plotter() = default;
  virtual ~Plotter() = default;

  //! Function for plotting particles
  /*!
  plots the particles of the particle array, takes an integer value and has no return value
  \param iteration an integer argument that sets the number of iterations
*/
  virtual auto plotParticles(ParticleContainer &particle_container, int) -> void = 0;
};

}// namespace simulator::io