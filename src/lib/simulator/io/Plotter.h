#pragma once

#include <iterator>

namespace simulator::io {

/*!
 * a abstract for different plotting methods
 */
template<const size_t DIMENSIONS>
class Plotter {
 public:
  Plotter() = default;
  virtual ~Plotter() = default;

  /**
   * Function for plotting particles
   * plots the particles of the particle container, takes an integer value and has no return value
   * \param iteration an integer argument that sets the number of iterations
   * */

  virtual auto plotParticles(container::ParticleContainer<DIMENSIONS> &particles, int) -> void = 0;
};

}// namespace simulator::io