#pragma once

#include <iterator>

namespace simulator::io {

/*!
 * a abstract for different plotting methods
 */
template<std::input_iterator I>
  requires(std::same_as<typename std::iterator_traits<I>::value_type, Particle>)
class Plotter {
 public:
  Plotter() = default;
  virtual ~Plotter() = default;

  //! Function for plotting particles
  /*!
  plots the particles of the particle array, takes an integer value and has no return value
  \param iteration an integer argument that sets the number of iterations
 */
  virtual auto plotParticles(I &particle_container, int) -> void = 0;
};

}// namespace simulator::io