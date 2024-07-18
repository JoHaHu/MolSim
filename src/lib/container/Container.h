#pragma once

#include "Particle.h"
#include "simulator/physics/MembraneForce.h"
namespace container {

template<const size_t DIMENSIONS>
class Container {
 public:
  virtual ~Container() = default;

  /**
  * @brief Applies a function to each particle in the container.
  * @param f Function to apply to each particle.
  */
  virtual void linear(std::function<void(Particles<DIMENSIONS> &, size_t)>) = 0;

  /**
  * @brief Applies a function to each particle in the container.
  * @param f Function to apply to each particle.
  */
  virtual void membrane(simulator::physics::MembraneForce &) = 0;

  /**
  * @brief Applies a force to each pair of particles in the container.
  *
  * @param f Function to apply to each pair of particles.
  */
  virtual void pairwise(bool parallel, bool vector) = 0;

  /**
   * Swap old and new forces
   * */
  virtual void swap_force() = 0;
  /**
  * @brief Returns the number of particles in the container.
  *
  * @return Size of the container.
  */
  virtual size_t size() = 0;

  /**
   * @brief applies boundary conditions, only effective when using linked cells
   * */
  virtual void boundary() = 0;

  /**
  * @brief Updates internal data structures after position recalculations.
  */
  virtual void refresh() = 0;

  /**
  * @brief Inserts a particle into the container.
  *
  * @param p Particle to insert.
  */
  virtual void insert(Particle<DIMENSIONS> p) = 0;
};

}// namespace container