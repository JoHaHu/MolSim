#pragma once

#include <array>
#include <cstdint>
#include <cstdlib>

namespace simulator::physics {

/**
 * @brief Base class for calculating forces in a simulation.
 */
class Force {
 public:
  virtual ~Force() = default;

#pragma omp declare simd simdlen(4) uniform(this, x1, y1, mass1, type1, correction, membrane1) linear(ref(x2, y2, mass2, type2, membrane2))
#pragma omp declare simd simdlen(8) uniform(this, x1, y1, mass1, type1, correction, membrane1) linear(ref(x2, y2, mass2, type2, membrane2))
  /**
   * @brief Calculates 2D force between particles.
   * @param x1 X-coordinate of the first particle.
   * @param y1 Y-coordinate of the first particle.
   * @param mass1 Mass of the first particle.
   * @param type1 Type of the first particle.
   * @param x2 X-coordinate of the second particle.
   * @param y2 Y-coordinate of the second particle.
   * @param mass2 Mass of the second particle.
   * @param type2 Type of the second particle.
   * @param result_x Resultant force in the x-direction.
   * @param result_y Resultant force in the y-direction.
   * @param correction Correction factors.
   */
  virtual void inline calculateForce_2D(
      double const &x1,
      double const &y1,
      double const &mass1,
      long const &type1,
      uint8_t const &membrane1,
      double &x2,
      double &y2,
      double &mass2,
      long &type2,
      uint8_t &membrane2,
      double &result_x,
      double &result_y,
      std::array<double, 2> &correction) = 0;

#pragma omp declare simd simdlen(4) uniform(this, x1, y1, z1, mass1, type1, correction, membrane1) linear(ref(x2, y2, z2, mass2, type2, membrane2))
#pragma omp declare simd simdlen(8) uniform(this, x1, y1, z1, mass1, type1, correction, membrane1) linear(ref(x2, y2, z2, mass2, type2, membrane2))
  /**
   * @brief Calculates 3D force between particles.
   * @param x1 X-coordinate of the first particle.
   * @param y1 Y-coordinate of the first particle.
   * @param z1 Z-coordinate of the first particle.
   * @param mass1 Mass of the first particle.
   * @param type1 Type of the first particle.
   * @param x2 X-coordinate of the second particle.
   * @param y2 Y-coordinate of the second particle.
   * @param z2 Z-coordinate of the second particle.
   * @param mass2 Mass of the second particle.
   * @param type2 Type of the second particle.
   * @param result_x Resultant force in the x-direction.
   * @param result_y Resultant force in the y-direction.
   * @param result_z Resultant force in the z-direction.
   * @param correction Correction factors.
   */
  virtual void inline calculateForce_3D(
      double const &x1,
      double const &y1,
      double const &z1,
      double const &mass1,
      long const &type1,
      uint8_t const &membrane1,
      double &x2,
      double &y2,
      double &z2,
      double &mass2,
      long &type2,
      uint8_t &membrane2,
      double &result_x,
      double &result_y,
      double &result_z,
      std::array<double, 3> &correction) = 0;

  /**
   * @brief Calculates the boundary force.
   * @param diff Difference in position.
   * @param type Type of the particle.
   * @return Calculated boundary force.
   */
  virtual double inline calculate_boundary_force(double diff, int type) = 0;
};

}// namespace simulator::physics
