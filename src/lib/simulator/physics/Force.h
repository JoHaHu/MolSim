#pragma once

#include <array>
#include <cstdlib>

namespace simulator::physics {

class Force {
 public:
  virtual ~Force() = default;

#pragma omp declare simd simdlen(4) uniform(this, x1, y1, mass1, type1, correction /*,membrane1*/) linear(ref(x2, y2, mass2, type2, membrane2))
#pragma omp declare simd simdlen(8) uniform(this, x1, y1, mass1, type1, correction, membrane1) linear(ref(x2, y2, mass2, type2, membrane2))
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

  virtual double inline calculate_boundary_force(double diff, int type) = 0;
};

}// namespace simulator::physics