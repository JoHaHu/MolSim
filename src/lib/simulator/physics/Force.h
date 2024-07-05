#pragma once

#include <array>
#include <cstdlib>

namespace simulator::physics {

class Force {
 public:
  virtual ~Force() = default;

#pragma omp declare simd simdlen(8) inbranch uniform(this, x1, y1, mass1, type1, correction) linear(ref(x2, y2, mass2, type2))
  virtual void inline calculateForce_2D(
      double x1,
      double y1,
      double mass1,
      long type1,
      double &x2,
      double &y2,
      double &mass2,
      long &type2,
      std::array<double, 2> &force,
      std::array<double, 2> &correction) = 0;

#pragma omp declare simd simdlen(8) inbranch uniform(this, x1, y1, z1, mass1, type1, correction) linear(ref(x2, y2, z2, mass2, type2))
  virtual void inline calculateForce_3D(
      double x1,
      double y1,
      double z1,
      double mass1,
      long type1,
      double &x2,
      double &y2,
      double &z2,
      double &mass2,
      long &type2,
      std::array<double, 3> &force,
      std::array<double, 3> &correction) = 0;

  virtual double inline calculate_boundary_force(double diff, int type) = 0;
};

}// namespace simulator::physics