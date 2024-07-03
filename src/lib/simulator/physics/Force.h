#pragma once

#include <array>
#include <cstdlib>

namespace simulator::physics {

template<const size_t DIMENSIONS>
class Force {
 public:
  virtual ~Force() = default;

#pragma omp declare simd simdlen(4) inbranch uniform(this, position1, mass1, type1) linear(ref(mass2, type2))
  virtual void inline calculateForce(
      std::array<double, DIMENSIONS> position1,
      double mass1,
      long type1,
      std::array<double, DIMENSIONS> position2,
      double &mass2,
      long &type2,
      std::array<double, DIMENSIONS> &force,
      std::array<double, DIMENSIONS> &correction) = 0;

  virtual double inline calculate_boundary_force(double diff, int type) = 0;
};

}// namespace simulator::physics