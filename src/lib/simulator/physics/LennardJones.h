#pragma once

#include "Particle.h"

namespace simulator::physics {
/**
 * Implements physics for simulations using lennard-jones as force Model
 * */
class LennardJones {
 public:
  explicit LennardJones(double cutoff = 3.0, double sigma = 1.0, double epsilon = 5.0);
  auto calculate_force(const double &position_x, const double &position_y, const double &position_z, const double &mass, const int &type,
                       const double &position_x2, const double &position_y2, const double &position_z2, const double &mass2, const int &type2) -> std::tuple<double, double, double>;

 private:
  double epsilon;
  double sigma;
  double cutoff;
};
}// namespace simulator::physics
