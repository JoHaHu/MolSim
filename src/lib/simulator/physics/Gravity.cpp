//
// Created by johannes on 14.05.24.
//

#include "Gravity.h"
#include "lib/utils/ArrayUtils.h"
auto Gravity::calculateF(const Particle &p1, const Particle &p2) -> std::array<double, 3> {
  const Container auto x_diff = p2.x - p1.x;

  auto norm = ArrayUtils::L2Norm(x_diff);
  const Container auto f = (p1.m * p2.m) / pow(norm, 3) * x_diff;
  return f;
}
