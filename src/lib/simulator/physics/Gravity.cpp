//
// Created by johannes on 14.05.24.
//

#include "Gravity.h"
#include "lib/utils/ArrayUtils.h"
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

auto Gravity::calculateF(const Particle &p1, const Particle &p2) -> std::array<double, 3> {
  spdlog::trace("Entering calculateF");

  spdlog::trace("Particle 1: mass = {}, position = ({}, {}, {})", p1.m, p1.x[0], p1.x[1], p1.x[2]);
  spdlog::trace("Particle 2: mass = {}, position = ({}, {}, {})", p2.m, p2.x[0], p2.x[1], p2.x[2]);

  const Container auto x_diff = p2.x - p1.x;

  spdlog::trace("Position difference: ({}, {}, {})", x_diff[0], x_diff[1], x_diff[2]);

  auto norm = ArrayUtils::L2Norm(x_diff);
  spdlog::trace("Norm of position difference: {}", norm);

  if (norm == 0) {
    spdlog::warn("Zero distance between particles encountered");
  }

  const Container auto f = (p1.m * p2.m) / pow(norm, 3) * x_diff;
  spdlog::trace("Calculated force: ({}, {}, {})", f[0], f[1], f[2]);

  spdlog::trace("Exiting calculateF");

  return f;
}

