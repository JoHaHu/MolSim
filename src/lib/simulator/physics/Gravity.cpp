#include "Gravity.h"
#include "lib/utils/ArrayUtils.h"
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
namespace simulator::physics {

auto Gravity::calculate_force(const Particle &particle1, const Particle &particle2) -> std::array<double, 3> {
  spdlog::trace("Entering calculate_force");

  spdlog::trace("Particle 1: mass = {}, position = ({}, {}, {})", particle1.mass, particle1.position[0], particle1.position[1], particle1.position[2]);
  spdlog::trace("Particle 2: mass = {}, position = ({}, {}, {})", particle2.mass, particle2.position[0], particle2.position[1], particle2.position[2]);

  const auto x_diff = particle2.position - particle1.position;

  spdlog::trace("Position difference: ({}, {}, {})", x_diff[0], x_diff[1], x_diff[2]);

  auto norm = ArrayUtils::L2Norm(x_diff);
  spdlog::trace("Norm of position difference: {}", norm);

  if (norm == 0) {
    spdlog::warn("Zero distance between particles encountered");
  }

  const auto f = (particle1.mass * particle2.mass) / pow(norm, 3) * x_diff;
  spdlog::trace("Calculated force: ({}, {}, {})", f[0], f[1], f[2]);

  spdlog::trace("Exiting calculate_force");

  return f;
}
}// namespace simulator::physics