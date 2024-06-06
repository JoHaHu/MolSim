#include "Gravity.h"
#include "utils/ArrayUtils.h"
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
namespace simulator::physics {

/**
 * @brief Calculates the gravitational force between two particles.
 *
 * Computes the gravitational force based on the masses and positions of two particles.
 *
 * @param particle1 The first particle.
 * @param particle2 The second particle.
 * @return std::array<double, 3> The calculated force vector.
 */
auto Gravity::calculate_force(const Particle &particle1, const Particle &particle2) -> std::array<double, 3> {
  SPDLOG_TRACE("Entering Gravity calculate_force");

  SPDLOG_TRACE("Particle 1: mass = {}, position = ({}, {}, {})", particle1.mass, particle1.position[0], particle1.position[1], particle1.position[2]);
  SPDLOG_TRACE("Particle 2: mass = {}, position = ({}, {}, {})", particle2.mass, particle2.position[0], particle2.position[1], particle2.position[2]);

  const auto x_diff = particle2.position - particle1.position;

  SPDLOG_TRACE("Position difference: ({}, {}, {})", x_diff[0], x_diff[1], x_diff[2]);

  const auto norm = ArrayUtils::L2Norm(x_diff);
  SPDLOG_TRACE("Norm of position difference: {}", norm);

  if (norm == 0) {
    SPDLOG_WARN("Zero distance between particles encountered");
    return {0.0, 0.0, 0.0};
  }

  const auto force = (particle1.mass * particle2.mass) / (norm * norm * norm) * x_diff;
  SPDLOG_TRACE("Calculated force: ({}, {}, {})", force[0], force[1], force[2]);

  SPDLOG_TRACE("Exiting Gravity calculate_force");

  return force;
}
}// namespace simulator::physics