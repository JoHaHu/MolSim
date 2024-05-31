#include "LennardJones.h"

#include <array>
#include <cmath>
#include <limits>

#include "spdlog/spdlog.h"
#include "utils/ArrayUtils.h"

namespace simulator::physics {
/**
       * @brief Calculates the Lennard Jones force between two particles.
       *
       * Computes the Lennard Jones forces based on the positions of two particles.
       *
       * @param particle1 The first particle.
       * @param particle2 The second particle.
       * @return std::array<double, 3> The calculated force vector.
       */
auto LennardJones::calculate_force(const Particle &particle1, const Particle &particle2) -> std::array<double, 3> {
  const auto epsilon = 5;
  const auto sigma = 1;

  SPDLOG_TRACE("Entering LennardJones calculate_force");

  SPDLOG_TRACE("Particle 1: mass = {}, position = ({}, {}, {})", particle1.mass, particle1.position[0], particle1.position[1], particle1.position[2]);
  SPDLOG_TRACE("Particle 2: mass = {}, position = ({}, {}, {})", particle2.mass, particle2.position[0], particle2.position[1], particle2.position[2]);

  const auto x_diff = particle2.position - particle1.position;

  SPDLOG_TRACE("Position difference: ({}, {}, {})", x_diff[0], x_diff[1], x_diff[2]);

  const auto norm = ArrayUtils::L2Norm(x_diff);
  SPDLOG_TRACE("Norm of position difference: {}", norm);

  if (norm == 0) {
    SPDLOG_WARN("Zero distance between particles encountered. This should not be possible.");
  }

  const auto sigma_over_norm_3 = (sigma / norm) * (sigma / norm) * (sigma / norm);
  const auto sigma_over_norm_6 = sigma_over_norm_3 * sigma_over_norm_3;
  const auto sigma_over_norm_12 = sigma_over_norm_6 * sigma_over_norm_6;

  const auto force = 24 * epsilon / (norm * norm) * (sigma_over_norm_6 - 2 * sigma_over_norm_12) * x_diff;

  SPDLOG_TRACE("Calculated force: ({}, {}, {})", force[0], force[1], force[2]);

  SPDLOG_TRACE("Exiting LennardJones calculate_force");

  return force;
}
}// namespace simulator::physics
