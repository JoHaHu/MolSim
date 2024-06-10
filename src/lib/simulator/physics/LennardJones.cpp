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
auto LennardJones::calculate_force(const double &position_x, const double &position_y, const double &position_z, const double &mass, const int &type, const double &position_x2, const double &position_y2, const double &position_z2, const double &mass2, const int &type2) -> std::tuple<double, double, double> {

  SPDLOG_TRACE("Entering LennardJones calculate_force");

  const auto x_diff = position_x2 - position_x;
  const auto y_diff = position_y2 - position_y;
  const auto z_diff = position_z2 - position_z;

  const auto diff = {x_diff, y_diff, z_diff};

  const auto norm = ArrayUtils::L2Norm(diff);
  SPDLOG_TRACE("Norm of position difference: {}", norm);

  if (norm == 0) {
    SPDLOG_WARN("Zero distance between particles encountered. This should not be possible.");
  }

  if (norm > cutoff) {
    return {0.0, 0.0, 0.0};
  }

  const auto sigma_over_norm_3 = (sigma / norm) * (sigma / norm) * (sigma / norm);
  const auto sigma_over_norm_6 = sigma_over_norm_3 * sigma_over_norm_3;

  const auto temp = 24 * epsilon / (norm * norm) * (1 - 2 * sigma_over_norm_6) * sigma_over_norm_6;
  const auto force_x = temp * x_diff;
  const auto force_y = temp * y_diff;
  const auto force_z = temp * z_diff;
  SPDLOG_TRACE("Exiting LennardJones calculate_force");

  return {force_x, force_y, force_z};
}

LennardJones::LennardJones(double cutoff, double sigma, double epsilon) : epsilon(epsilon), sigma(sigma), cutoff(cutoff) {
}

}// namespace simulator::physics
