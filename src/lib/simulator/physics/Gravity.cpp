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
auto Gravity::calculate_force(
    const double &position_x, const double &position_y, const double &position_z, const double &mass, const int &type,
    const double &position_x2, const double &position_y2, const double &position_z2, const double &mass2, const int &type2

    ) -> std::tuple<double, double, double> {
  SPDLOG_TRACE("Entering Gravity calculate_force");

  const auto x_diff = position_x2 - position_x;
  const auto y_diff = position_y2 - position_y;
  const auto z_diff = position_z2 - position_z;

  const auto diff = {x_diff, y_diff, z_diff};

  const auto norm = ArrayUtils::L2Norm(diff);
  SPDLOG_TRACE("Norm of position difference: {}", norm);

  if (norm == 0) {
    SPDLOG_DEBUG("Zero distance between particles encountered");
  }

  const auto temp = (mass * mass2) / (norm * norm * norm);
  const auto force_x = temp * x_diff;
  const auto force_y = temp * y_diff;
  const auto force_z = temp * z_diff;
  SPDLOG_TRACE("Exiting Gravity calculate_force");

  return {force_x, force_y, force_z};
}
}// namespace simulator::physics