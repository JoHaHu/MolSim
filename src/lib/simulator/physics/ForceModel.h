#pragma once

#include "Gravity.h"
#include "LennardJones.h"
#include <cstdint>
namespace simulator::physics {

/** <p> physics concept for force calculation </p>
* calculates the force for all particles,
*
* This is a concept since it allows the physics to be transparent for compiler optimizations.
* Trade off is binary size.
* \param T
* \param particle1
* \param particle2
*/
template<typename T>
concept physics = requires(T type, const double &value, const int &p_type) {
  { type.calculate_force(value, value, value, value, p_type, value, value, value, value, p_type) } -> std::convertible_to<std::tuple<double, double, double>>;
};

template<physics... p>
using force_model_variant = std::variant<p...>;

using force_model = force_model_variant<LennardJones, Gravity>;

static constexpr auto calculate_force(force_model &force, const double &position_x, const double &position_y, const double &position_z, const double &mass, const int &type,
                                      const double &position_x2, const double &position_y2, const double &position_z2, const double &mass2, const int &type2) -> std::tuple<double, double, double> {
  return std::visit(
      [&](auto &&force) { return force.calculate_force(position_x, position_y, position_z, mass, type,
                                                       position_x2, position_y2, position_z2, mass2, type2); }, force);
}

/**
 * A enum with supported Force Models
 * */
enum class ForceModel : uint8_t {
  Gravity,
  LennardJones
};
}// namespace simulator::physics
