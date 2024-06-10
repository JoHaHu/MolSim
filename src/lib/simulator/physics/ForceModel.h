#pragma once

#include "Gravity.h"
#include "LennardJones.h"
#include <cstdint>
#include <variant>
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
concept physics = requires(T type, Particle const &particle1, Particle const &particle2) {
  { type.calculate_force(particle1, particle2) } -> std::convertible_to<std::array<double, 3>>;
};

template<physics... p>
using force_model_variant = std::variant<p...>;

using force_model = force_model_variant<LennardJones, Gravity>;

static constexpr auto calculate_force(force_model &force, Particle &particle1, Particle &particle2) -> std::array<double, 3> {
  return std::visit(
      [&particle1, &particle2](auto &&force) { return force.calculate_force(particle1, particle2); }, force);
}

/**
 * A enum with supported Force Models
 * */
enum class ForceModel : uint8_t {
  Gravity,
  LennardJones
};
}// namespace simulator::physics
