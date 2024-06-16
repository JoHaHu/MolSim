#pragma once

#include<cstdint>
/**
 * @brief Enum for defining different types of boundary conditions.
 */
enum class BoundaryCondition : uint8_t {
  none,
  outflow,
  reflecting,
  periodic
};
