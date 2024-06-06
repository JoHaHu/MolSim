#pragma once

#include <cstdint>
namespace container {

  /**
 * @brief Enum for defining orientations used in boundary conditions.
 *
 * Represents the six possible orientations: front, back, left, right, bottom, and top.
 */
  enum class orientation : uint8_t {
    front,
    back,
    left,
    right,
    bottom,
    top
  };
}