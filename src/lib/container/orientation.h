#pragma once

#include <cstdint>
namespace container::boundary {

enum class orientation : uint8_t {
  front,
  back,
  left,
  right,
  bottom,
  top
};
}
