#pragma once

#include <cstdint>
namespace simulator::physics {
enum class ForceModel : uint8_t {
  Gravity,
  LennardJones
};
}
