#pragma once

#include "Gravity.h"
#include "LennardJones.h"
#include <cstdint>
namespace simulator::physics {



/**
 * A enum with supported Force Models
 * */
enum class ForceModel : uint8_t {
  Gravity,
  LennardJones
};
}// namespace simulator::physics
