#pragma once

#include "boundary.h"
#include "utils/ArrayUtils.h"
#include "utils/types.h"
#include <climits>
#include <cmath>
#include <span>

namespace container::index {

template<const size_t DIMENSIONS>
class Index {

 public:
  std::array<double, DIMENSIONS> domain{};
  std::array<BoundaryCondition, 2 * DIMENSIONS> bc{};
  std::array<double, DIMENSIONS> widths{};
  std::array<size_t, DIMENSIONS> dim{};
  std::array<long, DIMENSIONS> radius{};
  //  double diagonal;

  Index() = default;

  Index(const std::array<double, DIMENSIONS> &dm, const std::array<BoundaryCondition, 2 * DIMENSIONS> &bc, double cutoff) : domain(dm), bc(bc) {

    radius = [] {
      std::array<long, DIMENSIONS> temp = {};
      for (int i = 0; i < DIMENSIONS; ++i) {
        temp[i] = 1;
      }
      return temp;
    }();

    widths = [cutoff, this] {
      std::array<double, DIMENSIONS> temp = {};
      for (int i = 0; i < DIMENSIONS; ++i) {
        temp[i] = cutoff / radius[i];
      }
      return temp;
    }();

    dim = [this] {
      std::array<size_t, DIMENSIONS> temp = {};
      for (int i = 0; i < DIMENSIONS; ++i) {
        temp[i] = (size_t) std::ceil(domain[i] / widths[i]);
      }
      return temp;
    }();

    //    diagonal = distance(, widths);
  }

  auto in_cutoff_distance(std::array<size_t, DIMENSIONS> dim1, std::array<size_t, DIMENSIONS> dim2) -> bool {
    // TODO fix this
    return true;
  }

  auto distance(std::array<double, DIMENSIONS> dim1, std::array<double, DIMENSIONS> dim2) -> double {
    auto diff = dim1 - dim2;
    return ArrayUtils::L2Norm(diff);
  }

  auto calculate_correction(std::array<size_t, DIMENSIONS> from, std::array<long, DIMENSIONS> offset) -> std::array<double, DIMENSIONS> {
    std::array<double, DIMENSIONS> correction;
    for (int i = 0; i < DIMENSIONS; ++i) {
      correction[i] = 0.0;
      if (bc[i] == BoundaryCondition::periodic) {
        long diff = static_cast<long>(from[i]) + offset[i];
        if (diff >= static_cast<long>(dim[i])) {
          correction[i] = domain[i];
        } else if (diff < 0) {
          correction[i] = -domain[i];
        }
      }
    }
    return correction;
  }

  constexpr auto position_to_index(std::array<double, DIMENSIONS> position) -> size_t {
    std::array<size_t, DIMENSIONS> dimension;
    for (int i = 0; i < DIMENSIONS; ++i) {
      dimension[i] = (size_t) std::floor(position[i] / widths[i]);
    }
    return dimension_to_index(dimension);
  }

  auto dimension_to_index(std::array<size_t, DIMENSIONS> dimension) -> size_t {
    size_t result = 0;
    for (int i = 0; i < DIMENSIONS; ++i) {
      if (bc[i] == BoundaryCondition::periodic) {
        dimension[i] = dimension[i] - (dimension[i] / dim[i]) * dim[i];
      }
      if (dimension[i] >= dim[i]) {
        return ULONG_MAX;
      }
      result += dimension[i] * std::accumulate(dim.begin(), std::next(dim.begin(), i), 1, std::multiplies<>());
    }
    return result;
  }

  auto offset(std::array<size_t, DIMENSIONS> index, std::array<long, DIMENSIONS> off) -> size_t {

    std::array<size_t, DIMENSIONS> temp;
    for (int i = 0; i < DIMENSIONS; ++i) {
      temp[i] = 0;
      if (bc[i] == BoundaryCondition::periodic) {
        temp[i] += dim[i];
      }
      temp[i] += static_cast<long>(index[i]) + off[i];
    }
    return dimension_to_index(temp);
  }
};

}// namespace container::index