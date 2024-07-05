#pragma once

#include "boundary.h"
#include "utils/ArrayUtils.h"
#include <climits>
#include <cmath>
#include <span>

namespace container::index {

constexpr const std::array<size_t, 3> block_size = {2, 2, 1};

template<const size_t DIMENSIONS>
class Index {

 public:
  std::array<double, DIMENSIONS> domain{};
  std::array<BoundaryCondition, 2 * DIMENSIONS> bc{};
  std::array<size_t, DIMENSIONS> dim{};
  std::array<double, DIMENSIONS> widths{};
  std::array<long, DIMENSIONS> radius{};

  Index() = default;

  Index(const std::array<double, DIMENSIONS> &dm, const std::array<BoundaryCondition, 2 * DIMENSIONS> &bc, double cutoff)
      : domain(dm),
        bc(bc),
        dim([this, cutoff] {
          std::array<size_t, DIMENSIONS> temp = {};
          for (size_t i = 0; i < DIMENSIONS; ++i) {
            temp[i] = ((size_t) std::floor(domain[i] / cutoff));
          }
          return temp;
        }()),
        widths([this] {
          std::array<double, DIMENSIONS> temp = {};
          for (size_t i = 0; i < DIMENSIONS; ++i) {
            temp[i] = domain[i] / dim[i];
          }
          return temp;
        }()),
        radius([this, cutoff] {
          std::array<long, DIMENSIONS> temp = {};
          for (size_t i = 0; i < DIMENSIONS; ++i) {
            temp[i] = (size_t) std::ceil(cutoff / widths[i]);
          }
          return temp;
        }()) {
  }

  /**
   * The color of the block in the domain. Non overlapping partitioning of the domain in 4 colors to parallelize lockless
   * */
  auto color(std::array<size_t, DIMENSIONS> idx) -> size_t {
    auto color = 0;

    for (size_t d = DIMENSIONS - 1; d < DIMENSIONS; --d) {
      // width and modulo together defines the sidelength in one dimension of the block
      auto width = DIMENSIONS - d;
      auto modulo = width * (width + 1);

      //the shift distance
      auto shift = 2 + width;
      if (d == 0) {
        shift = 1;
      }
      /*
       * calculates the position in one dimension applies the shift from the previous dimensions
       * modulo is then applied to to restrict values to the range necessary for one period of the pattern
       * The division by the width is then used to determine which color it is
       * */
      color = shift * ((((idx[d]) / (radius[d] * block_size[d]) + color) % modulo) / width);
    }
    return color;
  }
  /**
   * the partition id of a block in the domain. A block contains 1x2x3 cells (in 3D) or 1x2 cells (in 2D)
   * */
  auto block_id(std::array<size_t, DIMENSIONS> idx) -> size_t {
    auto block_id = 0;
    auto dim_shift = 0;

    auto dim_size = 1;

    for (size_t d = 0; d < DIMENSIONS; ++d) {
      dim_size *= (size_t) std::ceil((double) dim[d] / (radius[d] * block_size[d]) / (DIMENSIONS - d));
    }

    for (size_t d = DIMENSIONS - 1; d < DIMENSIONS; --d) {

      dim_size /= (dim[d] / (radius[d] * block_size[d]) / (DIMENSIONS - d));

      auto width = DIMENSIONS - d;
      auto modulo = width * (width + 1);

      //the shift distance
      auto shift = 2 + width;

      auto period = ((idx[d]) / (radius[d] * block_size[d]) + dim_shift) / width;
      block_id += period * dim_size;
      dim_shift = shift * (((idx[d]) / (radius[d] * block_size[d]) + dim_shift) % modulo);
    }
    return block_id;
  }

  /**
   * For now this accepts all cell in the radius, but this can be improved with a proper calculation
   * */
  auto in_cutoff_distance(std::array<size_t, DIMENSIONS> dim1, std::array<size_t, DIMENSIONS> dim2) -> bool {
    // TODO
    return true;
  }

  auto distance(std::array<double, DIMENSIONS> dim1, std::array<double, DIMENSIONS> dim2) -> double {
    auto diff = dim1 - dim2;
    return ArrayUtils::L2Norm(diff);
  }

  auto calculate_correction(std::array<size_t, DIMENSIONS> from, std::array<long, DIMENSIONS> offset) -> std::array<double, DIMENSIONS> {
    std::array<double, DIMENSIONS> correction;
    for (size_t i = 0; i < DIMENSIONS; ++i) {
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
    for (size_t i = 0; i < DIMENSIONS; ++i) {
      dimension[i] = (size_t) std::floor(position[i] / widths[i]);
    }
    return dimension_to_index(dimension);
  }

  auto dimension_to_index(std::array<size_t, DIMENSIONS> dimension) -> size_t {
    size_t result = 0;
    for (size_t i = 0; i < DIMENSIONS; ++i) {
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
    for (size_t i = 0; i < DIMENSIONS; ++i) {
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