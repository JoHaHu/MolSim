#pragma once

#include "utils/ArrayUtils.h"
#include <array>
#include <cmath>
#include <concepts>
#include <functional>
#include <ranges>

namespace container::index {

/**
 * @brief Concept to define an Index with required methods.
 */
template<typename T>
concept Index =
    requires(T t, std::array<double, 3> pos, std::array<size_t, 3> dim, std::array<long, 3> offset, double width) {
  {T(pos, width)};
  { t.radius } -> std::convertible_to<std::array<size_t, 3>>;
  { t.dimension } -> std::convertible_to<std::array<size_t, 3>>;
  { t.width } -> std::convertible_to<std::array<double, 3>>;
  { t.boundary } -> std::convertible_to<std::array<double, 3>>;
  { t.position_to_index(pos) } -> std::convertible_to<size_t>;
  { t.dimension_to_index(dim) } -> std::convertible_to<size_t>;
  { t.offset(dim, offset) } -> std::convertible_to<size_t>;
  { t.min_distance(dim, dim) } -> std::convertible_to<double>;
};

/**
   * @brief Simple index iterating over dimensions in xyz-order.
   */
struct row_major_index {
 public:
  std::array<double, 3> boundary;
  std::array<double, 3> width;
  std::array<size_t, 3> dimension{};
  std::array<size_t, 3> radius{};

 public:
  row_major_index() = delete;

  explicit row_major_index(std::array<double, 3> boundary, double cutoff) : boundary(boundary), width({cutoff, cutoff, cutoff}) {
    auto [x, y, z] = boundary;

    dimension = {(size_t) std::ceil(x / width[0]),
                 (size_t) std::ceil(y / width[1]),
                 (size_t) std::ceil(z / width[2])};
    radius = {1, 1, 1};
  }

  /**
 * @brief Converts a position to a linear index.
 *
 * @param position The position to convert.
 * @return The linear index.
 */
  constexpr auto position_to_index(std::array<double, 3> position) -> size_t {
    auto [x, y, z] = position;

    std::array<size_t, 3> dim = {(size_t) std::floor(x / width[0]),
                                 (size_t) std::floor(y / width[1]),
                                 (size_t) std::floor(z / width[2])};

    return dimension_to_index(dim);
  }

  /**
 * @brief Converts 3D coordinates to a linear index.
 *
 * @param coords The 3D coordinates.
 * @return The linear index.
 */
  constexpr auto dimension_to_index(std::array<size_t, 3> coords) -> size_t {
    auto [x, y, z] = coords;
    return x + y * dimension[0] + z * dimension[0] * dimension[1];
  }

  /**
 * @brief Applies an offset to 3D coordinates and converts to a linear index.
 *
 * @param dim The 3D coordinates.
 * @param offset The offset to apply.
 * @return The linear index.
 */
  constexpr auto offset(std::array<size_t, 3> dim, std::array<long, 3> offset) -> size_t {
    auto off = std::array<size_t, 3>({(dim[0] + offset[0]),
                                      (dim[1] + offset[1]),
                                      (dim[2] + offset[2])});
    return dimension_to_index(off);
  }

  /**
 * @brief Calculates the minimum distance between two sets of 3D coordinates.
 *
 * @param dim1 The first set of coordinates.
 * @param dim2 The second set of coordinates.
 * @return The minimum distance.
 */
  constexpr auto min_distance(std::array<size_t, 3> dim1, std::array<size_t, 3> dim2) -> double {
    auto distances = std::vector<double>();

    for (auto [x, y, z] : std::views::cartesian_product(
             std::views::iota(0, 2),
             std::views::iota(0, 2),
             std::views::iota(0, 2))) {
      distances.emplace_back(
          distance({dim1[0] + x, dim1[1] + y, dim1[2] + z}, dim2));
    }

    double min = std::ranges::min(distances);
    return min;
  }

  /**
 * @brief Calculates the distance between two sets of 3D coordinates.
 *
 * @param dim1 The first set of coordinates.
 * @param dim2 The second set of coordinates.
 * @return The distance.
 */
  auto distance(std::array<size_t, 3> dim1, std::array<size_t, 3> dim2) -> double {
    auto diff = std::array<double, 3>({
        (double) std::abs(static_cast<long>(dim1[0]) - static_cast<long>(dim2[0])) * width[0],
        (double) std::abs(static_cast<long>(dim1[1]) - static_cast<long>(dim2[1])) * width[1],
        (double) std::abs(static_cast<long>(dim1[2]) - static_cast<long>(dim2[2])) * width[2],
    });
    return ArrayUtils::L2Norm(diff);
  }
};

static_assert(Index<row_major_index>);
/**
   * @brief Simple index iterating over dimensions in xyz-order with half-sized cells.
   */
struct half_index {
 public:
  std::array<double, 3> boundary;
  std::array<double, 3> width;
  std::array<size_t, 3> dimension;
  std::array<size_t, 3> radius;

 public:
  half_index() = delete;

  explicit half_index(std::array<double, 3> boundary, double cutoff) : boundary(boundary), width({cutoff / 2, cutoff / 2, cutoff / 2}) {
    auto [x, y, z] = boundary;

    dimension = {(size_t) std::ceil(x / width[0]),
                 (size_t) std::ceil(y / width[1]),
                 (size_t) std::ceil(z / width[2])};
    radius = {2, 2, 2};
  }

  constexpr auto position_to_index(std::array<double, 3> position) -> size_t {
    auto [x, y, z] = position;

    std::array<size_t, 3> dim = {(size_t) std::floor(x / width[0]),
                                 (size_t) std::floor(y / width[1]),
                                 (size_t) std::floor(z / width[2])};

    return dimension_to_index(dim);
  }

  constexpr auto dimension_to_index(std::array<size_t, 3> coords) -> size_t {
    auto [x, y, z] = coords;
    return x + y * dimension[0] + z * dimension[0] * dimension[1];
  }

  constexpr auto offset(std::array<size_t, 3> dim, std::array<long, 3> offset) -> size_t {
    auto off = std::array<size_t, 3>({(dim[0] + offset[0]),
                                      (dim[1] + offset[1]),
                                      (dim[2] + offset[2])});
    return dimension_to_index(off);
  }

  constexpr auto min_distance(std::array<size_t, 3> dim1, std::array<size_t, 3> dim2) -> double {
    auto distances = std::vector<double>();

    for (auto [x, y, z] : std::views::cartesian_product(
             std::views::iota(0, 2),
             std::views::iota(0, 2),
             std::views::iota(0, 2))) {
      distances.emplace_back(
          distance({dim1[0] + x, dim1[1] + y, dim1[2] + z}, dim2));
    }

    double min = std::ranges::min(distances);
    return min;
  }
  auto distance(std::array<size_t, 3> dim1, std::array<size_t, 3> dim2) -> double {
    auto diff = std::array<double, 3>({
        (double) std::abs(static_cast<long>(dim1[0]) - static_cast<long>(dim2[0])) * width[0],
        (double) std::abs(static_cast<long>(dim1[1]) - static_cast<long>(dim2[1])) * width[1],
        (double) std::abs(static_cast<long>(dim1[2]) - static_cast<long>(dim2[2])) * width[2],
    });
    return ArrayUtils::L2Norm(diff);
  }
};

static_assert(Index<half_index>);

/**
   * @brief Indexing scheme using a space-filling curve, specifically the Morton curve.
   */
class morton_index {
 public:
  std::array<double, 3> boundary;
  std::array<size_t, 3> dimension;
  std::array<double, 3> width;
  std::array<size_t, 3> radius;
  std::array<size_t, 3> log;

 public:
  morton_index() = delete;

  explicit morton_index(std::array<double, 3> boundary, double cutoff) : boundary(boundary) {
    auto [x, y, z] = boundary;
    dimension = {
        std::bit_ceil((size_t) std::ceil(boundary[0] / cutoff)),
        std::bit_ceil((size_t) std::ceil(boundary[1] / cutoff)),
        std::bit_ceil((size_t) std::ceil(boundary[2] / cutoff))};

    width = {
        x / (double) dimension[0],
        y / (double) dimension[1],
        z / (double) dimension[2],
    };

    radius = {
        static_cast<unsigned long>(std::ceil(cutoff / width[0])),
        static_cast<unsigned long>(std::ceil(cutoff / width[1])),
        static_cast<unsigned long>(std::ceil(cutoff / width[2]))};

    log = {
        static_cast<unsigned long>((std::bit_width(dimension[0] - 1))),
        static_cast<unsigned long>(std::bit_width(dimension[1] - 1)),
        static_cast<unsigned long>(std::bit_width(dimension[2] - 1)),
    };
  }

  constexpr auto position_to_index(std::array<double, 3> position) -> size_t {
    auto [x, y, z] = position;

    std::array<size_t, 3> dim = {(size_t) std::floor(x / width[0]),
                                 (size_t) std::floor(y / width[1]),
                                 (size_t) std::floor(z / width[2])};

    return dimension_to_index(dim);
  }

  constexpr auto dimension_to_index(std::array<size_t, 3> coords) -> size_t {
    auto [x, y, z] = coords;

    auto max_log = std::max({log[0], log[1], log[2]});
    size_t result = 0;
    auto i = 0;
    auto bit_position = 0;
    while (i < max_log) {
      if (i < log[0]) {
        result |= (x & (1 << i)) << (bit_position - i);
        bit_position++;
      }
      if (i < log[1]) {
        result |= (y & (1 << i)) << (bit_position - i);
        bit_position++;
      }
      if (i < log[2]) {
        result |= (z & (1 << i)) << (bit_position - i);
        bit_position++;
      }
      i++;
    }
    return result;
  }

  constexpr auto offset(std::array<size_t, 3> dim, std::array<long, 3> offset) -> size_t {
    auto off = std::array<size_t, 3>({(dim[0] + offset[0]),
                                      (dim[1] + offset[1]),
                                      (dim[2] + offset[2])});
    if (std::bit_width(off[0]) > log[0] || std::bit_width(off[1]) > log[1] || std::bit_width(off[2]) > log[2]) {
      return SIZE_MAX;
    }
    return dimension_to_index(off);
  }
  constexpr auto min_distance(std::array<size_t, 3> dim1, std::array<size_t, 3> dim2) -> double {
    auto distances = std::vector<double>();

    for (auto [x, y, z] : std::views::cartesian_product(
             std::views::iota(0, 2),
             std::views::iota(0, 2),
             std::views::iota(0, 2))) {
      distances.emplace_back(
          distance({dim1[0] + x, dim1[1] + y, dim1[2] + z}, dim2));
    }

    double min = std::ranges::min(distances);
    return min;
  }
  auto distance(std::array<size_t, 3> dim1, std::array<size_t, 3> dim2) -> double {
    auto diff = std::array<double, 3>({
        (double) std::abs(static_cast<long>(dim1[0]) - static_cast<long>(dim2[0])) * width[0],
        (double) std::abs(static_cast<long>(dim1[1]) - static_cast<long>(dim2[1])) * width[1],
        (double) std::abs(static_cast<long>(dim1[2]) - static_cast<long>(dim2[2])) * width[2],
    });
    return ArrayUtils::L2Norm(diff);
  }
};

static_assert(Index<morton_index>);

/**
   * @brief Test scheme using a space-filling curve, specifically the Morton curve.
   */
class power_index {
 public:
  std::array<double, 3> boundary;
  std::array<size_t, 3> dimension;
  std::array<double, 3> width;
  std::array<size_t, 3> radius;

 public:
  power_index() = delete;

  explicit power_index(std::array<double, 3> boundary, double cutoff) : boundary(boundary) {
    auto [x, y, z] = boundary;
    dimension = {
        std::bit_ceil((size_t) std::ceil(boundary[0] / cutoff)),
        std::bit_ceil((size_t) std::ceil(boundary[1] / cutoff)),
        std::bit_ceil((size_t) std::ceil(boundary[2] / cutoff))};

    width = {
        x / (double) dimension[0],
        y / (double) dimension[1],
        z / (double) dimension[2],
    };

    radius = {
        static_cast<unsigned long>(std::ceil(cutoff / width[0])),
        static_cast<unsigned long>(std::ceil(cutoff / width[1])),
        static_cast<unsigned long>(std::ceil(cutoff / width[2]))};
  }

  constexpr auto position_to_index(std::array<double, 3> position) -> size_t {
    auto [x, y, z] = position;

    std::array<size_t, 3> dim = {(size_t) std::floor(x / width[0]),
                                 (size_t) std::floor(y / width[1]),
                                 (size_t) std::floor(z / width[2])};

    return dimension_to_index(dim);
  }

  constexpr auto dimension_to_index(std::array<size_t, 3> coords) -> size_t {
    auto [x, y, z] = coords;
    return x + y * dimension[0] + z * dimension[0] * dimension[1];
  }

  constexpr auto offset(std::array<size_t, 3> dim, std::array<long, 3> offset) -> size_t {
    auto off = std::array<size_t, 3>({(dim[0] + offset[0]),
                                      (dim[1] + offset[1]),
                                      (dim[2] + offset[2])});
    return dimension_to_index(off);
  }

  constexpr auto min_distance(std::array<size_t, 3> dim1, std::array<size_t, 3> dim2) -> double {
    auto distances = std::vector<double>();

    for (auto [x, y, z] : std::views::cartesian_product(
             std::views::iota(0, 2),
             std::views::iota(0, 2),
             std::views::iota(0, 2))) {
      distances.emplace_back(
          distance({dim1[0] + x, dim1[1] + y, dim1[2] + z}, dim2));
    }

    double min = std::ranges::min(distances);
    return min;
  }
  auto distance(std::array<size_t, 3> dim1, std::array<size_t, 3> dim2) -> double {
    auto diff = std::array<double, 3>({
        (double) std::abs(static_cast<long>(dim1[0]) - static_cast<long>(dim2[0])) * width[0],
        (double) std::abs(static_cast<long>(dim1[1]) - static_cast<long>(dim2[1])) * width[1],
        (double) std::abs(static_cast<long>(dim1[2]) - static_cast<long>(dim2[2])) * width[2],
    });
    return ArrayUtils::L2Norm(diff);
  }
};

static_assert(Index<power_index>);

}// namespace container::index