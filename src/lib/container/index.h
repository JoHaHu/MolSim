#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <functional>
#include <ranges>

namespace container::index {

template<typename T>
concept Index =
    requires(T t, std::array<double, 3> pos, std::array<size_t, 3> dim, std::array<long, 3> offset, size_t index, double width) {
      { T(pos, width) };
      { t.radius } -> std::convertible_to<std::array<size_t, 3>>;
      { t.dimension } -> std::convertible_to<std::array<size_t, 3>>;
      { t.width } -> std::convertible_to<std::array<double, 3>>;
      { t.boundary } -> std::convertible_to<std::array<double, 3>>;
      { t.position_to_index(pos) } -> std::convertible_to<size_t>;
      { t.dimension_to_index(dim) } -> std::convertible_to<size_t>;
      { t.offset(index, offset) } -> std::convertible_to<size_t>;
    };

/**
 * A simple index iterating over dimensions in xyz-order
 * */
struct simple_index : public std::ranges::view_interface<simple_index> {

  using View =
      std::ranges::cartesian_product_view<
          std::ranges::iota_view<size_t, size_t>,
          std::ranges::iota_view<size_t, size_t>,
          std::ranges::iota_view<size_t, size_t>>;

 public:
  std::array<double, 3> boundary;
  std::array<double, 3> width;
  std::array<size_t, 3> dimension;
  std::array<size_t, 3> radius;

 public:
  simple_index() = delete;

  explicit simple_index(std::array<double, 3> boundary, double cutoff) : boundary(boundary), width({cutoff, cutoff, cutoff}) {
    auto [x, y, z] = boundary;

    dimension = {(size_t) std::ceil(x / width[0]),
                 (size_t) std::ceil(y / width[1]),
                 (size_t) std::ceil(z / width[2])};
    radius = {
        (size_t) std::ceil(cutoff / width[0]),
        (size_t) std::ceil(cutoff / width[1]),
        (size_t) std::ceil(cutoff / width[2])};
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

  constexpr auto offset(size_t index, std::array<long, 3> offset) -> size_t {
    auto [x, y, z] = offset;
    return index + x + y * dimension[0] + z * dimension[0] * dimension[1];
  }
};

static_assert(Index<simple_index>);

///**
// * A indexing scheme using space filling curves, specifically a compact version of the 3 dimensional Hilbert curve with different order each side
// * */
//class hilbert_index : public std::ranges::view_interface<hilbert_index> {
// private:
// public:
//};

//static_assert(Index<HilbertIndex>);

}// namespace container::index