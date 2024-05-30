#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <functional>
#include <ranges>

namespace container::index {

template<typename T>
concept Index =
    std::ranges::forward_range<T>
    && std::same_as<std::ranges::range_value_t<T>, std::tuple<size_t, size_t, size_t>>
    && requires(T t, std::array<double, 3> pos, std::array<size_t, 3> dim, std::array<long, 3> offset, size_t index, double width) {
         { T(pos, width) };
         { t.dimension() } -> std::convertible_to<std::array<size_t, 3>>;
         { t.position_to_index(pos) } -> std::convertible_to<size_t>;
         { t.dimension_to_index(dim) } -> std::convertible_to<size_t>;
         { t.offset(index, offset) } -> std::convertible_to<size_t>;
         { t.index_to_dimension(index) } -> std::convertible_to<std::array<size_t, 3>>;
       };

//class dim {
// public:
//  dim(std::array<double, 3> domain, double cutoff) : domain(domain) {
//  }
//  std::array<size_t, 3> dim;
//  double cell_width;
//
// private:
//  std::array<double, 3> domain;
//};

/**
 * A simple index iterating over dimensions in xyz-order
 * */
class simple_index : public std::ranges::view_interface<simple_index> {

  using View =
      std::ranges::cartesian_product_view<
          std::ranges::iota_view<size_t, size_t>,
          std::ranges::iota_view<size_t, size_t>,
          std::ranges::iota_view<size_t, size_t>>;

 private:
  std::array<double, 3> bounds;
  double width;
  View view;
  std::array<size_t, 3> dim;

 public:
  simple_index() = delete;

  explicit simple_index(std::array<double, 3> boundary, double width) : bounds(boundary), width(width) {
    auto [x, y, z] = boundary;

    dim = {(size_t) std::ceil(x / width) + 2,
           (size_t) std::ceil(y / width) + 2,
           (size_t) std::ceil(z / width) + 2};

    view = std::views::cartesian_product(
        std::views::iota(0UL, dim[0] - 2),
        std::views::iota(0UL, dim[1] - 2),
        std::views::iota(0UL, dim[2] - 2));
  }

  auto begin() -> auto {
    return view.begin();
  }
  auto end() -> auto {
    return view.end();
  }

  auto dimension() -> std::array<size_t, 3> {
    return dim;
  }

  auto boundary() -> auto {
    return bounds;
  }
  constexpr auto position_to_index(std::array<double, 3> position) -> size_t {
    auto [x, y, z] = position;

    std::array<size_t, 3> dimension = {(size_t) std::floor(x / width),
                                       (size_t) std::floor(y / width),
                                       (size_t) std::floor(z / width)};

    return dimension_to_index(dimension);
  }

  constexpr auto dimension_to_index(std::array<size_t, 3> coords) -> size_t {
    auto [x, y, z] = coords;
    return (x + 1) + ((y + 1) * dim[0]) + ((z + 1) * dim[0] * dim[1]);
  }
  constexpr auto index_to_dimension(size_t index) -> std::array<size_t, 3> {
    auto z = index / (dim[0] * dim[1]);
    auto y = (index - z * (dim[0] * dim[1])) / (dim[1]);
    auto x = (index - z * (dim[0] * dim[1]) - y * dim[0]);
    return {x - 1, y - 1, z - 1};
  }

  constexpr auto offset(size_t index, std::array<long, 3> offset) -> size_t {
    auto [x, y, z] = offset;
    return index + x + y * dim[0] + z * dim[0] * dim[1];
  }
  auto max_index() -> size_t {
    return dim[0] * dim[1] * dim[2];
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