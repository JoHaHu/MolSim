#pragma once

#include <array>
#include <concepts>
#include <functional>
#include <ranges>

namespace container::index {

template<typename T>
concept Index =
    std::ranges::forward_range<T>
    && std::same_as<std::ranges::range_value_t<T>, std::tuple<size_t, size_t, size_t>>
    && requires(T, std::array<size_t, 3> dim) {
         { T(dim) };
       };

//class dimension {
// public:
//  dimension(std::array<double, 3> domain, double cutoff) : domain(domain) {
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
  View view;

 public:
  simple_index() = delete;

  explicit simple_index(std::array<size_t, 3> dim)
      : view(std::views::cartesian_product(
            std::views::iota(0UL, dim[0]),
            std::views::iota(0UL, dim[1]),
            std::views::iota(0UL, dim[2]))) {
  }

  auto begin() -> auto {
    return view.begin();
  }
  auto end() -> auto {
    return view.end();
  }
};

static_assert(Index<simple_index>);

/**
 * A indexing scheme using space filling curves, specifically a compact version of the 3 dimensional Hilbert curve with different order each side
 * */
class hilbert_index : public std::ranges::view_interface<hilbert_index> {
 private:
 public:
};

//static_assert(Index<HilbertIndex>);

}// namespace container::index