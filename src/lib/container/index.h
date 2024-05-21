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

/**
 * A simple index iterating over dimensions in xyz-order
 * */
class SimpleIndex : public std::ranges::view_interface<SimpleIndex> {

  using View =
      std::ranges::cartesian_product_view<
          std::ranges::iota_view<size_t>,
          std::ranges::iota_view<size_t>,
          std::ranges::iota_view<size_t>>;

 private:
  View view;

 public:
  SimpleIndex() = delete;

  explicit SimpleIndex(std::array<size_t, 3> dim)
      : view(std::views::cartesian_product(
            std::views::iota(dim[0]),
            std::views::iota(dim[1]),
            std::views::iota(dim[2]))) {
  }

  auto begin() -> auto {
    return view.begin();
  }
  auto end() -> auto {
    return view.begin();
  }
};

static_assert(Index<SimpleIndex>);

/**
 * A indexing scheme using space filling curves, specifically the 3 dimensional Hilbert curve
 *
 * */
class HilbertIndex : public std::ranges::view_interface<HilbertIndex> {
 private:
 public:
};

//static_assert(Index<HilbertIndex>);

}// namespace container::index