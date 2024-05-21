#pragma once

#include "Particle.h"
#include "combination.h"
#include "container/index.h"
#include <array>
#include <concepts>
#include <list>
#include <memory>
#include <ranges>
#include <vector>

#include "range/v3/view/concat.hpp"

namespace container {

class Cell {
 public:
  Cell() = default;
  explicit Cell(std::vector<Particle> &&particles) : particles(std::move(particles)){};

  auto begin() -> std::vector<Particle>::iterator {
    return particles.begin();
  }

  auto end() -> std::vector<Particle>::iterator {
    return particles.end();
  }

 private:
  std::vector<Particle> particles;
};

static_assert(std::ranges::forward_range<Cell>);

class LinkedCellFunctor {
 private:
  // SAFETY: the LinkedCells this reference comes from must always outlive this functor
  std::vector<Cell> &cells;
  std::array<size_t, 3> dim;

 public:
  LinkedCellFunctor() = delete;
  explicit LinkedCellFunctor(std::vector<Cell> &cells, std::array<size_t, 3> dim) : cells(cells), dim(dim) {}

  constexpr auto to_index(size_t x, size_t y, size_t z) -> size_t {
    return dim[0] * x + dim[1] * y + dim[2] * z;
  }

  auto operator()(std::tuple<size_t, size_t, size_t> idx) -> auto {
    const auto [x, y, z] = idx;

    const auto cell_index = to_index(x, y, z);
    auto &cell = cells[cell_index];

    auto combinations = std::ranges::ref_view(cell) | combination;
    // TODO check bounds
    // All cells with a coordinate that is +1 bigger. (Newton 3rd law)
    auto prods = std::vector{
        std::views::cartesian_product(cells[to_index(x + 1, y - 1, z)], cell),
        std::views::cartesian_product(cells[to_index(x + 1, y, z)], cell),
        std::views::cartesian_product(cells[to_index(x + 1, y + 1, z)], cell),
        std::views::cartesian_product(cells[to_index(x, y + 1, z)], cell),
    };

    for (int dx : std::views::iota(-1, 1)) {
      for (int dy : std::views::iota(-1, 1)) {
        prods.emplace_back(std::views::cartesian_product(cells[to_index(x + dx, y + dy, z + 1)], cell));
      }
    }

    auto outer = prods | std::views::join;

    return ranges::views::concat(combinations, outer);
  }
};

template<index::Index I>
class LinkedCells : public std::ranges::view_interface<LinkedCells<I>> {

  // unfortunately cant be auto-inferred
  using view = std::ranges::join_view<std::ranges::transform_view<I, LinkedCellFunctor>>;

 public:
  LinkedCells() = default;
  explicit LinkedCells(const std::array<size_t, 3> &dim) : dim(dim) {
    cells.reserve(dim[0] * dim[1] * dim[2]);
  };

  auto begin() -> std::ranges::iterator_t<view> {
    return iter.begin();
  }
  auto end() -> std::ranges::sentinel_t<view> {
    return iter.end();
  }

 private:
  std::array<size_t, 3> dim{};
  std::vector<Cell> cells;
  I index = I(dim);
  LinkedCellFunctor cell_functor = LinkedCellFunctor(cells, dim);

  view iter = index
      | std::views::transform(cell_functor)
      | std::views::join;
};

//static_assert(std::ranges::forward_range<LinkedCells<index::SimpleIndex>>);
}// namespace container
