#pragma once

#include "Particle.h"
#include "arena.h"
#include "combination.h"
#include "container.h"
#include "index.h"
#include "utils/ArrayUtils.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <list>
#include <memory>
#include <ranges>
#include <utility>
#include <vector>

#include "range/v3/view/concat.hpp"
#include "shared.h"

namespace container {

enum class boundary_condition : std::uint8_t {
  outflow,
  reflecting,
  //  periodic,
  none
};

enum class cell_type : std::uint8_t {
  inner,
  inner_and_halo,
  halo
};

template<index::Index I>
class linked_cell;

class cell {
  template<index::Index I>
  friend class linked_cell;

 public:
  using particle_vector = std::ranges::transform_view<std::ranges::ref_view<std::vector<std::reference_wrapper<arena<Particle>::entry>>>, decltype(&arena<Particle>::entry::data)>;
  using product_range = std::vector<std::ranges::cartesian_product_view<
      particle_vector,
      particle_vector>>;

  using pairwise_range = ranges::concat_view<container::combination_view<particle_vector>, std::ranges::join_view<std::ranges::owning_view<product_range>>>;

  cell() = default;
  explicit cell(std::vector<std::reference_wrapper<arena<Particle>::entry>> &&particles, cell_type type, size_t idx) : particles(std::move(particles)), type(type), idx(idx) {};

  auto linear() -> auto {
    return particles
        | std::views::transform(&container::arena<Particle>::entry::data);
  }
  auto pairwise() -> auto & {
    return *range;
  };

  constexpr auto is_boundary() -> bool {
    return type == cell_type::inner_and_halo || type == cell_type::halo;
  }

  auto insert(arena<Particle>::entry &particle) {
    particles.emplace_back(particle);
  }

  auto clear() {
    particles.clear();
  }

 private:
  std::vector<std::reference_wrapper<arena<Particle>::entry>> particles;
  cell_type type = cell_type::inner;
  std::optional<pairwise_range> range;
  size_t idx{};
  std::array<boundary_condition, 6> boundary = {
      boundary_condition::none,
      boundary_condition::none,
      boundary_condition::none,
      boundary_condition::none,
      boundary_condition::none,
      boundary_condition::none};
};

template<index::Index I>
class linked_cell {
  friend class cell;

 public:
  linked_cell() = delete;
  explicit linked_cell(const std::array<double, 3> &domain, double cutoff, boundary_condition bc, unsigned int number_of_particles)
      : arena(number_of_particles),
        index(I(domain, cutoff)),
        cutoff(cutoff) {
    auto dim = index.dimension();
    cells.reserve(dim[0] * dim[1] * dim[2]);

    for (size_t i = 0; i < dim[0] * dim[1] * dim[2]; ++i) {
      this->cells.emplace_back(std::vector<std::reference_wrapper<container::arena<Particle>::entry>>(), cell_type::inner, i);
    }

    for (auto i : this->index) {
      auto [x, y, z] = i;
      auto &c = this->cells[this->index.dimension_to_index({x, y, z})];
      if (x == 0 || y == 0 || z == 0 || x == dim[0] - 1 || y == dim[1] - 1 || z == dim[2] - 1) {
        c.type = cell_type::inner_and_halo;
      }
      c.range = create_range(*this, {x, y, z});
    }
  };

  static auto create_range(linked_cell &lc, std::tuple<size_t, size_t, size_t> idx) -> auto {
    auto [x, y, z] = idx;
    auto cell_idx = lc.index.dimension_to_index({x, y, z});

    cell &cell = lc.cells[cell_idx];
    if (cell.type == cell_type::inner || cell.type == cell_type::inner_and_halo) {

      auto cartesian_products = cell::product_range();
      cartesian_products.reserve(13);
      for (int x = -1; x <= 1; ++x) {
        for (int y = -1; y <= 1; ++y) {
          cartesian_products.emplace_back(cell.linear(), lc.cells[lc.index.offset(cell_idx, {x, y, 1})].linear());
        }
      }
      cartesian_products.emplace_back(cell.linear(), lc.cells[lc.index.offset(cell_idx, {1, -1, 0})].linear());
      cartesian_products.emplace_back(cell.linear(), lc.cells[lc.index.offset(cell_idx, {1, 0, 0})].linear());
      cartesian_products.emplace_back(cell.linear(), lc.cells[lc.index.offset(cell_idx, {1, 1, 0})].linear());
      cartesian_products.emplace_back(cell.linear(), lc.cells[lc.index.offset(cell_idx, {0, 1, 0})].linear());

      auto joined = std::move(cartesian_products) | std::views::join;
      return std::optional(cell::pairwise_range(cell.linear() | combination, std::move(joined)));
    }

    return std::optional<cell::pairwise_range>();
  }

  auto boundary() -> auto {
    return cells
        | std::views::filter(&cell::is_boundary);
  }

  auto linear() -> auto {
    return arena.range_entries()
        | std::views::transform(&container::arena<Particle>::entry::data);
  }

  auto pairwise() -> auto {

    return index
        | std::views::transform([this](std::tuple<size_t, size_t, size_t> idx) -> auto & {
             auto [x, y, z] = idx;
             auto cell_idx = index.dimension_to_index({x, y, z});
             cell &cell = cells[cell_idx];
             return cell.pairwise();
           })
        | std::views::join;
  }

  auto insert(Particle &particle) {
    auto &inserted = arena.emplace_back(particle);
    insert_into_cell(inserted);
  }

  auto size() {
    return arena.size();
  }

  auto fix_positions() {
    std::ranges::for_each(cells, &cell::clear);
    std::ranges::for_each(arena.range_entries(), [this](container::arena<Particle>::entry &p) {
      insert_into_cell(p);
    });
    spdlog::trace("fixed positions");
  }

 private:
  auto insert_into_cell(arena<Particle>::entry &particle) {

    auto idx = index.position_to_index(particle.data.position);
    if (idx < index.max_index()) {
      cells[idx].insert(particle);
    } else {
      spdlog::warn("a particle has positions that is out of bounds {} {} {}", particle.data.position[0], particle.data.position[1], particle.data.position[2]);
    }
  }

  container::arena<Particle> arena;
  std::vector<cell> cells;
  I index;
  double cutoff;
};

}// namespace container
