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
  using particle_vector = std::ranges::ref_view<std::vector<std::reference_wrapper<Particle>>>;
  using product_range = std::vector<std::ranges::cartesian_product_view<
      particle_vector,
      particle_vector>>;

  using pairwise_range = ranges::concat_view<container::combination_view<particle_vector>, std::ranges::join_view<std::ranges::owning_view<product_range>>>;

  cell() = default;
  explicit cell(std::vector<std::reference_wrapper<Particle>> &&particles, cell_type type, size_t idx) : particles(std::move(particles)), type(type), idx(idx) {};

  auto linear() -> auto {
    return std::ranges::ref_view(particles);
  }
  auto pairwise() -> auto {
    return range;
  };

  constexpr auto is_boundary() -> bool {
    return type == cell_type::inner_and_halo || type == cell_type::halo;
  }

  auto insert(Particle &particle) {
    particles.emplace_back(particle);
  }

  auto clear() {
    particles.clear();
  }

 private:
  std::vector<std::reference_wrapper<Particle>> particles;
  cell_type type = cell_type::inner;
  shared_view<pairwise_range> range;
  size_t idx{};
  std::array<boundary_condition, 6> boundary = {
      boundary_condition::none,
      boundary_condition::none,
      boundary_condition::none,
      boundary_condition::none,
      boundary_condition::none,
      boundary_condition::none};

  template<index::Index I>
  static auto create_range(std::shared_ptr<linked_cell<I>> lc, std::tuple<size_t, size_t, size_t> idx) -> auto {
    auto [x, y, z] = idx;
    auto cell_idx = lc->index.dimension_to_index({x, y, z});
    auto cartesian_products = product_range();

    cell &cell = lc->cells[cell_idx];
    if (cell.type == cell_type::inner || cell.type == cell_type::inner_and_halo) {
      for (int x = -1; x <= 1; ++x) {
        for (int y = -1; y <= 1; ++y) {
          auto prod = std::views::cartesian_product(cell.linear(), lc->cells[lc->index.offset(cell_idx, {x, y, 1})].linear());
          cartesian_products.emplace_back(prod);
        }
      }
      cartesian_products.emplace_back(cell.linear(), lc->cells[lc->index.offset(cell_idx, {1, -1, 0})].linear());
      cartesian_products.emplace_back(cell.linear(), lc->cells[lc->index.offset(cell_idx, {1, 0, 0})].linear());
      cartesian_products.emplace_back(cell.linear(), lc->cells[lc->index.offset(cell_idx, {1, 1, 0})].linear());
      cartesian_products.emplace_back(cell.linear(), lc->cells[lc->index.offset(cell_idx, {0, 1, 0})].linear());
    }

    auto joined = std::move(cartesian_products) | std::views::join;
    return std::make_shared<pairwise_range>(cell.linear() | combination, std::move(joined));
  }
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
  };

  static auto make_linked_cell(const std::array<double, 3> &domain, double cutoff, boundary_condition bc, unsigned long number_of_particles) -> std::shared_ptr<linked_cell<I>> {
    auto ptr = std::make_shared<linked_cell<I>>(domain, cutoff, bc, number_of_particles);

    auto dim = ptr->index.dimension();
    for (size_t i = 0; i < dim[0] * dim[1] * dim[2]; ++i) {
      ptr->cells.emplace_back(std::vector<std::reference_wrapper<Particle>>(), cell_type::inner, i);
    }

    for (auto i : ptr->index) {
      auto [x, y, z] = i;
      auto &c = ptr->cells[ptr->index.dimension_to_index({x, y, z})];
      if (x == 0 || y == 0 || z == 0 || x == dim[0] - 1 || y == dim[1] - 1 || z == dim[2] - 1) {
        c.type = cell_type::inner_and_halo;
      }
      c.range = shared_view(cell::create_range(ptr, {x, y, z}));
    }

    return ptr;
  }

  auto boundary() -> auto {
    return cells
        | std::views::filter(&cell::is_boundary);
  }

  auto linear() -> auto {
    return arena.range();
  }

  auto pairwise() -> auto {

    return index
        | std::views::transform([this](std::tuple<size_t, size_t, size_t> idx) {
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
    std::ranges::for_each(linear(), [this](Particle &p) {
      insert_into_cell(p);
    });
    // Always remove all out of bounds particle
    // TODO fix this if we implement periodic boundary conditions
    std::ranges::for_each(arena.range_entries(), [this](container::arena<Particle>::entry &entry) {
      auto &p = entry.data;
      const auto [pos_x, pos_y, pos_z] = p.position;
      const auto [bound_x, bound_y, bound_z] = index.boundary();

      const auto condition = pos_x < 0 || pos_y < 0 || pos_z < 0 || pos_x > bound_x || pos_y > bound_y || pos_z > bound_z;
      entry.active = !condition;
    });

    spdlog::trace("fixed positions");
  }

 private:
  auto insert_into_cell(Particle &particle) {

    auto idx = index.position_to_index(particle.position);
    if (idx < index.max_index()) {
      cells[idx].insert(particle);
    } else {
      spdlog::warn("a particle has positions that is out of bounds {} {} {}", particle.position[0], particle.position[1], particle.position[2]);
    }
  }

  container::arena<Particle> arena;
  std::vector<cell> cells;
  I index;
  double cutoff;
};

}// namespace container
