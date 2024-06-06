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

#include "orientation.h"
#include "range/v3/view/concat.hpp"
#include "spdlog/spdlog.h"

namespace container {

enum class boundary_condition : std::uint8_t {
  outflow,
  reflecting,
  //  periodic,
  none
};

enum class cell_type : std::uint8_t {
  inner,
  boundary,
};

template<index::Index I>
class linked_cell;

/**
 * A cell in the linked cell datastructure
 * */
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
  explicit cell(std::vector<std::reference_wrapper<arena<Particle>::entry>> &&particles,
                cell_type type,
                std::array<size_t, 3> idx,
                std::array<double, 3> widths) : particles(std::move(particles)), type(type), idx(idx), widths(widths) {};
  /**
   * a linear range over the particles
   * */
  auto linear() -> auto {
    return particles
        | std::views::transform(&container::arena<Particle>::entry::data);
  }

  /**
   * a pairwise range over the particles, uses a cached range initialized at the start of the program
   * */
  auto pairwise() -> auto & {
    return *range;
  };

  constexpr auto is_boundary() -> bool {
    return type == cell_type::boundary;
  }

  /**
   * inserts a new particle into the cell
   * */
  auto insert(arena<Particle>::entry &particle) {
    particles.emplace_back(particle);
  }

  /**
   * clears the cell
   * */
  auto clear() {
    particles.clear();
  }

 public:
  std::vector<std::reference_wrapper<arena<Particle>::entry>> particles;
  cell_type type = cell_type::inner;
  std::optional<pairwise_range> range;
  std::array<size_t, 3> idx{};
  std::array<double, 3> widths{};
  std::array<boundary_condition, 6> boundary = {
      boundary_condition::none,
      boundary_condition::none,
      boundary_condition::none,
      boundary_condition::none,
      boundary_condition::none,
      boundary_condition::none};
};

/**
 * The linked cell container class
 * the Index template parameter can be used to change the iteration order of the LinkedCell.
 * */
template<index::Index I>
class linked_cell {
  friend class cell;

 public:
  linked_cell() = delete;
  explicit linked_cell(const std::array<double, 3> &domain, double cutoff, std::array<boundary_condition, 6> bc, unsigned int number_of_particles, double sigma)
      : arena(number_of_particles),
        bc(bc),
        index(I(domain, cutoff)),
        sigma(sigma),
        cutoff(cutoff) {
    auto dim = index.dimension;
    auto widths = index.width;
    cells.reserve(dim[0] * dim[1] * dim[2]);

    auto range = std::views::cartesian_product(
        std::views::iota(0UL, dim[0]),
        std::views::iota(0UL, dim[1]),
        std::views::iota(0UL, dim[2]));

    for (auto [x, y, z] : range) {
      this->cells.emplace_back(std::vector<std::reference_wrapper<container::arena<Particle>::entry>>(), cell_type::inner, std::array<size_t, 3>({x, y, z}), widths);
    }

    // sort cells by order defined by the index
    std::ranges::sort(cells, [this](auto &cell1, auto &cell2) -> bool {
      return index.dimension_to_index(cell1.idx) < index.dimension_to_index(cell2.idx);
    });

    for (auto &c : cells) {
      auto [x, y, z] = c.idx;
      if (x == 0) {
        c.type = cell_type::boundary;
        c.boundary[(size_t) boundary::orientation::left] = bc[(size_t) boundary::orientation::left];
      }
      if (x == dim[0] - 1) {
        c.type = cell_type::boundary;
        c.boundary[(size_t) boundary::orientation::right] = bc[(size_t) boundary::orientation::right];
      }
      if (y == 0) {
        c.type = cell_type::boundary;
        c.boundary[(size_t) boundary::orientation::bottom] = bc[(size_t) boundary::orientation::bottom];
      }
      if (y == dim[1] - 1) {
        c.type = cell_type::boundary;
        c.boundary[(size_t) boundary::orientation::top] = bc[(size_t) boundary::orientation::top];
      }
      if (z == 0) {
        c.type = cell_type::boundary;
        c.boundary[(size_t) boundary::orientation::back] = bc[(size_t) boundary::orientation::back];
      }
      if (z == dim[2] - 1) {
        c.type = cell_type::boundary;
        c.boundary[(size_t) boundary::orientation::front] = bc[(size_t) boundary::orientation::front];
      }

      c.range = create_range(*this, {x, y, z});
    }
  };

  /**
   * Creates the range for a cell, containing the neighbours needed to calculate forces and respects Newtons 3. Law
   * Can be calculated once at the start and then be reused
   * */
  static auto create_range(linked_cell &lc, std::array<size_t, 3> idx) -> auto {
    auto cell_idx = lc.index.dimension_to_index(idx);

    cell &cell = lc.cells[cell_idx];

    auto [radius_x, radius_y, radius_z] = lc.index.radius;

    auto cartesian_products = cell::product_range();

    for (long x = -(long) radius_x; x <= (long) radius_x; ++x) {
      for (long y = -(long) radius_y; y <= (long) radius_y; ++y) {
        for (long z = -(long) radius_z; z <= (long) radius_z; ++z) {
          if (x != 0 || y != 0 || z != 0) {
            if (z > 0 || (x > 0 && z >= 0) || (x >= 0 && z >= 0 && y > 0)) {
              auto other_index = lc.index.offset(idx, {x, y, z});
              // checks that cell is not out of bounds and not the same cell
              // TODO calculate min distance between the two cells to skip if not reachable within cutoff distance
              if (other_index < lc.cells.size() && other_index != cell_idx && lc.index.min_distance(idx, lc.cells[other_index].idx) <= lc.cutoff) {
                cartesian_products.emplace_back(cell.linear(), lc.cells[other_index].linear());
              }
            }
          }
        }
      }
    }
    cartesian_products.shrink_to_fit();
    auto joined = std::move(cartesian_products) | std::views::join;
    return std::optional(cell::pairwise_range(cell.linear() | combination, std::move(joined)));
  }

  /**
   * a range over the boundary cells
   * */
  auto boundary() -> auto {
    return cells
        | std::views::filter(&cell::is_boundary);
  }

  /**
   * returns a linear range of Particles, stored in arena entries
   * */
  auto linear() -> auto {
    return arena.range_entries()
        | std::views::transform(&container::arena<Particle>::entry::data);
  }

  /**
   * returns a range with pairs of Particles
   * */
  auto pairwise() -> auto {

    return cells
        | std::views::transform(&cell::pairwise)
        | std::views::join;
  }

  /**
   * inserts a new particle into the arena and the into the cells
   * */
  auto insert(Particle &particle) {
    auto &inserted = arena.emplace_back(particle);
    insert_into_cell(inserted);
  }

  /**
   * Returns the number of particles in the arena. O(n), because it gets counted each time
   * */
  auto size() {
    return arena.size();
  }

  /**
   * Clear all cells and reinserts all particles into the cells.
   * Better than updating positions when changes are made, because the vectors get resized/reorderd a lot when only single elements get removed
   * */
  auto fix_positions() {
    std::ranges::for_each(cells, &cell::clear);
    std::ranges::for_each(arena.range_entries(), [this](container::arena<Particle>::entry &p) {
      insert_into_cell(p);
    });
    SPDLOG_TRACE("fixed positions");
  }

 private:
  /**
   * Inserts a reference to a particle in a arena into it's cell
   * */
  auto insert_into_cell(arena<Particle>::entry &particle) {

    auto idx = index.position_to_index(particle.data.position);
    if (idx < cells.size()) {
      cells[idx].insert(particle);
    } else {
      spdlog::warn("a particle has positions that is out of bounds {} {} {}", particle.data.position[0], particle.data.position[1], particle.data.position[2]);
    }
  }

  container::arena<Particle> arena;
  std::vector<cell> cells;
  std::array<boundary_condition, 6> bc;

 public:
  I index;
  double sigma;
  double cutoff;
};

/**
 * a function to apply all up to 6 boundary conditions to a all boundary cell of a linked container
 * */
template<index::Index I>
static void calculate_boundary_condition(linked_cell<I> &lc,
                                         std::function<void(std::tuple<Particle &, Particle &>)> const &force_calculation

) {

  std::ranges::for_each(lc.boundary(), [&lc, &force_calculation](cell &cell) {
    for (auto [side, b] : std::views::enumerate(cell.boundary)) {
      auto o = boundary::orientation(side);
      switch (b) {
        case boundary_condition::outflow:
          std::ranges::for_each(cell.particles, [&lc, &o](auto &e) { outflow(lc, e, o); });
          break;
        case boundary_condition::reflecting:
          std::ranges::for_each(cell.linear(), [&lc, &o, &force_calculation](auto &p) { reflecting(lc, p, o, force_calculation); });
          break;
        case boundary_condition::none: break;
      }
    }
  });
}

}// namespace container
