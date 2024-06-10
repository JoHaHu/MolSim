#pragma once

#include "Particle.h"
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

#include "boundary.h"
#include "orientation.h"
#include "range/v3/view/concat.hpp"
#include "spdlog/spdlog.h"

namespace container {

/**
 * @brief Enum for defining different types of cells.
 */
enum class cell_type : std::uint8_t {
  inner,
  boundary,
};

template<index::Index I>
class linked_cell;

/**
   * @brief A cell in the linked cell data structure.
   *
   * Holds particles and supports operations on them.
   */
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
  std::array<BoundaryCondition, 6> boundary = {
      BoundaryCondition::none,
      BoundaryCondition::none,
      BoundaryCondition::none,
      BoundaryCondition::none,
      BoundaryCondition::none,
      BoundaryCondition::none};
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
  explicit linked_cell(const std::array<double, 3> &domain, double cutoff, std::array<BoundaryCondition, 6> bc, unsigned int number_of_particles, double sigma)
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
        c.boundary[(size_t) orientation::left] = bc[(size_t) orientation::left];
      }
      if (x == dim[0] - 1) {
        c.type = cell_type::boundary;
        c.boundary[(size_t) orientation::right] = bc[(size_t) orientation::right];
      }
      if (y == 0) {
        c.type = cell_type::boundary;
        c.boundary[(size_t) orientation::bottom] = bc[(size_t) orientation::bottom];
      }
      if (y == dim[1] - 1) {
        c.type = cell_type::boundary;
        c.boundary[(size_t) orientation::top] = bc[(size_t) orientation::top];
      }
      if (z == 0) {
        c.type = cell_type::boundary;
        c.boundary[(size_t) orientation::back] = bc[(size_t) orientation::back];
      }
      if (z == dim[2] - 1) {
        c.type = cell_type::boundary;
        c.boundary[(size_t) orientation::front] = bc[(size_t) orientation::front];
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
  auto insert(Particle &&particle) {
    auto &inserted = arena.emplace_back(std::move(particle));
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
  std::array<BoundaryCondition, 6> bc;

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
      auto o = orientation(side);
      switch (b) {
        case BoundaryCondition::outflow:
          std::ranges::for_each(cell.particles, [&lc, &o](auto &e) { outflow<I>(lc, e, o); });
          break;
        case BoundaryCondition::reflecting:
          std::ranges::for_each(cell.linear(), [&lc, &o, &force_calculation](auto &p) { reflecting<I>(lc, p, o, force_calculation); });
          break;
        case BoundaryCondition::none: break;
      }
    }
  });
}

/**
 * function to apply outflow boundary condition
 * */
template<index::Index I>
static void outflow(linked_cell<I> &lc, arena<Particle>::entry &entry, orientation o) {

  auto &p = entry.data;
  const auto [pos_x, pos_y, pos_z] = p.position;
  const auto [bound_x, bound_y, bound_z] = lc.index.boundary;

  bool condition = false;
  switch (o) {
    case orientation::front:
      condition |= pos_z > bound_z;
      break;
    case orientation::back:
      condition |= pos_z < 0;
      break;
    case orientation::left:
      condition |= pos_x < 0;
      break;
    case orientation::right:
      condition |= pos_x > bound_x;
      break;
    case orientation::bottom:
      condition |= pos_y < 0;
      break;
    case orientation::top:
      condition |= pos_y > bound_y;
      break;
  }
  entry.active &= !condition;
}

/**
 * function to apply reflecting boundary condition
 * */

template<index::Index I>
static void reflecting(linked_cell<I> &lc, auto &p, orientation o, auto fc) {
  auto [pos_x, pos_y, pos_z] = p.position;
  auto [bound_x, bound_y, bound_z] = lc.index.boundary;
  auto diff = lc.index.boundary - p.position;
  auto [diff_x, diff_y, diff_z] = diff;

  const double sigma = lc.sigma;
  const double distance = std::pow(2, 1 / 6.0) * sigma;

  switch (o) {
    case orientation::front:
      if (diff_z <= distance && diff_z > 0) {
        auto ghost = Particle({pos_x, pos_y, bound_z + diff_z}, {0, 0, 0}, 0, 0);
        fc({p, ghost});
      }
      break;
    case orientation::back:
      if (pos_z <= distance && pos_z > 0) {
        auto ghost = Particle({pos_x, pos_y, -pos_z}, {0, 0, 0}, 0, 0);
        fc({p, ghost});
      }
      break;

    case orientation::left:
      if (pos_x <= distance && pos_x > 0) {
        auto ghost = Particle({-pos_x, pos_y, pos_z}, {0, 0, 0}, 0, 0);
        fc({p, ghost});
      }
      break;
    case orientation::right:
      if (diff_x <= distance && diff_x > 0) {
        auto ghost = Particle({bound_x + diff_x, pos_y, pos_z}, {0, 0, 0}, 0, 0);
        fc({p, ghost});
      }
      break;
    case orientation::bottom:
      if (pos_y <= distance && pos_y > 0) {
        auto ghost = Particle({pos_x, -pos_y, pos_z}, {0, 0, 0}, 0, 0);
        fc({p, ghost});
      }
      break;
    case orientation::top:
      if (diff_y <= distance && diff_y > 0) {
        auto ghost = Particle({pos_x, bound_y + diff_y, pos_z}, {0, 0, 0}, 0, 0);
        fc({p, ghost});
      }
      break;
  }
}

}// namespace container