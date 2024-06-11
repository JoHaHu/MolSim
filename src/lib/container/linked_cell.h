#pragma once

#include "Particle.h"
#include "container.h"
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
#include "spdlog/spdlog.h"

namespace container {

/**
 * @brief Enum for defining different types of cells.
 */
enum class cell_type : std::uint8_t {
  inner,
  boundary,
};

class LinkedCell;

/**
   * @brief A cell in the linked cell data structure.
   *
   * Holds particles and supports operations on them.
   */
class Cell {
  friend class LinkedCell;

 private:
  struct ParticlePage {
    size_t index;
    double_mask mask;
  };

 public:
  explicit Cell(Particles &particles, cell_type type, std::array<size_t, 3> idx, std::vector<size_t> neighbours) : particles(particles), type(type), neighbours(std::move(neighbours)), idx(idx) {};

  /**
   * a pairwise range over the particles, uses a cached range initialized at the start of the program
   * */
  template<typename Callable>
  auto pairwise(Callable c) {
    for (auto &page : particle_pages) {

      for (int offset = 0; offset < double_mask::size(); ++offset) {
        auto p1 = particles.load_vectorized_single(page.index + offset);
        if (stdx::any_of(p1.active>0 )) {

        }
      }
    }
  };

  constexpr auto is_boundary() -> bool {
    return type == cell_type::boundary;
  }

  /**
   * inserts a new particle into the cell
   * */
  auto insert(size_t index, double_mask mask) {
    auto found = std::find_if(particle_pages.begin(), particle_pages.end(), [index](auto el) { return el.index == index; });
    if (found != particle_pages.end()) {
      *found = ParticlePage(index, found->mask | mask);
    }
  }

  /**
   * clears the cell
   * */
  auto clear() {
    particle_pages.clear();
  }

 public:
  Particles &particles;
  std::vector<ParticlePage> particle_pages{};
  cell_type type = cell_type::inner;
  std::vector<size_t> neighbours{};//todo correction vector
  std::array<size_t, 3> idx{};
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
class LinkedCell {
  friend class cell;

 public:
  LinkedCell() = delete;
  explicit LinkedCell(const std::array<double, 3> &domain, double cutoff, std::array<BoundaryCondition, 6> bc, double sigma)
      : particles(Particles()),
        bc(bc),
        sigma(sigma),
        cutoff(cutoff) {
    widths = {cutoff,
              cutoff,
              cutoff};

    dim = {
        (size_t) std::ceil(domain[0] / widths[0]),
        (size_t) std::ceil(domain[1] / widths[1]),
        (size_t) std::ceil(domain[2] / widths[2])};
    cells.reserve(dim[0] * dim[1] * dim[2]);

    auto range = std::views::cartesian_product(
        std::views::iota(0UL, dim[0]),
        std::views::iota(0UL, dim[1]),
        std::views::iota(0UL, dim[2]));

    for (auto [x, y, z] : range) {
      this->cells.emplace_back(particles, cell_type::inner, std::array<size_t, 3>({x, y, z}), std::vector<size_t>());
    }

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

      c.neighbours = create_neighbours({x, y, z});
    }
  };

  /**
   * Creates the range for a cell, containing the neighbours needed to calculate forces and respects Newtons 3. Law
   * Can be calculated once at the start and then be reused
   * */
  auto create_neighbours(std::array<size_t, 3> idx) -> std::vector<size_t> {
    auto cell_idx = dimension_to_index(idx);
    auto [radius_x, radius_y, radius_z] = std::array<size_t, 3>{1, 1, 1};
    auto neighbours = std::vector<size_t>();

    for (long x = -(long) radius_x; x <= (long) radius_x; ++x) {
      for (long y = -(long) radius_y; y <= (long) radius_y; ++y) {
        for (long z = -(long) radius_z; z <= (long) radius_z; ++z) {
          if (x != 0 || y != 0 || z != 0) {
            if (z > 0 || (x > 0 && z >= 0) || (x >= 0 && z >= 0 && y > 0)) {
              auto other_index = offset(idx, {x, y, z});
              // checks that cell is not out of bounds and not the same cell
              if (other_index < cells.size() && other_index != cell_idx) {
                neighbours.emplace_back(other_index);
              }
            }
          }
        }
      }
    }
    neighbours.shrink_to_fit();
    return neighbours;
  }

  constexpr auto position_to_index(std::array<double, 3> position) -> size_t {
    auto [x, y, z] = position;

    std::array<size_t, 3> dim = {(size_t) std::floor(x / widths[0]),
                                 (size_t) std::floor(y / widths[1]),
                                 (size_t) std::floor(z / widths[2])};

    return dimension_to_index(dim);
  }
  auto dimension_to_index(std::array<size_t, 3> dimension) -> size_t {
    auto [x, y, z] = dimension;
    return x + y * dim[0] + z * dim[0] * dim[1];
  }

  auto offset(std::array<size_t, 3> index, std::array<long, 3> off) -> size_t {
    auto [x, y, z] = index;
    auto [off_x, off_y, off_z] = off;
    return dimension_to_index({x + off_x, y + off_y, z + off_z});
  }

  /**
   * a range over the boundary cells
   * */
  template<typename Callable>
  auto boundary(Callable f) -> auto {
    for (auto &cell : cells) {
      if (cell.is_boundary()) {
        cell.pairwise(f);
      }
    }
  }

  /**
   * returns a range with pairs of Particles
   * */

  template<typename Callable>
  auto pairwise(Callable f) -> auto {
    for (auto &cell : cells) {
      cell.pairwise(f);
    }
  }

  /**
   * inserts a new particle into the arena and the into the cells
   * */
  auto insert(Particle particle) {
    auto inserted = particles.insert_particle(particle);
    insert_into_cell(inserted);
  }

  /**
   * Returns the number of particles in the arena. O(n), because it gets counted each time
   * */
  auto size() {
    return particles.size;
  }

  /**
   * Clear all cells and reinserts all particles into the cells.
   * Better than updating positions when changes are made, because the vectors get resized/reorderd a lot when only single elements get removed
   * */
  auto fix_positions() {
    std::ranges::for_each(cells, &Cell::clear);

    for (int i = 0; i < particles.size; ++i) {
      insert_into_cell(i);
    }

    SPDLOG_TRACE("fixed positions");
  }

  Particles particles;

 private:
  /**
   * Inserts a reference to a particle in a arena into it's cell
   * */
  auto insert_into_cell(size_t particle) -> void {

    auto cell = position_to_index({particles.position_x[particle], particles.position_y[particle], particles.position_z[particle]});
    if (cell < cells.size()) {
      auto offset = particle % double_mask ::size();
      auto mask = double_mask();
      mask[offset] = true;
      cells[cell].insert(particle - offset, mask);
    } else {
      spdlog::warn("a particle has positions that is out of bounds ");
    }
  }

  std::vector<Cell> cells;
  std::array<BoundaryCondition, 6> bc;
  std::array<double, 3> widths;
  std::array<size_t, 3> dim;

 public:
  double sigma;
  double cutoff;
};

///**
// * a function to apply all up to 6 boundary conditions to a all boundary cell of a linked container
// * */
//
//static void calculate_boundary_condition(linked_cell &lc,
//                                         std::function<void(std::tuple<Particle &, Particle &>)> const &force_calculation
//
//) {
//
//  std::ranges::for_each(lc.boundary(), [&lc, &force_calculation](cell &cell) {
//    for (auto [side, b] : std::views::enumerate(cell.boundary)) {
//      auto o = orientation(side);
//      switch (b) {
//        case BoundaryCondition::outflow:
//          std::ranges::for_each(cell.particles, [&lc, &o](auto &e) { outflow<I>(lc, e, o); });
//          break;
//        case BoundaryCondition::reflecting:
//          std::ranges::for_each(cell.linear(), [&lc, &o, &force_calculation](auto &p) { reflecting<I>(lc, p, o, force_calculation); });
//          break;
//        case BoundaryCondition::none: break;
//      }
//    }
//  });
//}
//
///**
// * function to apply outflow boundary condition
// * */
//
//static void outflow(linked_cell &lc, arena<Particle>::entry &entry, orientation o) {
//
//  auto &p = entry.data;
//  const auto [pos_x, pos_y, pos_z] = p.position;
//  const auto [bound_x, bound_y, bound_z] = lc.index.boundary;
//
//  bool condition = false;
//  switch (o) {
//    case orientation::front:
//      condition |= pos_z > bound_z;
//      break;
//    case orientation::back:
//      condition |= pos_z < 0;
//      break;
//    case orientation::left:
//      condition |= pos_x < 0;
//      break;
//    case orientation::right:
//      condition |= pos_x > bound_x;
//      break;
//    case orientation::bottom:
//      condition |= pos_y < 0;
//      break;
//    case orientation::top:
//      condition |= pos_y > bound_y;
//      break;
//  }
//  entry.active &= !condition;
//}
//
///**
// * function to apply reflecting boundary condition
// * */
//
//static void reflecting(linked_cell &lc, auto &p, orientation o, auto fc) {
//  auto [pos_x, pos_y, pos_z] = p.position;
//  auto [bound_x, bound_y, bound_z] = lc.index.boundary;
//  auto diff = lc.index.boundary - p.position;
//  auto [diff_x, diff_y, diff_z] = diff;
//
//  const double sigma = lc.sigma;
//  const double distance = std::pow(2, 1 / 6.0) * sigma;
//
//  switch (o) {
//    case orientation::front:
//      if (diff_z <= distance && diff_z > 0) {
//        auto ghost = Particle({pos_x, pos_y, bound_z + diff_z}, {0, 0, 0}, 0, 0);
//        fc({p, ghost});
//      }
//      break;
//    case orientation::back:
//      if (pos_z <= distance && pos_z > 0) {
//        auto ghost = Particle({pos_x, pos_y, -pos_z}, {0, 0, 0}, 0, 0);
//        fc({p, ghost});
//      }
//      break;
//
//    case orientation::left:
//      if (pos_x <= distance && pos_x > 0) {
//        auto ghost = Particle({-pos_x, pos_y, pos_z}, {0, 0, 0}, 0, 0);
//        fc({p, ghost});
//      }
//      break;
//    case orientation::right:
//      if (diff_x <= distance && diff_x > 0) {
//        auto ghost = Particle({bound_x + diff_x, pos_y, pos_z}, {0, 0, 0}, 0, 0);
//        fc({p, ghost});
//      }
//      break;
//    case orientation::bottom:
//      if (pos_y <= distance && pos_y > 0) {
//        auto ghost = Particle({pos_x, -pos_y, pos_z}, {0, 0, 0}, 0, 0);
//        fc({p, ghost});
//      }
//      break;
//    case orientation::top:
//      if (diff_y <= distance && diff_y > 0) {
//        auto ghost = Particle({pos_x, bound_y + diff_y, pos_z}, {0, 0, 0}, 0, 0);
//        fc({p, ghost});
//      }
//      break;
//  }
//}

}// namespace container