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

  struct Neighbour {
    size_t c;
    // TODO correction vector
    explicit Neighbour(size_t neighbour_index) : c(neighbour_index) {}
  };

 private:
 public:
  explicit Cell(cell_type type, std::array<size_t, 3> idx) : type(type), idx(idx) {};

  constexpr auto is_boundary() const -> bool {
    return type == cell_type::boundary;
  }

  auto size() -> size_t {
    return end_index - start_index;
  }

 public:
  size_t start_index;
  size_t end_index;
  cell_type type = cell_type::inner;
  std::vector<Neighbour> neighbours{};
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
        domain(domain),
        bc(bc),
        sigma(sigma),
        cutoff(cutoff) {
    widths = {cutoff / 2,
              cutoff / 2,
              cutoff / 2};

    dim = {
        (size_t) std::ceil(domain[0] / widths[0]),
        (size_t) std::ceil(domain[1] / widths[1]),
        (size_t) std::ceil(domain[2] / widths[2])};
    cells.reserve(dim[0] * dim[1] * dim[2]);

    // Iterate in reverse order to get correct index
    for (size_t z = 0; z < dim[2]; ++z) {
      for (size_t y = 0; y < dim[1]; ++y) {
        for (size_t x = 0; x < dim[0]; ++x) {
          this->cells.emplace_back(cell_type::inner, std::array<size_t, 3>({x, y, z}));
        }
      }
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
  auto create_neighbours(std::array<size_t, 3> idx) -> std::vector<Cell::Neighbour> {
    auto cell_idx = dimension_to_index(idx);
    auto [radius_x, radius_y, radius_z] = std::array<size_t, 3>{2, 2, 2};
    auto neighbours = std::vector<Cell::Neighbour>();

    for (long x = -(long) radius_x; x <= (long) radius_x; ++x) {
      for (long y = -(long) radius_y; y <= (long) radius_y; ++y) {
        for (long z = -(long) radius_z; z <= (long) radius_z; ++z) {
          if (x != 0 || y != 0 || z != 0) {
            if (z > 0 || (x > 0 && z >= 0) || (x >= 0 && z >= 0 && y > 0)) {
              auto other_index = offset(idx, {x, y, z});
              // checks that cell is not out of bounds and not the same cell
              if (other_index < cells.size() && other_index != cell_idx && min_distance(idx, cells[other_index].idx) <= cutoff) {
                neighbours.emplace_back(other_index);
              }
            }
          }
        }
      }
    }
    neighbours.shrink_to_fit();
    std::sort(neighbours.begin(), neighbours.end(), [](auto &n1, auto &n2) {
      return n1.c < n2.c;
    });
    return neighbours;
  }

  auto min_distance(std::array<size_t, 3> dim1, std::array<size_t, 3> dim2) -> double {
    auto distances = std::vector<double>();

    for (auto [x, y, z] : std::views::cartesian_product(
             std::views::iota(0, 2),
             std::views::iota(0, 2),
             std::views::iota(0, 2))) {
      distances.emplace_back(
          distance({dim1[0] + x, dim1[1] + y, dim1[2] + z}, dim2));
    }

    double min = std::ranges::min(distances);
    return min;
  }
  /**
 * @brief Calculates the reflecting_distance between two sets of 3D coordinates.
 *
 * @param dim1 The first set of coordinates.
 * @param dim2 The second set of coordinates.
 * @return The reflecting_distance.
 */
  auto distance(std::array<size_t, 3> dim1, std::array<size_t, 3> dim2) -> double {
    auto diff = std::array<double, 3>({
        (double) std::abs(static_cast<long>(dim1[0]) - static_cast<long>(dim2[0])) * widths[0],
        (double) std::abs(static_cast<long>(dim1[1]) - static_cast<long>(dim2[1])) * widths[1],
        (double) std::abs(static_cast<long>(dim1[2]) - static_cast<long>(dim2[2])) * widths[2],
    });
    return ArrayUtils::L2Norm(diff);
  }

  constexpr auto position_to_index(std::array<double, 3> position) -> size_t {
    auto [x, y, z] = position;

    std::array<size_t, 3> dimension = {(size_t) std::floor(x / widths[0]),
                                       (size_t) std::floor(y / widths[1]),
                                       (size_t) std::floor(z / widths[2])};

    return dimension_to_index(dimension);
  }
  auto dimension_to_index(std::array<size_t, 3> dimension) -> size_t {
    auto [x, y, z] = dimension;
    if (x >= dim[0] || y >= dim[1] || z >= dim[2]) {
      return ULONG_MAX;
    }
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

    for (size_t index = 0; index < particles.size; index++) {
      if (particles.active[index]) [[likely]] {
        const auto cell_idx = particles.cell[index];
        const auto &cell = cells[cell_idx];
        if (cell.is_boundary()) {
          for (auto [side, bc] : std::ranges::enumerate_view(cell.boundary)) {

            auto o = orientation(side);
            switch (bc) {
              case BoundaryCondition::outflow:
                outflow(index, o);
                break;
              case BoundaryCondition::reflecting:
                reflecting(index, o, f);
                break;
              case BoundaryCondition::none:
              case BoundaryCondition::periodic: break;
            }
          }
        }
      }
    }
  }

  void outflow(size_t index, orientation o) {

    const auto pos_x = particles.position_x[index];
    const auto pos_y = particles.position_y[index];
    const auto pos_z = particles.position_z[index];
    const auto [bound_x, bound_y, bound_z] = domain;

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
    particles.active[index] &= static_cast<int>(!condition);
  }

  /**
   * function to apply reflecting boundary condition
   * */
  template<typename Callable>
  void reflecting(size_t index, orientation o, Callable f) {
    const auto pos_x = particles.position_x[index];
    const auto pos_y = particles.position_y[index];
    const auto pos_z = particles.position_z[index];
    const auto [bound_x, bound_y, bound_z] = domain;

    auto diff_x = bound_x - pos_x;
    auto diff_y = bound_y - pos_y;
    auto diff_z = bound_z - pos_z;

    switch (o) {
      case orientation::front:
        if (diff_z <= reflecting_distance && diff_z > 0) {
          std::array<double, 3> ghost = {pos_x, pos_y, bound_z + diff_z};
          f(particles, index, ghost);
        }
        break;
      case orientation::back:
        if (pos_z <= reflecting_distance && pos_z > 0) {
          std::array<double, 3> ghost = {pos_x, pos_y, -pos_z};
          f(particles, index, ghost);
        }
        break;

      case orientation::left:
        if (pos_x <= reflecting_distance && pos_x > 0) {
          std::array<double, 3> ghost = {-pos_x, pos_y, pos_z};
          f(particles, index, ghost);
        }
        break;
      case orientation::right:
        if (diff_x <= reflecting_distance && diff_x > 0) {
          std::array<double, 3> ghost = {bound_x + diff_x, pos_y, pos_z};
          f(particles, index, ghost);
        }
        break;
      case orientation::bottom:
        if (pos_y <= reflecting_distance && pos_y > 0) {
          std::array<double, 3> ghost = {pos_x, -pos_y, pos_z};
          f(particles, index, ghost);
        }
        break;
      case orientation::top:
        if (diff_y <= reflecting_distance && diff_y > 0) {
          std::array<double, 3> ghost = {pos_x, bound_y + diff_y, pos_z};
          f(particles, index, ghost);
        }
        break;
    }
  }

  /**
   * a pairwise range over the particles, uses a cached range initialized at the start of the program
   * */
  // TODO align all accesses by applying proper masks
  template<typename Callable>
  auto pairwise(Callable c) {
    for (size_t index = 0; index < particles.size; ++index) {
      if (particles.active[index]) [[likely]] {
        const auto temp = particles.cell[index];
        auto &cell = cells[temp];
        auto p1 = particles.load_vectorized_single(index);

        // Same cell particles
        size_t i = 0;
        while (index + 1 + (i * double_v::size()) < cell.end_index) {
          auto active_mask = index_vector + 1 + (i * double_v::size()) < (cell.size());
          auto p2 = particles.load_vectorized(index + 1 + (i * double_v::size()));
          auto mask = stdx::static_simd_cast<double_v>(active_mask) && p2.active;
          c(p1, p2, mask);
          particles.store_force_vector(p2, index + 1 + (i * double_v::size()), mask);
          i++;
        }

        // Neighbour cells
        for (auto &neighbour : cell.neighbours) {
          size_t i = 0;
          auto &neighbour_cell = cells[neighbour.c];
          size_t n_start = neighbour_cell.start_index;
          size_t n_end = neighbour_cell.end_index;
          while (n_start + (i * double_v::size()) < n_end) {
            auto active_mask = index_vector + (i * double_v::size()) < (neighbour_cell.size());
            auto p2 = particles.load_vectorized(n_start + (i * double_v::size()));
            auto mask = stdx::static_simd_cast<double_v>(active_mask) && p2.active;
            c(p1, p2, mask);
            particles.store_force_vector(p2, n_start + (i * double_v::size()), mask);
            i++;
          }
        }
        particles.store_force_single(p1, index);
      } else [[unlikely]] {
        // Stop on first inactive encountered particle
        break;
        //        SPDLOG_WARN("Particle crossed boundary not handled before next loop");
      }
    }
  }

  /**
   * inserts a new particle into the arena and the into the cells
   * */
  auto insert(Particle particle) {
    particles.insert_particle(particle);
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

    for (auto &cell : cells) {
      cell.start_index = 0;
      cell.end_index = 0;
    }

    particles.sort([this](std::array<double, 3> idx) -> size_t { return position_to_index(idx); });

    size_t start = 0;
    size_t cell = particles.cell[0];

    for (size_t i = 0; i < particles.size; ++i) {
      if (particles.cell[i] >= cells.size()) {
        particles.active[i] = false;
      }
      if (particles.cell[i] != cell) {
        if (cell < cells.size()) {
          cells[cell].start_index = start;
          cells[cell].end_index = i;
        }
        start = i;
        cell = particles.cell[i];
      }
    }

    SPDLOG_TRACE("fixed positions");
  }

  static auto init_index_vector() noexcept -> size_v {
    size_v index_vector_tmp = 0;
    for (size_t i = 0; i < size_v::size(); ++i) {
      index_vector_tmp[i] = i;
    }
    return index_vector_tmp;
  }

  Particles particles;

 private:
  std::vector<Cell> cells;
  std::array<double, 3> domain;
  std::array<BoundaryCondition, 6> bc;
  std::array<double, 3> widths;
  std::array<size_t, 3> dim;
  size_v index_vector = init_index_vector();
  double sigma;
  double cutoff;
  double reflecting_distance = std::pow(2, 1 / 6.0) * sigma;
};

}// namespace container