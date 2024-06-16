#pragma once

#include "Particle.h"
#include "utils/ArrayUtils.h"

#include <algorithm>
#include <array>
#include <concepts>
#include <memory>
#include <ranges>
#include <utility>
#include <vector>

#include "boundary.h"
#include "index.h"
#include "spdlog/spdlog.h"

namespace container {

/**
 * @brief Enum for defining different types of cells.
 */
enum class cell_type : std::uint8_t {
  /**
   * inner cells without boundaries
   * */
  inner,
  /**
   * boundary cells
   * */
  boundary,
};

/**
 * A struct representing one neighbour. The corrections is a vector that needs to be applied to positions in this cell when using periodic boundary conditions
 * */
template<const size_t DIMENSIONS>
struct Neighbour {
  /**
   * index of the represented neighbour cell
   * */
  size_t cell;
  /**
   * correction vector to apply to account for periodic boundary conditions
   * */
  std::array<double, DIMENSIONS> correction;

  // rotation to apply to particles if periodic boundaries are not parallel to each other
  //  std::array<double, 3> rotation;

  explicit Neighbour(size_t neighbour_index, std::array<double, DIMENSIONS> correction) : cell(neighbour_index), correction(correction) {}
};

/**
   * @brief A cell in the linked cell data structure.
   *
   * Holds particles and supports operations on them.
   */
template<const size_t DIMENSIONS>
class Cell {

 private:
 public:
  explicit Cell(cell_type type, std::array<size_t, DIMENSIONS> idx) : type(type), idx(idx) {};

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
  std::vector<Neighbour<DIMENSIONS>> neighbours{};
  std::array<size_t, DIMENSIONS> idx{};
  std::array<BoundaryCondition, 2 * DIMENSIONS> boundary = initialize_boundary();

  static auto initialize_boundary() -> std::array<BoundaryCondition, 2 * DIMENSIONS> {
    std::array<BoundaryCondition, 2 * DIMENSIONS> bc;
    for (int i = 0; i < 2 * DIMENSIONS; ++i) {
      bc[i] = BoundaryCondition::none;
    }
    return bc;
  }
};

/**
 * The linked cell container class
 * the Index template parameter can be used to change the iteration order of the LinkedCell.
 * */
template<const size_t DIMENSIONS>
class LinkedCell {

 public:
  LinkedCell() = delete;

  constexpr void recursive_fill(std::array<size_t, DIMENSIONS> &coords, size_t depth) {
    if (depth == 0) {
      this->cells.emplace_back(cell_type::inner, std::array<size_t, DIMENSIONS>(coords));
    } else {
      for (size_t i = 0; i < index.dim[depth - 1]; ++i) {
        coords[depth - 1] = i;
        recursive_fill(coords, depth - 1);
      }
    }
  }

  explicit LinkedCell(const std::array<double, DIMENSIONS> &domain, double cutoff, std::array<BoundaryCondition, 2 * DIMENSIONS> bc, double sigma)
      : particles(Particles<DIMENSIONS>()),
        index(container::index::Index<DIMENSIONS>(domain, bc, cutoff)),
        sigma(sigma),
        cutoff(cutoff) {

    cells.reserve(std::ranges::fold_left(index.dim, 1.0, std::multiplies<>()));

    std::array<size_t, DIMENSIONS> temp{};
    recursive_fill(temp, DIMENSIONS);

    for (auto &c : cells) {
      for (int i = 0; i < DIMENSIONS; ++i) {
        if (c.idx[i] == 0) {
          c.type = cell_type::boundary;
          c.boundary[i] = bc[i];
        }
        if (c.idx[i] == index.dim[i] - 1) {
          c.type = cell_type::boundary;
          c.boundary[i + DIMENSIONS] = bc[i + DIMENSIONS];
        }
      }

      c.neighbours = create_neighbours(c.idx);
    }
  };

  constexpr void recursive_fill_neighbours(size_t depth, std::vector<Neighbour<DIMENSIONS>> &neighbours, std::array<size_t, DIMENSIONS> idx, std::array<long, DIMENSIONS> &offset) {

    if (depth == 0) {
      auto cell_idx = index.dimension_to_index(idx);

      if (std::ranges::any_of(offset, [](auto a) { return a != 0; })) {

        if (std::ranges::any_of(std::views::iota(0UL, DIMENSIONS), [&](auto i) {
              bool condition = offset[i] > 0;

              for (int j = i + 1; j < DIMENSIONS; ++j) {
                condition &= offset[j] >= 0;
              }
              return condition;
            })) {

          auto other_index = index.offset(idx, offset);
          auto correction = index.calculate_correction(idx, offset);
          // checks that cell is not out of bounds and not the same cell
          if (other_index < cells.size() && other_index != cell_idx && index.in_cutoff_distance(idx, cells[other_index].idx)) {
            neighbours.emplace_back(other_index, correction);
          }
        }
      }

    } else {
      for (long i = -index.radius[depth - 1]; i <= index.radius[depth - 1]; ++i) {
        offset[depth - 1] = i;
        recursive_fill_neighbours(depth - 1, neighbours, idx, offset);
      }
    }
  }

  /**
   * Creates the range for a cell, containing the neighbours needed to calculate forces and respects Newtons 3. Law
   * Can be calculated once at the start and then be reused
   * */
  constexpr auto create_neighbours(std::array<size_t, DIMENSIONS> idx) -> std::vector<Neighbour<DIMENSIONS>> {

    auto neighbours = std::vector<Neighbour<DIMENSIONS>>();
    std::array<long, DIMENSIONS> temp;
    recursive_fill_neighbours(DIMENSIONS, neighbours, idx, temp);
    neighbours.shrink_to_fit();
    std::ranges::sort(neighbours, [](auto &n1, auto &n2) {
      return n1.cell < n2.cell;
    });
    return neighbours;
  }

  /**
   * a range over the boundary cells
   * */
  template<typename Callable>
  constexpr auto boundary(Callable f) -> auto {

    for (size_t idx = 0; idx < particles.size; idx++) {
      if (particles.active[idx]) [[likely]] {
        const auto cell_idx = particles.cell[idx];
        const auto &cell = cells[cell_idx];
        if (cell.is_boundary()) {
          for (auto [side, bc] : std::ranges::enumerate_view(cell.boundary)) {

            size_t axis = side % DIMENSIONS;
            bool start_of_axis = side / DIMENSIONS == 0;
            switch (bc) {
              case BoundaryCondition::outflow:
                outflow(idx, axis, start_of_axis);
                break;
              case BoundaryCondition::reflecting:
                reflecting(idx, axis, start_of_axis, f);
                break;
              case BoundaryCondition::periodic:
                periodic(idx, axis);
                break;
              case BoundaryCondition::none: break;
            }
          }
        }
      }
    }
  }
  /**
   *
   * @param start_of_axis true is start of axis, false is end of axis
   * */
  constexpr void outflow(size_t idx, size_t axis, bool start_of_axis) {
    bool condition = false;
    if (start_of_axis) {
      const auto &pos = particles.positions[axis][idx];
      condition |= pos < 0;
    } else {
      const auto &pos = particles.positions[axis][idx];
      condition |= pos > index.domain[axis];
    }
    particles.active[idx] &= static_cast<int>(!condition);
  }

  /**
   * Assumes particles are stable and are never multiple periods outside the boundary
   *
   * */
  constexpr void periodic(size_t idx, size_t axis) {
    const auto &pos = particles.positions[axis][idx];
    if (pos < 0) {
      particles.positions[axis][idx] = index.domain[axis] - pos;
    } else if (pos > index.domain[axis]) {
      particles.positions[axis][idx] = pos - index.domain[axis];
    }
  }

  /**
   * function to apply reflecting boundary condition
   * */
  template<typename Callable>
  constexpr void reflecting(size_t idx, size_t axis, bool start_of_axis, Callable f) {
    const auto pos = particles.positions[axis][idx];
    double bound = index.domain[axis];
    auto diff = bound - pos;
    if (start_of_axis) {
      if (pos <= reflecting_distance && pos > 0) {
        particles.forces[axis][idx] -= f(2 * pos);
      }
    } else {
      if (diff <= reflecting_distance && diff > 0) {
        particles.forces[axis][idx] += f(2 * diff);
      }
    }
  }

  /**
   * a pairwise range over the particles, uses a cached range initialized at the start of the program
   * */
  // TODO align all accesses by applying proper masks
  template<typename Callable>
  constexpr auto pairwise(Callable c) {
    for (size_t idx = 0; idx < particles.size; ++idx) {
      if (particles.active[idx]) [[likely]] {
        const auto temp = particles.cell[idx];
        auto &cell = cells[temp];
        auto p1 = particles.load_vectorized_single(idx);

        // Same cell particles
        size_t i = 0;
        while (idx + 1 + (i * double_v::size()) < cell.end_index) {
          auto active_mask = index_vector + 1 + (i * double_v::size()) < (cell.size());
          auto p2 = particles.load_vectorized(idx + 1 + (i * double_v::size()));
          auto mask = stdx::static_simd_cast<double_v>(active_mask) && p2.active;
          c(p1, p2, mask, empty_correction);
          particles.store_force_vector(p2, idx + 1 + (i * double_v::size()));
          i++;
        }

        // Neighbour cells
        std::array<double_v, DIMENSIONS> correction;
        for (Neighbour<DIMENSIONS> &neighbour : cell.neighbours) {
          size_t i = 0;
          Cell<DIMENSIONS> &neighbour_cell = cells[neighbour.cell];
          size_t n_start = neighbour_cell.start_index;
          size_t n_end = neighbour_cell.end_index;

          for (int i = 0; i < DIMENSIONS; ++i) {
            correction[i] = double_v(neighbour.correction[i]);
          }

          while (n_start + (i * double_v::size()) < n_end) {
            auto active_mask = index_vector + (i * double_v::size()) < (neighbour_cell.size());
            auto p2 = particles.load_vectorized(n_start + (i * double_v::size()));
            auto mask = stdx::static_simd_cast<double_v>(active_mask) && p2.active;
            c(p1, p2, mask, correction);
            particles.store_force_vector(p2, n_start + (i * double_v::size()));
            i++;
          }
        }
        particles.store_force_single(p1, idx);
      } else [[unlikely]] {
        // Stop on first inactive encountered particle
        SPDLOG_WARN("Particle crossed boundary not handled before next loop");
        break;
      }
    }
  }

  /**
   * inserts a new particle into the arena and the into the cells
   * */
  constexpr auto insert(Particle<DIMENSIONS> particle) {
    particles.insert_particle(particle);
  }

  /**
   * Returns the number of particles in the arena. O(n), because it gets counted each time
   * */
  constexpr auto size() {
    return particles.size;
  }

  /**
   * Clear all cells and reinserts all particles into the cells.
   * Better than updating positions when changes are made, because the vectors get resized/reorderd a lot when only single elements get removed
   * */
  constexpr auto fix_positions() {

    for (auto &cell : cells) {
      cell.start_index = 0;
      cell.end_index = 0;
    }

    particles.sort([this](size_t idx) -> size_t {
      std::array<double, DIMENSIONS> temp;
      for (int i = 0; i < DIMENSIONS; ++i) {
        temp[i] = particles.positions[i][idx];
      }
      return index.position_to_index(temp);
    });

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
      if (i == particles.size - 1) {
        cells[cell].start_index = start;
        cells[cell].end_index = particles.size;
      }
    }

    SPDLOG_TRACE("fixed positions");
  }

  constexpr static auto init_index_vector() noexcept -> size_v {
    size_v index_vector_tmp = 0;
    for (size_t i = 0; i < size_v::size(); ++i) {
      index_vector_tmp[i] = i;
    }
    return index_vector_tmp;
  }

  Particles<DIMENSIONS> particles;

 private:
  /**
   * Cells with start and end offset of particles and the neighbour cells
   * */
  std::vector<Cell<DIMENSIONS>> cells;
  /**
   * an index helper
   * */
  container::index::Index<DIMENSIONS> index;
  /**
   * a vector (AVX) with stepwise incremented
   * */
  size_v index_vector = init_index_vector();
  double sigma;
  double cutoff;
  double reflecting_distance = std::pow(2, 1 / 6.0) * sigma;

  constexpr static auto init_empty_correction() -> std::array<double_v, DIMENSIONS> {
    std::array<double_v, DIMENSIONS> temp{};
    for (int i = 0; i < DIMENSIONS; ++i) {
      temp[i] = double_v(0.0);
    }
    return temp;
  }

  std::array<double_v, DIMENSIONS> empty_correction = init_empty_correction();
};

}// namespace container