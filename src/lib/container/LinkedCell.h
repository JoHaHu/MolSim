#pragma once

#include "Particle.h"
#include "simulator/physics/Force.h"
#include "utils/ArrayUtils.h"

#include <algorithm>
#include <array>
#include <concepts>
#include <memory>
#include <omp.h>
#include <utility>
#include <vector>

#include "boundary.h"
#include "index.h"
#include "spdlog/spdlog.h"

#include "Container.h"
#include "range/v3/algorithm.hpp"
#include "range/v3/view.hpp"
#include "simulator/physics/Force.h"

namespace container {

/**
 * @brief Enum for defining different types of cells.
 */
enum class CellType : std::uint8_t {
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
 public:
  /**
   * index of the represented neighbour cell
   * */
  size_t cell{};
  /**
   * correction vector to apply to account for periodic boundary conditions
   * */
  std::array<double, DIMENSIONS> correction{};

  // rotation to apply to particles if periodic boundaries are not parallel to each other
  //  std::array<double, 3> rotation;

  Neighbour() = default;

  explicit Neighbour(size_t neighbour_index, std::array<double, DIMENSIONS> correction) : cell(neighbour_index), correction(correction) {}
};

/**
 * A block referencing the start and end index of the particles of it's cells
 * */
struct Block {
  size_t start_index;
  size_t end_index;
  Block(size_t start_index, size_t end_index) : start_index(start_index), end_index(end_index) {}
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
  explicit Cell(CellType type, std::array<size_t, DIMENSIONS> idx) : type(type), idx(idx) {};

  constexpr auto is_boundary() const -> bool {
    return type == CellType::boundary;
  }

  auto size() const -> size_t {
    return end_index - start_index;
  }

 public:
  /**
   * start index of the cell
   * */
  size_t start_index;
  /**
   * end  index of the particles
   * */
  size_t end_index;

  CellType type = CellType::inner;
  /**
   * All neighbours
   * */
  std::vector<Neighbour<DIMENSIONS>> neighbours{};
  std::array<size_t, DIMENSIONS> idx{};
  std::array<BoundaryCondition, 2 * DIMENSIONS> boundary = initialize_boundary();

  /**
   * The color of the cell, representing in which partition of the parralellization it runs
   * */
  size_t cell_color = 0;
  /**
   * the block id with neighbouring cells with same color
   * */
  size_t block = 0;

  static auto initialize_boundary() -> std::array<BoundaryCondition, 2 * DIMENSIONS> {
    std::array<BoundaryCondition, 2 * DIMENSIONS> bc;
    for (size_t i = 0; i < 2 * DIMENSIONS; ++i) {
      bc[i] = BoundaryCondition::none;
    }
    return bc;
  }
};

/**
 * The linked cell container class
 * */
template<const size_t DIMENSIONS>
class LinkedCell {

 public:
  LinkedCell() = delete;

  constexpr void recursive_fill(std::array<size_t, DIMENSIONS> &coords, size_t depth) {
    if (depth == 0) {
      this->cells.emplace_back(CellType::inner, std::array<size_t, DIMENSIONS>(coords));
    } else {
      for (size_t i = 0; i < index.dim[depth - 1]; ++i) {
        coords[depth - 1] = i;
        recursive_fill(coords, depth - 1);
      }
    }
  }

  explicit LinkedCell(
      F &&force,
      const std::array<double, DIMENSIONS> &domain,
      double cutoff,
      std::array<BoundaryCondition, 2 * DIMENSIONS> bc,
      std::vector<double> &sigma)
      : particles(Particles<DIMENSIONS>()),
        force(std::move(force)),
        index(container::index::Index<DIMENSIONS>(domain, bc, cutoff)),
        sigma(sigma) {

    size_t temp_size = 1.0;
    for (size_t i = 0; i < DIMENSIONS; ++i) {
      temp_size *= index.dim[i];
    }

    reflecting_distance = std::vector<double>();
    for (size_t i = 0; i < sigma.size(); ++i) {
      reflecting_distance.emplace_back(std::pow(2, 1 / 6.0) * sigma[i]);
    }

    std::array<size_t, DIMENSIONS> temp{};
    recursive_fill(temp, DIMENSIONS);

    for (auto &c : cells) {
      for (size_t i = 0; i < DIMENSIONS; ++i) {
        if (c.idx[i] == 0) {
          c.type = CellType::boundary;
          c.boundary[i] = bc[i];
        }
        if (c.idx[i] == index.dim[i] - 1) {
          c.type = CellType::boundary;
          c.boundary[i + DIMENSIONS] = bc[i + DIMENSIONS];
        }
      }
      c.cell_color = index.color(c.idx);
      c.block = index.block_id(c.idx);

      c.neighbours = create_neighbours(c.idx);
    }
  };

  /**
   * calculates all neighbours and the correction vector if a axis is a periodic axis.
   * The correction vector gets added to all calculations between particles of those 2 cells and
   * accounts for the boundary conditions.
   *
   * */
  constexpr void recursive_fill_neighbours(size_t depth, std::vector<Neighbour<DIMENSIONS>> &neighbours, std::array<size_t, DIMENSIONS> idx, std::array<long, DIMENSIONS> &offset) {

    if (depth == 0) {
      auto cell_idx = index.dimension_to_index(idx);

      if (ranges::any_of(offset, [](auto a) { return a != 0; })) {

        if (ranges::any_of(ranges::iota_view(0UL, DIMENSIONS), [&](auto i) {
              bool condition = offset[i] > 0;

              for (size_t j = i + 1; j < DIMENSIONS; ++j) {
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
   * a applies the boundary conditions
   * */
  void boundary() override {

    // Keep the particle size outside the loop to improve performance and to ensure that when deleting particles, still all particles are interated
    const size_t particle_size = particles.size;
    for (size_t idx = 0; idx < particle_size; idx++) {
      if (particles.active[idx]) [[likely]] {
        const auto cell_idx = particles.cell[idx];
        const auto &cell = cells[cell_idx];
        if (cell.is_boundary()) {
          for (size_t i = 0; i < DIMENSIONS * 2; ++i) {
            auto side = i;
            auto bc = cell.boundary[i];

            size_t axis = side % DIMENSIONS;
            bool start_of_axis = side / DIMENSIONS == 0;
            switch (bc) {
              case BoundaryCondition::outflow:
                outflow(idx, axis, start_of_axis);
                break;
              case BoundaryCondition::reflecting:
                reflecting(idx, axis, start_of_axis);
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
    if (!particles.active[idx]) {
      particles.size--;
    }
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
  constexpr void reflecting(size_t idx, size_t axis, bool start_of_axis) {
    const auto pos = particles.positions[axis][idx];
    double bound = index.domain[axis];
    auto diff = bound - pos;
    if (start_of_axis) {
      if (2 * pos < reflecting_distance[particles.type[idx]] && pos > 0) {
        particles.forces[axis][idx] -= force.calculate_boundary_force(2 * pos, particles.type[idx]);
      }
    } else {
      if (2 * diff < reflecting_distance[particles.type[idx]] && diff > 0) {
        particles.forces[axis][idx] += force.calculate_boundary_force(2 * diff, particles.type[idx]);
      }
    }
  }

  /**
   * a pairwise range over the particles, uses a cached range initialized at the start of the program
   * Linear scan over the particles.  They are lined out contiguously in memory.
   * for each particles it's cell is loaded and the forces are calculated between it and all within the same cell.
   * After that the same is done for the neighbour cells.
   * All calculations are done vectorized.
   * After that this process is repeated for the next particle.
   *
   * Prerequisite for this is that all particles are sorted by cell.
   *
   * */

  void pairwise(bool parallel, bool vectorize) override {

    for (auto &blocks : colored_blocks) {
#pragma omp parallel for default(none) shared(blocks, vectorize) schedule(dynamic) if (parallel)
      for (auto block : blocks) {
        for (size_t idx = block.start_index; idx < block.end_index; ++idx) {
          if (particles.active[idx] == 1) [[likely]] {
            const auto temp = particles.cell[idx];
            auto &cell = cells[temp];

            const size_t end_index = cell.end_index;

            const auto x1 = particles.positions[0][idx];
            const auto y1 = particles.positions[1][idx];
            double z1 = 0;
            if constexpr (DIMENSIONS == 3) {
              z1 = particles.positions[2][idx];
            }
            const auto mass1 = particles.mass[idx];
            const auto type1 = particles.type[idx];

            double force_sum[DIMENSIONS] = {0};

            size_t other_idx = idx + 1;
// Same cell particles
#pragma omp simd linear(other_idx) simdlen(8) reduction(+ : force_sum[ : DIMENSIONS]) if (vectorize)
            for (other_idx = idx + 1; other_idx < end_index; ++other_idx) {

              if (particles.active[other_idx] == 1) {
                double result_x = 0;
                double result_y = 0;
                double result_z = 0;

                if constexpr (DIMENSIONS == 2) {
                  auto correction = std::array<double, 2>({0, 0});
                  force.calculateForce_2D(
                      x1, y1, mass1, type1,
                      particles.positions[0][other_idx],
                      particles.positions[1][other_idx],
                      particles.mass[other_idx],
                      particles.type[other_idx],
                      result_x,
                      result_y,
                      correction);
                } else {
                  auto correction = std::array<double, 3>({0, 0, 0});
                  force.calculateForce_3D(
                      x1, y1, z1, mass1, type1,
                      particles.positions[0][other_idx],
                      particles.positions[1][other_idx],
                      particles.positions[2][other_idx],
                      particles.mass[other_idx],
                      particles.type[other_idx],
                      result_x,
                      result_y,
                      result_z,
                      correction);
                }
                particles.forces[0][other_idx] -= result_x;
                particles.forces[1][other_idx] -= result_y;
                force_sum[0] += result_x;
                force_sum[1] += result_y;
                if constexpr (DIMENSIONS > 2) {
                  particles.forces[2][other_idx] -= result_z;
                  force_sum[2] += result_z;
                }
              }
            }

            // Neighbour cells
            for (Neighbour<DIMENSIONS> &neighbour : cell.neighbours) {

              Cell<DIMENSIONS> &neighbour_cell = cells[neighbour.cell];
              size_t n_start = neighbour_cell.start_index;
              size_t n_end = neighbour_cell.end_index;

#pragma omp simd linear(other_idx) simdlen(8) reduction(+ : force_sum[ : DIMENSIONS]) if (vectorize)
              for (other_idx = n_start; other_idx < n_end; ++other_idx) {
                if (particles.active[other_idx] == 1) {
                  double result_x = 0;
                  double result_y = 0;
                  double result_z = 0;

                  if constexpr (DIMENSIONS == 2) {
                    force.calculateForce_2D(
                        x1, y1, mass1, type1,
                        particles.positions[0][other_idx],
                        particles.positions[1][other_idx],
                        particles.mass[other_idx],
                        particles.type[other_idx],
                        result_x,
                        result_y,
                        neighbour.correction);
                  } else {
                    force.calculateForce_3D(
                        x1, y1, z1, mass1, type1,
                        particles.positions[0][other_idx],
                        particles.positions[1][other_idx],
                        particles.positions[2][other_idx],
                        particles.mass[other_idx],
                        particles.type[other_idx],
                        result_x,
                        result_y,
                        result_z,
                        neighbour.correction);
                  }

                  particles.forces[0][other_idx] -= result_x;
                  particles.forces[1][other_idx] -= result_y;
                  force_sum[0] += result_x;
                  force_sum[1] += result_y;

                  if constexpr (DIMENSIONS > 2) {
                    particles.forces[2][other_idx] -= result_z;
                    force_sum[2] += result_z;
                  }
                }
              }
            }

            for (size_t d = 0; d < DIMENSIONS; ++d) {
              particles.forces[d][idx] += force_sum[d];
            }

          } else [[unlikely]] {
            // Stop on first inactive encountered particle
            SPDLOG_WARN("Particle crossed boundary not handled before next loop");
            //        break;
          }
        }
      }
    }
  }

  /**
   * inserts a new particle into the arena and the into the cells
   * */
  void insert(Particle<DIMENSIONS> particle) override {
    particles.insert_particle(particle);
  }

  /**
   * Returns the number of particles in the arena. O(n), because it gets counted each time
   * */
  size_t size() override {
    return particles.size;
  }

  void linear(std::function<void(Particles<DIMENSIONS> &, size_t)> function) override {
    for (size_t i = 0; i < particles.size; ++i) {
      function(particles, i);
    }
  }

  void swap_force() override {
    particles.swap_force();
  }

  /**
   * calculates the new cells of each particle and sort the Particles based on that.
   * After that all cell ranges get reconstructed.
   *
   * */
  void refresh() override {

    for (auto &cell : cells) {
      cell.start_index = 0;
      cell.end_index = 0;
    }
    for (auto &block : colored_blocks) {
      block.clear();
    }

    particles.sort([this](size_t idx) -> std::tuple<size_t, size_t, size_t> {
      std::array<double, DIMENSIONS> temp;
      for (size_t d = 0; d < DIMENSIONS; ++d) {
        temp[d] = particles.positions[d][idx];
      }
      auto cell = index.position_to_index(temp);
      return std::tuple(
          cells[cell].cell_color,
          cells[cell].block,
          cell);
    });

    size_t start = 0;
    size_t cell = particles.cell[0];

    const size_t particle_size = particles.size;
    for (size_t i = 0; i < particle_size; ++i) {

      if (particles.cell[i] >= cells.size()) {
        particles.active[i] = false;
        SPDLOG_WARN("Particle crossed boundary not handled before next loop");
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
    cells[cell].start_index = start;
    cells[cell].end_index = particle_size;

    /*
     * Rebuilds blocks
     * */
    auto r = ranges::zip_view(particles.color, particles.block, ranges::views::iota(0UL, particle_size))
        | ranges::views::chunk_by([](auto a, auto b) {
               return std::get<0>(a) == std::get<0>(b) && std::get<1>(a) == std::get<1>(b);
             });
    ranges::for_each(r, [this](auto chunk) {
      auto color = std::get<0>(chunk[0]);
      size_t start = std::get<2>(chunk[0]);
      size_t end = start + static_cast<size_t>(chunk.size());
      colored_blocks[color].emplace_back(Block(start, end));
    });

    SPDLOG_TRACE("fixed positions");
  }

 private:
  Particles<DIMENSIONS> particles;
  F force;

  /**
   * Cells with start and end offset of particles and the neighbour cells
   * */
  std::vector<Cell<DIMENSIONS>> cells;
  /**
   * an index helper
   * */
  container::index::Index<DIMENSIONS> index;

  std::vector<double> sigma;
  std::vector<double> reflecting_distance;
  /**
   * A vector of blocks for 4 different colors
   * */
  std::array<std::vector<Block>, 4> colored_blocks{};
};

}// namespace container