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
#include "spdlog/spdlog.h"

namespace container {

enum class boundary_condition : std::uint8_t {
  outflow,
  reflecting,
  //  periodic,
  none
};

enum class orientation : std::uint8_t {
  front,
  back,
  left,
  right,
  bottom,
  top
};

enum class cell_type : std::uint8_t {
  inner,
  boundary,
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
  explicit cell(std::vector<std::reference_wrapper<arena<Particle>::entry>> &&particles,
                cell_type type,
                std::array<size_t, 3> idx,
                std::array<double, 3> widths) : particles(std::move(particles)), type(type), idx(idx), widths(widths) {};

  auto linear() -> auto {
    return particles
        | std::views::transform(&container::arena<Particle>::entry::data);
  }
  auto pairwise() -> auto & {
    return *range;
  };

  constexpr auto is_boundary() -> bool {
    return type == cell_type::boundary;
  }

  auto insert(arena<Particle>::entry &particle) {
    particles.emplace_back(particle);
  }

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

template<index::Index I>
class linked_cell {
  friend class cell;

 public:
  linked_cell() = delete;
  explicit linked_cell(const std::array<double, 3> &domain, double cutoff, std::array<boundary_condition, 6> bc, unsigned int number_of_particles, double sigma)
      : arena(number_of_particles),
        bc(bc),
        index(I(domain, cutoff)),
        sigma(sigma) {
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

  static auto create_range(linked_cell &lc, std::array<size_t, 3> idx) -> auto {
    auto cell_idx = lc.index.dimension_to_index(idx);

    cell &cell = lc.cells[cell_idx];

    auto [radius_x, radius_y, radius_z] = lc.index.radius;

    auto cartesian_products = cell::product_range();

    for (long x = -(long) radius_x; x <= (long) radius_x; ++x) {
      for (long y = -(long) radius_y; y <= (long) radius_y; ++y) {
        for (long z = -(long) radius_z; z <= (long) radius_z; ++z) {
          if (x != 0 || y != 0 || z != 0) {
            if (z > 0 || x > 0 || (x >= 0 && z >= 0 && y > 0)) {
              auto other_index = lc.index.offset(cell_idx, {x, y, z});
              // checks that cell is not out of bounds and not the same cell
              if (other_index < lc.cells.size() && other_index != cell_idx) {
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

  auto boundary() -> auto {
    return cells
        | std::views::filter(&cell::is_boundary);
  }

  auto linear() -> auto {
    return arena.range_entries()
        | std::views::transform(&container::arena<Particle>::entry::data);
  }

  auto pairwise() -> auto {

    return cells
        | std::views::transform(&cell::pairwise)
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
    SPDLOG_TRACE("fixed positions");
  }

 private:
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
};

}// namespace container
