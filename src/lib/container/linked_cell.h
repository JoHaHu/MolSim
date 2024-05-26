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

#include "range/v3/view/concat.hpp"

namespace container {

//template<typename T>
//concept boundary_condition = requires();

enum class cell_type : std::uint8_t {
  inner,
  boundary,
  halo_and_boundary,
  halo
};

class cell {
 private:
  std::vector<std::shared_ptr<Particle>> particles;

 public:
  cell_type type;
  cell() = default;
  explicit cell(
      std::vector<std::shared_ptr<Particle>> &&particles,
      cell_type type) : particles(std::move(particles)),
                        type(type){};

  auto linear() -> auto {
    return std::ranges::ref_view(particles);
  }
  constexpr auto is_halo() -> bool {
    return type == cell_type::halo_and_boundary || type == cell_type::halo;
  }

  constexpr auto is_boundary() -> bool {
    return type == cell_type::halo_and_boundary || type == cell_type::boundary;
  }
  auto insert(const std::shared_ptr<Particle> &particle) {
    particles.emplace_back(particle);
  }
  auto clear() {
    particles.clear();
  }
};

/**
 * \param dim
 * */
template<index::Index I>
class linked_cell {

 public:
  linked_cell() = delete;
  explicit linked_cell(const std::array<double, 3> &domain, double cutoff) : index(I(domain, cutoff)), cutoff(cutoff) {

    auto dim = index.dimension();

    cells.reserve((dim[0]) * (dim[1]) * (dim[2]));

    for (size_t i = 0; i < dim[0] * dim[1] * dim[2]; ++i) {
      cells.emplace_back(std::vector<std::shared_ptr<Particle>>(), cell_type::inner);
    }

    for (auto i : index) {
      auto [x, y, z] = i;
      if (x == 0 || y == 0 || z == 0 || x == dim[0] - 1 || y == dim[1] - 1 || z == dim[2] - 1) {
        cells[x + y * dim[0] + z * dim[0] * dim[1]].type = cell_type::halo_and_boundary;
      }
    }
  };

  auto halo() -> auto {
    return cells
        | std::views::filter(&cell::is_halo)
        | std::views::transform(&cell::linear)
        | std::views::join;
  }

  auto boundary() -> auto {
    return cells
        | std::views::filter(&cell::is_boundary)
        | std::views::transform(&cell::linear)
        | std::views::join;
  }

  auto linear() -> auto {
    return std::ranges::ref_view(store);
  }

  auto pairwise() -> auto {

    return index
        | std::views::transform([this](std::tuple<size_t, size_t, size_t> idx) {
             auto [x, y, z] = idx;
             auto cell_idx = index.dimension_to_index({x, y, z});
             // TODO make this vector a member of cell and initialize when creating the cell,
             //  because the position of a cell and therefore it's neighbours do not change after initialization
             auto cartesian_products = std::vector<std::ranges::cartesian_product_view<
                 std::ranges::ref_view<std::vector<std::shared_ptr<Particle>>>,
                 std::ranges::ref_view<std::vector<std::shared_ptr<Particle>>>>>();

             cell &cell = cells[cell_idx];
             if (cell.type == cell_type::inner) {
              for (int x = -1; x <= 1; ++x) {
                for (int y = -1; y <= 1; ++y) {
                  auto prod = std::views::cartesian_product(cell.linear(), cells[index.offset(cell_idx, {x, y, 1})].linear());
                  cartesian_products.emplace_back(prod);
                }
              }
               cartesian_products.emplace_back(cell.linear(), cells[index.offset(cell_idx, {1, -1, 0})].linear());
               cartesian_products.emplace_back(cell.linear(), cells[index.offset(cell_idx, {1, 0, 0})].linear());
               cartesian_products.emplace_back(cell.linear(), cells[index.offset(cell_idx, {1, 1, 0})].linear());
               cartesian_products.emplace_back(cell.linear(), cells[index.offset(cell_idx, {0, 1, 0})].linear());
             }

             auto joined = std::move(cartesian_products) | std::views::join;
             return ranges::concat_view(cell.linear() | combination, std::move(joined));
           })
        | std::views::join;
  }

  auto insert(Particle &particle) {
    auto shared = std::make_shared<Particle>(particle);
    store.emplace_back(shared);
    insert_shared(shared);
  }

  auto size() {
    return linear().size();
  }

  // TODO better fixup than completly replace
  auto fix_positions() {
    std::ranges::for_each(cells, &cell::clear);
    std::ranges::for_each(linear(), [this](std::shared_ptr<Particle> &p) {
      insert_shared(p);
    });
  }

 private:
  auto insert_shared(std::shared_ptr<Particle> shared) {

    auto idx = index.position_to_index(shared->position);
    //    if (idx[0] < index.dimension()[0] && idx[1] < index.dimension()[1] && idx[2] < index.dimension()[2]) {
    // TODO proper handling of boundary conditions
    if (idx < index.max_index()){
      cells[idx].insert(shared);
    }
    //    } else {
    //      spdlog::warn("a particle has positions that is out of bounds {} {} {}", shared->position[0], shared->position[1], shared->position[2]);
    //    }
  }

  // TODO use an arena allocator for all shared_ptr so that they will be sequential in memory
  std::vector<std::shared_ptr<Particle>> store{};
  std::vector<cell> cells;
  I index;
  double cutoff;
};

}// namespace container
