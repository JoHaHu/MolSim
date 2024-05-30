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
                        type(type) {};

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

template<index::Index I>
class linked_cell_functor;

template<index::Index I>
class linked_cell {
  friend class linked_cell_functor<I>;

 public:
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

    pairwise_functor = std::optional<linked_cell_functor<I>>();
  };

  static auto make_linked_cell(const std::array<double, 3> &domain, double cutoff) -> std::shared_ptr<linked_cell<I>> {
    auto ptr = std::make_shared<linked_cell<I>>(domain, cutoff);
    auto range = linked_cell_functor<I>(ptr);
    ptr->pairwise_functor = std::optional(range);
    return ptr;
  }

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
    if (!pairwise_view.has_value()) {
      pairwise_view = std::optional(
          index
          | std::views::transform(*pairwise_functor)
          | std::views::join);
    }

    return *pairwise_view;
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
  linked_cell() = delete;

  auto insert_shared(std::shared_ptr<Particle> shared) {

    auto idx = index.position_to_index(shared->position);
    //    if (idx[0] < index.dimension()[0] && idx[1] < index.dimension()[1] && idx[2] < index.dimension()[2]) {
    // TODO proper handling of boundary conditions
    if (idx < index.max_index()) {
      cells[idx].insert(shared);
    } else {
      spdlog::warn("a particle has positions that is out of bounds {} {} {}", shared->position[0], shared->position[1], shared->position[2]);
    }
  }

  // TODO use an arena allocator for all shared_ptr so that they will be sequential in memory
  std::vector<std::shared_ptr<Particle>> store{};
  std::vector<cell> cells;
  I index;
  double cutoff;
  std::optional<linked_cell_functor<I>> pairwise_functor;
  std::optional<std::ranges::join_view<std::ranges::transform_view<I, linked_cell_functor<I>>>> pairwise_view{};
};

template<index::Index I>
class linked_cell_functor {
 public:
  using particle_vector = std::ranges::ref_view<std::vector<std::shared_ptr<Particle>>>;
  using product_range = std::vector<std::ranges::cartesian_product_view<
      particle_vector,
      particle_vector>>;

  explicit linked_cell_functor(std::shared_ptr<linked_cell<I>> lc) : lc(lc) {}

  auto operator()(std::tuple<size_t, size_t, size_t> idx) -> ranges::concat_view<container::combination_view<particle_vector>, std::ranges::join_view<std::ranges::owning_view<product_range>>> {
    auto [x, y, z] = idx;
    auto cell_idx = lc->index.dimension_to_index({x, y, z});
    auto cartesian_products = product_range();

    cell &cell = lc->cells[cell_idx];
    if (cell.type == cell_type::inner) {
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
    return ranges::concat_view(cell.linear() | combination, std::move(joined));
  }

 private:
  std::shared_ptr<linked_cell<I>> lc;
};

}// namespace container
