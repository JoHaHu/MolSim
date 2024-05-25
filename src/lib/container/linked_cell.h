#pragma once

#include "Particle.h"
#include "combination.h"
#include "container.h"
#include "index.h"
#include "utils/ArrayUtils.h"

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
  std::vector<Particle> particles;

 public:
  cell_type type;
  cell() = default;
  explicit cell(
      std::vector<Particle> &&particles,
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
  auto insert(Particle &particle) {
    particles.emplace_back(particle);
  }
};

/**
 * \param dim
 * */
template<index::Index I>
class linked_cell {

 public:
  linked_cell() = delete;
  explicit linked_cell(const std::array<double, 3> &domain, double cutoff) : domain(domain), cutoff(cutoff) {
    {
      double x_rem = NAN;
      double x = NAN;
      x_rem = std::modf(domain[0] / cutoff, &x);
      x = x + 1;
      if (x_rem == 0.0) {
        x = x + 1;
      }
      double y_rem = NAN;
      double y = NAN;
      y_rem = std::modf(domain[1] / cutoff, &y);
      y = y + 1;
      if (y_rem == 0.0) {
        y = y + 1;
      }

      double z_rem = NAN;
      double z = NAN;
      z_rem = std::modf(domain[2] / cutoff, &z);
      z = z + 1;
      if (z_rem == 0.0) {
        z = z + 1;
      }

      dim = {(size_t) x,
             (size_t) y,
             (size_t) z};
    }
    cells.reserve((dim[0]) * (dim[1]) * (dim[2]));

    for (size_t i = 0; i < dim[0] * dim[1] * dim[2]; ++i) {
      cells.emplace_back(std::vector<Particle>(), cell_type::inner);
    }
    index = I(dim);
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
    return cells
        | std::views::transform(&cell::linear)
        | std::views::join;
  }

  auto pairwise() -> auto {

    return index
        | std::views::transform([this](std::tuple<size_t, size_t, size_t> idx) {
             auto [x, y, z] = idx;
             auto cell_idx = x + dim[0] * y + z * dim[0] * dim[1];
             // TODO make this vector a member of cell and initialize when creating the cell,
             //  because the position of a cell and therefore it's neighbours do not change
             auto cartesian_products = std::vector<std::ranges::cartesian_product_view<
                 std::ranges::ref_view<std::vector<Particle>>,
                 std::ranges::ref_view<std::vector<Particle>>>>();

             cell &cell = cells.at(cell_idx);
             if (cell.type == cell_type::inner) {
               for (int x = -1; x <= 1; ++x) {
                 for (int y = -1; y <= 1; ++y) {
                   auto prod = std::views::cartesian_product(cell.linear(), cells.at(cell_idx + x + y * dim[1] + dim[1] * dim[2]).linear());
                   cartesian_products.emplace_back(prod);
                 }
               }
             }
             auto joined = std::move(cartesian_products) | std::views::join;
             return ranges::concat_view(cell.linear() | combination, std::move(joined));
           })
        | std::views::join;
  }

  auto insert(Particle &particle) {
    auto [x, y, z] = particle.position;
    size_t x_idx = std::ceil(x / cutoff);
    size_t y_idx = std::ceil(y / cutoff);
    size_t z_idx = std::ceil(z / cutoff);
    cells[x_idx + y_idx * dim[0] + z_idx * dim[0] * dim[1]].insert(particle);
  }

 private:
  std::array<double, 3> domain{};
  std::array<size_t, 3> dim{};
  std::vector<cell> cells;
  I index = I(dim);
  double cutoff;
};

}// namespace container
