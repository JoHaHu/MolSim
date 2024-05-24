#pragma once

#include "Particle.h"
#include "combination.h"
#include "container.h"
#include "index.h"
#include <math.h>

#include <array>
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
 public:
  cell() = default;
  explicit cell(std::vector<Particle> &&particles, cell_type type) : particles(std::move(particles)), type(type){};

  auto linear() -> auto {
    return std::ranges::ref_view(particles);
  }

  auto pairwise() -> auto {
    return particles | combination;
  }

  constexpr auto halo() -> bool {
    return type == cell_type::halo_and_boundary || type == cell_type::halo;
  }

  constexpr auto boundary() -> bool {
    return type == cell_type::halo_and_boundary || type == cell_type::boundary;
  }

 private:
  std::vector<Particle> particles;
  cell_type type;
};

/**
 * \param dim
 * */
template<index::Index I>
class linked_cells {

 public:
  linked_cells() = delete;
  explicit linked_cells(const std::array<double, 3> &domain, const double cutoff) : domain(domain) {

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

    cells.reserve((dim[0]) * (dim[1]) * (dim[2]));
  };

  auto halo() -> auto {
    return cells
        | std::views::filter(&cell::halo)
        | std::views::transform(&cell::linear)
        | std::views::join;
  }

  auto boundary() -> auto {
    return cells
        | std::views::filter(&cell::boundary)
        | std::views::transform(&cell::linear)
        | std::views::join;
  }

  auto linear() -> auto {
    return cells
        | std::views::transform(&cell::linear)
        | std::views::join;
  }

  auto pairwise() -> auto {
    return cells
        | std::views::transform(&cell::pairwise)
        | std::views::join;
  }

 private:
  std::array<double, 3> domain{};
  std::array<size_t, 3> dim{};
  std::vector<cell> cells;
  I index = I(dim);
};

}// namespace container
