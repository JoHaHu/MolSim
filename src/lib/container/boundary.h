#pragma once

#include "container.h"
#include "linked_cell.h"
#include <variant>

namespace container::boundary {

static void outflow(auto &lc, arena<Particle>::entry &entry, orientation o) {

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
  entry.active = !condition;
}

static void reflecting(auto &lc, auto &p, orientation o, auto fc) {
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

template<index::Index I>
static void calculate_boundary_condition(linked_cell<I> &lc,
                                         std::function<void(std::tuple<Particle &, Particle &>)> const &force_calculation

) {

  std::ranges::for_each(lc.boundary(), [&lc, &force_calculation](cell &cell) {
    for (auto [side, b] : std::views::enumerate(cell.boundary)) {
      auto o = orientation(side);
      switch (b) {
        case boundary_condition::outflow:
          std::ranges::for_each(cell.particles, [&lc, &o](auto &e) { outflow(lc, e, o); });
          break;
        case boundary_condition::reflecting:
          std::ranges::for_each(cell.linear(), [&lc, &o, &force_calculation](auto &p) { reflecting(lc, p, o, force_calculation); });
          break;
        case boundary_condition::none: break;
      }
    }
  });
}

}// namespace container::boundary