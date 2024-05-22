#pragma once

#include "LinkedCell.h"
#include "Particle.h"
#include "utils/variants.h"
#include <ranges>
#include <vector>
namespace container {

template<typename T>
concept container =
    std::ranges::forward_range<T>
    && std::ranges::sized_range<T>
    && std::same_as<typename std::ranges::range_value_t<T>, Particle>;

template<container... PC>
using particle_container_variant = std::variant<PC...>;

struct particle_container : public std::ranges::view_interface<particle_container> {
 public:
  using variant = particle_container_variant<std::vector<Particle>, container::LinkedCells<index::SimpleIndex>>;
  using pair_range = std::variant<combination_view<std::ranges::ref_view<std::vector<Particle>>>, LinkedCells<index::SimpleIndex>::view>;

  explicit particle_container(const variant &&var) : var(var) {}

  auto begin() -> auto {
    return std::visit([](auto &container) { return container.begin(); }, var);
  }
  auto end() -> auto {
    return std::visit([](auto &container) { return container.end(); }, var);
  }

  auto pairwise() -> pair_range {
    return std::visit(overloads{
                          [](std::vector<Particle> &container) { return pair_range(container | combination); },
                          [](LinkedCells<index::SimpleIndex> &container) { return pair_range({container.pairwise()}); },
                      },
                      var);
  }

 private:
  variant var;
};

static_assert(container<particle_container>);
}// namespace container