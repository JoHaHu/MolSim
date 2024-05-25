#pragma once

#include "Particle.h"
#include "linked_cell.h"
#include "utils/variants.h"
#include <ranges>
#include <vector>
namespace container {

using particle_container_variant = std::variant<std::vector<Particle>, container::linked_cell<index::simple_index>>;

struct particle_container {
 public:
  using linear_variant = std::variant<
      std::ranges::ref_view<std::vector<Particle>>,
      decltype(std::declval<container::linked_cell<index::simple_index>>().linear())>;

  using pairwise_variant = std::variant<
      container::combination_view<std::ranges::ref_view<std::vector<Particle>>>,
      decltype(std::declval<container::linked_cell<index::simple_index>>().pairwise())>;
  using halo_variant = std::variant<
      container::combination_view<std::ranges::ref_view<std::vector<Particle>>>,
      decltype(std::declval<container::linked_cell<index::simple_index>>().halo())>;
  using boundary_variant = std::variant<
      container::combination_view<std::ranges::ref_view<std::vector<Particle>>>,
      decltype(std::declval<container::linked_cell<index::simple_index>>().boundary())>;

  explicit particle_container(const particle_container_variant &&var) : var(var) {}

  auto linear() -> linear_variant {
    return std::visit<linear_variant>(overloads{
                                          [](std::vector<Particle> &container) { return std::ranges::ref_view(container); },
                                          [](linked_cell<index::simple_index> &container) { return container.linear(); },
                                      },
                                      var);
  }

  auto pairwise() -> pairwise_variant {
    return std::visit<pairwise_variant>(overloads{
                                            [](std::vector<Particle> &container) { return container | combination; },
                                            [](linked_cell<index::simple_index> &container) { return container.pairwise(); }},
                                        var);
  }

  auto size() -> int {
    return 0;
  }

  void insert(Particle &&p) {
  }

 private:
  particle_container_variant var;
};

}// namespace container