#pragma once

#include "Particle.h"
#include "linked_cell.h"
#include "utils/variants.h"
#include <ranges>
#include <vector>
namespace container {

using particle_container_variant = std::variant<std::vector<Particle>, std::shared_ptr<container::linked_cell<index::simple_index>>>;

struct particle_container {
 public:
  using linear_variant = std::variant<
      std::ranges::ref_view<std::vector<Particle>>>;

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
    return std::visit<linear_variant>(overloaded{
                                          [](std::vector<Particle> &container) { return std::ranges::ref_view(container); },
                                          [](std::shared_ptr<linked_cell<index::simple_index>> &container) { return container->linear(); },
                                      },
                                      var);
  }

  auto pairwise() -> pairwise_variant {
    return std::visit<pairwise_variant>(overloaded{
                                            [](std::vector<Particle> &container) { return container | combination; },
                                            [](std::shared_ptr<linked_cell<index::simple_index>> &container) { return container->pairwise(); }},
                                        var);
  }

  auto size() -> size_t {
    return std::visit(overloaded{
                          [](std::vector<Particle> &c) { return c.size(); },
                          [](std::shared_ptr<linked_cell<index::simple_index>> &c) { return c->size(); },
                      },
                      var);
  }

  void insert(Particle &p) {
    std::visit(overloaded{
                   [&p](std::vector<Particle> &container) { container.emplace_back(p); },
                   [&p](std::shared_ptr<linked_cell<index::simple_index>> &container) { container->insert(p); }},
               var);
  }

  void refresh() {
    std::visit(overloaded{
                   [](std::vector<Particle> &container) {},
                   [](std::shared_ptr<linked_cell<index::simple_index>> &container) { container->fix_positions(); }},
               var);
  }

 private:
  particle_container_variant var;
};

}// namespace container