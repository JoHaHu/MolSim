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
      std::ranges::ref_view<std::vector<Particle>>,
      decltype(std::declval<container::linked_cell<index::simple_index>>().linear())>;

  using pairwise_variant = std::variant<
      container::combination_view<std::ranges::ref_view<std::vector<Particle>>>,
      decltype(std::declval<container::linked_cell<index::simple_index>>().pairwise())>;
  using halo_variant = std::variant<
      container::combination_view<std::ranges::ref_view<std::vector<Particle>>>,
      decltype(std::declval<container::linked_cell<index::simple_index>>().halo())>;
  using boundary_variant = std::variant<
      std::ranges::empty_view<std::reference_wrapper<Particle>>,
      decltype(std::declval<container::linked_cell<index::simple_index>>().boundary())>;
  using ghost_variant = std::variant<
      std::ranges::empty_view<std::tuple<std::reference_wrapper<Particle>, Particle>>,
      decltype(std::declval<container::linked_cell<index::simple_index>>().ghosts())>;

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

  auto boundary() -> boundary_variant {
    return std::visit<boundary_variant>(overloaded{
                                            [](std::vector<Particle> &container) { return std::ranges::empty_view<std::reference_wrapper<Particle>>(); },
                                            [](std::shared_ptr<linked_cell<index::simple_index>> &container) { return container->boundary(); },
                                        },
                                        var);
  }
  auto ghosts() -> ghost_variant {
    return std::visit<ghost_variant>(overloaded{
                                         [](std::vector<Particle> &container) { return std::ranges::empty_view<std::tuple<std::reference_wrapper<Particle>, Particle>>(); },
                                         [](std::shared_ptr<linked_cell<index::simple_index>> &container) { return container->ghosts(); },
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