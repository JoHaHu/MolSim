#pragma once

#include "Particle.h"
#include "boundary.h"
#include "index.h"
#include "linked_cell.h"
#include "utils/variants.h"
#include <ranges>
#include <vector>
namespace container {

using particle_container_variant = std::variant<std::vector<Particle>, linked_cell<index::simple_index>>;

struct particle_container {
 public:
  explicit particle_container(particle_container_variant &&var) : var(std::move(var)) {}

  auto linear(std::function<void(Particle &)> const &f) {
    std::visit(overloaded{
                   [&f](std::vector<Particle> &container) { std::ranges::for_each(container, f); },
                   [&f](linked_cell<index::simple_index> &container) { std::ranges::for_each(container.linear(), f); },
               },
               var);
  }

  auto pairwise(std::function<void(std::tuple<Particle &, Particle &>)> const &f) {
    std::visit(overloaded{
                   [&f](std::vector<Particle> &container) { std::ranges::for_each(container | combination, f); },
                   [&f](linked_cell<index::simple_index> &container) { std::ranges::for_each(container.pairwise(), f); }},
               var);
  }

  auto size() -> size_t {
    return std::visit(overloaded{
                          [](std::vector<Particle> &c) { return c.size(); },
                          [](linked_cell<index::simple_index> &c) { return c.size(); },
                      },
                      var);
  }

  auto boundary(std::function<void(std::tuple<Particle &, Particle &>)> const &f) {
    std::visit(overloaded{
                   [](std::vector<Particle> &container) {},
                   [&f](linked_cell<index::simple_index> &container) {
                     boundary::calculate_boundary_condition(container, f);
                   },
               },
               var);
  }

  void insert(Particle &p) {
    std::visit(overloaded{
                   [&p](std::vector<Particle> &container) { container.emplace_back(p); },
                   [&p](linked_cell<index::simple_index> &container) { container.insert(p); }},
               var);
  }

  void refresh() {
    std::visit(overloaded{
                   [](std::vector<Particle> &container) {},
                   [](linked_cell<index::simple_index> &container) { container.fix_positions(); }},
               var);
  }

 private:
  particle_container_variant var;
};

}// namespace container