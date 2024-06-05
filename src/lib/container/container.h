#pragma once

#include "Particle.h"
#include "boundary.h"
#include "index.h"
#include "linked_cell.h"
#include "utils/variants.h"
#include <ranges>
#include <vector>
namespace container {

using particle_container_variant = std::variant<std::vector<Particle>,
                                                linked_cell<index::row_major_index>,
                                                linked_cell<index::half_index>,
                                                linked_cell<index::morton_index>,
                                                linked_cell<index::power_index>>;

struct particle_container {
 public:
  explicit particle_container(particle_container_variant &&var) : var(std::move(var)) {}

  /**
   * Applies a function to particle in the container
   * */
  auto linear(std::function<void(Particle &)> const &f) {
    std::visit(overloaded{
                   [&f](std::vector<Particle> &container) { std::ranges::for_each(container, f); },
                   [&f](auto &container) { std::ranges::for_each(container.linear(), f); },
               },
               var);
  }

  /**
   * Applies a function to pairs in the container
   * */
  auto pairwise(std::function<void(std::tuple<Particle &, Particle &>)> const &f) {
    std::visit(overloaded{
                   [&f](std::vector<Particle> &container) { std::ranges::for_each(container | combination, f); },
                   [&f](auto &container) { std::ranges::for_each(container.pairwise(), f); }},
               var);
  }

  /**
 * Return the size of the container
 * */
  auto size() -> size_t {
    return std::visit(overloaded{
                          [](std::vector<Particle> &c) { return c.size(); },
                          [](auto &c) { return c.size(); },
                      },
                      var);
  }

  /**
   * applies boundary conditions
   * */
  auto boundary(std::function<void(std::tuple<Particle &, Particle &>)> const &f) {
    std::visit(overloaded{
                   [](std::vector<Particle> &container) {},
                   [&f](auto &container) {
                     boundary::calculate_boundary_condition(container, f);
                   },
               },
               var);
  }

  /**
   * Inserts boundary conditions into the container
   * */
  void insert(Particle &p) {
    std::visit(overloaded{
                   [&p](std::vector<Particle> &container) { container.emplace_back(p); },
                   [&p](auto &container) { container.insert(p); }},
               var);
  }

  /**
   * updates internal datastructures. Should be called after every recalculation of positions
   * */
  void refresh() {
    std::visit(overloaded{
                   [](std::vector<Particle> &container) {},
                   [](auto &container) { container.fix_positions(); }},
               var);
  }

 private:
  particle_container_variant var;
};

}// namespace container