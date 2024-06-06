#pragma once

#include "Particle.h"
#include "index.h"
#include "linked_cell.h"
#include "utils/variants.h"
#include <ranges>
#include <vector>
namespace container {

 /**
 * @brief A variant-based particle container supporting multiple underlying data structures.
 */
using particle_container_variant = std::variant<std::vector<Particle>,
                                                linked_cell<index::row_major_index>,
                                                linked_cell<index::half_index>,
                                                linked_cell<index::morton_index>,
                                                linked_cell<index::power_index>>;

struct particle_container {
 public:
  explicit particle_container(particle_container_variant &&var) : var(std::move(var)) {}


 /**
  * @brief Applies a function to each particle in the container.
  *
  * @param f Function to apply to each particle.
  */
  auto linear(std::function<void(Particle &)> const &f) {
    std::visit(overloaded{
                   [&f](std::vector<Particle> &container) { std::ranges::for_each(container, f); },
                   [&f](auto &container) { std::ranges::for_each(container.linear(), f); },
               },
               var);
  }

 /**
  * @brief Applies a function to each pair of particles in the container.
  *
  * @param f Function to apply to each pair of particles.
  */
  auto pairwise(std::function<void(std::tuple<Particle &, Particle &>)> const &f) {
    std::visit(overloaded{
                   [&f](std::vector<Particle> &container) { std::ranges::for_each(container | combination, f); },
                   [&f](auto &container) { std::ranges::for_each(container.pairwise(), f); }},
               var);
  }

 /**
  * @brief Returns the number of particles in the container.
  *
  * @return Size of the container.
  */
  auto size() -> size_t {
    return std::visit(overloaded{
                          [](std::vector<Particle> &c) { return c.size(); },
                          [](auto &c) { return c.size(); },
                      },
                      var);
  }

 /**
  * @brief Applies boundary conditions to the container.
  *
  * @param f Function to apply for boundary conditions.
  */
  auto boundary(std::function<void(std::tuple<Particle &, Particle &>)> const &f) {
    std::visit(overloaded{
                   [](std::vector<Particle> &container) {},
                   [&f](auto &container) {
                     calculate_boundary_condition(container, f);
                   },
               },
               var);
  }

 /**
  * @brief Inserts a particle into the container.
  *
  * @param p Particle to insert.
  */
  void insert(Particle &&p) {
    std::visit(overloaded{
                   [&p](std::vector<Particle> &container) { container.emplace_back(p); },
                   [&p](auto &container) { container.insert(std::move(p)); }},
               var);
  }

 /**
  * @brief Updates internal data structures after position recalculations.
  */
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