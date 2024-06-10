#pragma once

#include "Particle.h"
#include "index.h"
#include "utils/variants.h"
#include <ranges>
#include <vector>
namespace container {

using particle_container_variant = std::variant<Particles>;

/*!
 * @brief A variant-based particle container supporting multiple underlying data structures.
 * @image html submission/worksheet3/media/Benchmark.png
 */
struct particle_container {

 public:
  explicit particle_container(particle_container_variant &&var) : var(std::move(var)) {}

  /**
  * @brief Applies a function to each particle in the container.
  *
  * @param f Function to apply to each particle.
  */
  auto linear(std::function<void(Particles &, size_t)> const &f) {
    std::visit(
        [&f](Particles &container) {
          std::ranges::for_each(container.linear(), [&f, &container](size_t index) { f(container, index); });
        },
        var);
  }

  auto swap_force() {
    std::visit([](Particles &p) {
      p.swap_force();
    },
               var);
  }

  /**
  * @brief Applies a function to each pair of particles in the container.
  *
  * @param f Function to apply to each pair of particles.
  */
  auto pairwise(std::function<void(Particles &p, std::tuple<size_t, size_t>)> const &f) {
    std::visit(
        [&](Particles &container) { std::ranges::for_each(container.pairwise(), [&](auto index) { f(container, index); }); },
        var);
  }

  /**
  * @brief Returns the number of particles in the container.
  *
  * @return Size of the container.
  */
  auto size() -> size_t {
    return std::visit(
        [](auto &c) { return c.size; },
        var);
  }

  /**
  * @brief Applies boundary conditions to the container.
  *
  * @param f Function to apply for boundary conditions.
  */
  auto boundary(std::function<void(Particles &, std::tuple<size_t, size_t>)> const &f) {
    std::visit(
        [](Particles &container) {},
        var);
  }

  /**
  * @brief Inserts a particle into the container.
  *
  * @param p Particle to insert.
  */
  void insert(Particle p) {
    std::visit(
        [p](Particles &container) { container.insert_particle(p); },
        var);
  }

  /**
  * @brief Updates internal data structures after position recalculations.
  */
  void refresh() {
    std::visit(
        [](Particles &container) {},
        var);
  }

 private:
  particle_container_variant var;
};

}// namespace container