#pragma once

#include "Particle.h"
#include "container/linked_cell.h"
#include "utils/variants.h"
#include <ranges>
#include <vector>
namespace container {

using particle_container_variant = std::variant<Particles, LinkedCell>;

/*!
 * @brief A variant-based particle container supporting multiple underlying data structures.
 * @image html submission/worksheet3/media/Benchmark.png
 */
struct particle_container {

 public:
  explicit particle_container(particle_container_variant &&var) : var(std::move(var)) {}

  auto particles() -> Particles & {
    return std::visit(overloaded{
                          [](Particles &container) -> auto & {
                            return container;
                          },
                          [](LinkedCell &c) -> auto & {
                            return c.particles;
                          }},
                      var);
  }

  /**
  * @brief Applies a function to each particle in the container.
  *
  * @param f Function to apply to each particle.
  */
  template<typename C>
  auto linear(C f) {
    Particles &p = particles();
    for (size_t index = 0; index < p.size; ++index) {
      f(p, index);
    };
  }

  auto swap_force() {
    particles().swap_force();
  }

  /**
  * @brief Applies a function to each pair of particles in the container.
  *
  * @param f Function to apply to each pair of particles.
  */
  template<typename C>
  auto pairwise(C f) {
    std::visit(overloaded{
                   [&](LinkedCell &lc) { lc.pairwise(f); },
                   [&](Particles &p) {
                     for (size_t index = 0; index < p.size - 1; ++index) {
                       auto p1 = p.load_vectorized_single(index);
                       size_v index_vector_tmp = 0;
                       for (size_t i = 0; i < size_v::size(); ++i) {
                         index_vector_tmp[i] = i;
                       }
                       const auto index_vector = index_vector_tmp;
                       size_t i = 0;
                       while (index + 1 + (i * double_v::size()) < p.size) {
                         auto active_mask = index_vector + 1 + (i * double_v::size()) < p.size;
                         auto p2 = p.load_vectorized(index + 1 + (i * double_v::size()));
                         active_mask = active_mask && p2.active > 0;
                         auto mask = stdx::static_simd_cast<double_v>(active_mask);
                         f(p1, p2, mask);
                         p.store_force_vector(p2, index + 1 + (i * double_v::size()), mask);
                         i++;
                       }
                       p.store_force_single(p1, index);
                     }
                   }},
               var);
  }

  /**
  * @brief Returns the number of particles in the container.
  *
  * @return Size of the container.
  */
  auto size() -> size_t {
    return particles().size;
  }

  //  /**
  //  * @brief Applies boundary conditions to the container.
  //  *
  //  * @param f Function to apply for boundary conditions.
  //  */
  //  auto boundary(std::function<void(Particles &, std::tuple<size_t, size_t>)> const &f) {
  //    std::visit(
  //        [](Particles &container) {},
  //        var);
  //  }

  /**
  * @brief Inserts a particle into the container.
  *
  * @param p Particle to insert.
  */
  void insert(Particle p) {
    std::visit(overloaded{
                   [p](Particles &container) { container.insert_particle(p); },
                   [p](LinkedCell &lc) {
                     lc.insert(p);
                   }},
               var);
  }

  /**
  * @brief Updates internal data structures after position recalculations.
  */
  void refresh() {
    std::visit(overloaded{
                   [](Particles &container) {},
                   [](LinkedCell &lc) { lc.fix_positions(); }},
               var);
  }

 private:
  particle_container_variant var;
};

}// namespace container