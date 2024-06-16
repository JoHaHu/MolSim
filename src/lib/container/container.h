#pragma once

#include "Particle.h"
#include "container/linked_cell.h"
#include "utils/variants.h"
#include <ranges>
#include <vector>
namespace container {

template<const size_t DIMENSIONS>
using ParticleContainerVariant = std::variant<Particles<DIMENSIONS>,
                                              LinkedCell<DIMENSIONS>>;

/*!
 * @brief A variant-based particle container supporting multiple underlying data structures.
 * @image html submission/worksheet3/media/Benchmark.png
 */

template<const size_t DIMENSIONS>
struct ParticleContainer {

 public:
  explicit ParticleContainer(ParticleContainerVariant<DIMENSIONS> &&var) : var(std::move(var)) {
  }

  auto particles() -> Particles<DIMENSIONS> & {
    return std::visit(overloaded{
                          [](Particles<DIMENSIONS> &container) -> auto & {
                            return container;
                          },
                          [](LinkedCell<DIMENSIONS> &c) -> auto & {
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
    Particles<DIMENSIONS> &p = particles();
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
                   [&](LinkedCell<DIMENSIONS> &lc) { lc.pairwise(f); },
                   [&](Particles<DIMENSIONS> &p) {
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
                         auto mask = stdx::static_simd_cast<double_v>(active_mask) && p2.active;
                         std::array<double_v, DIMENSIONS> correction;
                         for (int j = 0; j < DIMENSIONS; ++j) {
                           correction[i] = double_v(0.0);
                         }
                         f(p1, p2, mask, correction);
                         p.store_force_vector(p2, index + 1 + (i * double_v::size()));
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

  /**
  * @brief Inserts a particle into the container.
  *
  * @param p Particle to insert.
  */
  void insert(Particle<DIMENSIONS> p) {
    std::visit(overloaded{
                   [p](Particles<DIMENSIONS> &container) { container.insert_particle(p); },
                   [p](LinkedCell<DIMENSIONS> &lc) {
                     lc.insert(p);
                   }},
               var);
  }
  /**
   * @brief applies boundary conditions, only effective when using linked cells
   * */
  template<typename Callable>
  void boundary(Callable f) {
    std::visit(overloaded{
                   [](Particles<DIMENSIONS> &container) {},
                   [&f](LinkedCell<DIMENSIONS> &lc) { lc.boundary(f); }},
               var);
  }

  /**
  * @brief Updates internal data structures after position recalculations.
  */
  void refresh() {
    std::visit(overloaded{
                   [](Particles<DIMENSIONS> &container) {},
                   [](LinkedCell<DIMENSIONS> &lc) { lc.fix_positions(); }},
               var);
  }

 private:
  ParticleContainerVariant<DIMENSIONS> var;
};

}// namespace container