#pragma once

#include "Container.h"
#include "simulator/physics/Force.h"
#include <cmath>

template<const size_t DIMENSIONS>
struct ParticleContainer final : public container::Container<DIMENSIONS> {
 private:
  Particles<DIMENSIONS> particles;
  std::unique_ptr<simulator::physics::Force> force;

 public:
  ~ParticleContainer() override = default;

  void linear(std::function<void(Particles<DIMENSIONS> &, size_t)> function) override {
    for (size_t i = 0; i < particles.size; ++i) {
      function(particles, i);
    }
  }
  void swap_force() override {
    particles.swap_force();
  }
  size_t size() override {
    return particles.size;
  }
  void insert(Particle<DIMENSIONS> p) override {
    particles.insert_particle(p);
  }
  void boundary() override {}
  void refresh() override {}

  auto pairwise(bool parallel, bool vector) {

    const auto particle_size = particles.size;
    for (size_t index_1 = 0; index_1 < particle_size - 1; ++index_1) {
      if (particles.active[index_1]) {
        auto particle_force = std::array<double, DIMENSIONS>();

#pragma omp simd reduction(+ : particle_force[ : DIMENSIONS])
        for (size_t index_2 = index_1 + 1; index_2 < particle_size; ++index_2) {
          if (particles.active[index_1]) {

            std::array<double, DIMENSIONS> result;
            force->calculateForce_2D(particles, index_1, index_2, result, {0, 0, 0});

            for (size_t d = 0; d < DIMENSIONS; ++d) {
              particles.forces[d][index_2] -= result[d];
              particle_force[d] += result[d];
            }
          }
        }

        for (size_t d = 0; d < DIMENSIONS; ++d) {
          particles.forces[d][index_1] += particle_force[d];
        }
      }
    }
  }
};
