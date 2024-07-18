#pragma once

#include "Container.h"
#include "simulator/physics/Force.h"
#include <cmath>
#include <memory>

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

  void pairwise(bool parallel, bool vector) override {
    const auto particle_size = particles.size;
    for (size_t index_1 = 0; index_1 < particle_size - 1; ++index_1) {
      if (particles.active[index_1]) {
        auto particle_force = std::array<double, DIMENSIONS>();

#pragma omp simd reduction(+ \
                           : particle_force[:DIMENSIONS])
        for (size_t index_2 = index_1 + 1; index_2 < particle_size; ++index_2) {
          if (particles.active[index_2]) {
            std::array<double, DIMENSIONS> result = {0.0};
            std::array<double, 2> correction_2d = {0.0, 0.0};
            std::array<double, 3> correction_3d = {0.0, 0.0, 0.0};

            if (DIMENSIONS == 2) {
              force->calculateForce_2D(
                  particles.positions[0][index_1], particles.positions[1][index_1],
                  particles.mass[index_1], particles.type[index_1],
                  particles.positions[0][index_2], particles.positions[1][index_2],
                  particles.mass[index_2], particles.type[index_2],
                  result[0], result[1], correction_2d);
            } else if (DIMENSIONS == 3) {
              force->calculateForce_3D(
                  particles.positions[0][index_1], particles.positions[1][index_1], particles.positions[2][index_1],
                  particles.mass[index_1], particles.type[index_1],
                  particles.positions[0][index_2], particles.positions[1][index_2], particles.positions[2][index_2],
                  particles.mass[index_2], particles.type[index_2],
                  result[0], result[1], result[2], correction_3d);
            }

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
