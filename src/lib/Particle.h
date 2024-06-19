#pragma once

#include "utils/types.h"
#include <array>
#include <execution>
#include <experimental/simd>
#include <ranges>
#include <vector>

/**
 * the Particle class representing the simulated particles
 * */
template<const size_t DIMENSIONS>
class Particle {
 public:
  /**
   * Position of the particle
   */
  std::array<double, DIMENSIONS> position{};

  /**
   * Velocity of the particle
   */
  std::array<double, DIMENSIONS> velocity{};

  /**
   * Force effective on this particle
   */
  std::array<double, DIMENSIONS> force{};

  /**
   * Force which was effective on this particle
   */
  std::array<double, DIMENSIONS> old_force{};

  /**
   * Mass of this particle
   */
  double mass{};

  /**
   * Type of the particle. Use it for whatever you want (e.g. to separate
   * molecules belonging to different bodies, matters, and so on)
   */
  int type;

  Particle(
      // for visualization, we need always 3 coordinates
      // -> in case of 2d, we use only the first and the second
      std::array<double, DIMENSIONS> x_arg, std::array<double, DIMENSIONS> v_arg, double m_arg,
      int type = 0) : position(x_arg), velocity(v_arg), mass(m_arg), type(type)

  {
    for (int i = 0; i < DIMENSIONS; ++i) {
      force[i] = 0.0;
      old_force[i] = 0.0;
    }
  }
};

template<const size_t DIMENSIONS>
struct VectorizedParticle {
  std::array<double_v, DIMENSIONS> position{};
  std::array<double_v, DIMENSIONS> force{};
  double_v mass{};
  long_v type;
  double_mask active;

  VectorizedParticle(
      const std::array<double_v, DIMENSIONS> &position,
      const std::array<double_v, DIMENSIONS> &velocity,
      const std::array<double_v, DIMENSIONS> &force,
      const std::array<double_v, DIMENSIONS> &oldForce,
      const double_v &mass, const long_v &type,
      const double_mask &active) : position(position), force(force), mass(mass), type(type), active(active) {}
};

template<const size_t DIMENSIONS>
class Particles {
 public:
  std::array<std::vector<double>, DIMENSIONS> positions{};

  std::array<std::vector<double>, DIMENSIONS> velocities{};

  std::array<std::vector<double>, DIMENSIONS> forces{};

  std::array<std::vector<double>, DIMENSIONS> old_forces{};

  std::vector<double> mass{};

  std::vector<long> type{};
  std::vector<uint8_t> active{};
  std::vector<size_t> cell{};

#if DEBUG
  std::vector<size_t> ids{};
#endif

  auto store_force_single(VectorizedParticle<DIMENSIONS> &p1, size_t index) {
    for (int i = 0; i < DIMENSIONS; ++i) {
      forces[i][index] = p1.force[i][0];
    }
  }

  auto store_force_vector(VectorizedParticle<DIMENSIONS> &p, size_t index) {
    for (int i = 0; i < DIMENSIONS; ++i) {
      p.force[i].copy_to(&forces[i][index], stdx::element_aligned);
    }
  }

  auto load_vectorized_single(size_t index) -> VectorizedParticle<DIMENSIONS> {

    std::array<double_v, DIMENSIONS> position;
    for (int i = 0; i < DIMENSIONS; ++i) {
      position[i] = double_v(positions[i][index]);
    }
    std::array<double_v, DIMENSIONS> velocity;
    for (int i = 0; i < DIMENSIONS; ++i) {
      velocity[i] = double_v(velocities[i][index]);
    }
    std::array<double_v, DIMENSIONS> force;
    for (int i = 0; i < DIMENSIONS; ++i) {
      force[i] = double_v(forces[i][index]);
    }
    std::array<double_v, DIMENSIONS> old_force;
    for (int i = 0; i < DIMENSIONS; ++i) {
      old_force[i] = double_v(old_forces[i][index]);
    }

    return VectorizedParticle(
        position,
        velocity,
        force,
        old_force,
        double_v(mass[index]),
        long_v(type[index]),
        double_mask(static_cast<bool>(active[index])));
  }

  auto load_vectorized(size_t index) -> VectorizedParticle<DIMENSIONS> {
    std::array<double_v, DIMENSIONS> position_vector;
    for (int i = 0; i < DIMENSIONS; ++i) {
      position_vector[i] = double_v(&positions[i][index], stdx::element_aligned);
    }
    std::array<double_v, DIMENSIONS> velocity_vector;
    for (int i = 0; i < DIMENSIONS; ++i) {
      velocity_vector[i] = double_v(&velocities[i][index], stdx::element_aligned);
    }
    std::array<double_v, DIMENSIONS> force_vector;
    for (int i = 0; i < DIMENSIONS; ++i) {
      force_vector[i] = double_v(&forces[i][index], stdx::element_aligned);
    }
    std::array<double_v, DIMENSIONS> old_force_vector;
    for (int i = 0; i < DIMENSIONS; ++i) {
      old_force_vector[i] = double_v(&old_forces[i][index], stdx::element_aligned);
    }

    auto mass_vector = double_v();
    mass_vector.copy_from(&mass[index], stdx::element_aligned);

    auto type_vector = long_v(&type[index], stdx::element_aligned);

    auto active_vector = stdx::static_simd_cast<double_mask>(size_v(&active[index], stdx::element_aligned) > 0);

    return VectorizedParticle(position_vector,
                              velocity_vector,
                              force_vector,
                              old_force_vector,
                              mass_vector,
                              type_vector,
                              active_vector);
  }

  auto insert_particle(Particle<DIMENSIONS> p) {

    for (int i = 0; i < DIMENSIONS; ++i) {
      positions[i].emplace_back(p.position[i]);
    }
    for (int i = 0; i < DIMENSIONS; ++i) {
      velocities[i].emplace_back(p.velocity[i]);
    }
    for (int i = 0; i < DIMENSIONS; ++i) {
      forces[i].emplace_back(p.force[i]);
    }
    for (int i = 0; i < DIMENSIONS; ++i) {
      old_forces[i].emplace_back(p.old_force[i]);
    }

    mass.emplace_back(p.mass);
    type.emplace_back(p.type);
    active.emplace_back(true);
    cell.emplace_back(0);
#if DEBUG
    ids.emplace_back(size);
#endif
    size++;
  }
  size_t size{};

  template<typename Callable>
  void sort(Callable c) {

    for (int i = 0; i < size; ++i) {
      cell[i] = c(i);
    }
    // This is ugly but i haven't found a better way
    std::apply([&](auto &...ps) {
      return std::apply([&](auto &...vs) {
        return std::apply([&](auto &...fs) {
          return std::apply([&](auto &...ofs) {
            std::ranges::sort(std::ranges::zip_view(
                                  cell,
                                  ps...,
                                  vs...,
                                  fs...,
                                  ofs...,
                                  mass, type, active
#if DEBUG
                                  ,
                                  ids
#endif
                                  ),
                              [](auto tuple1, auto tuple2) {
                                return std::get<0>(tuple1) < std::get<0>(tuple2);
                              });
          },
                            old_forces);
        },
                          forces);
      },
                        velocities);
    },
               positions);
  }

  void
  swap_force() {

    std::swap(old_forces, forces);

    for (int d = 0; d < DIMENSIONS; ++d) {
      for (size_t i = 0; i < size; i++) {
        forces[d][i] = 0;
      }
    }
  }
};
