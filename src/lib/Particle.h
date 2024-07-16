#pragma once

#include "range/v3/algorithm.hpp"
#include "range/v3/view/zip.hpp"
#include <array>
#include <execution>
#include <experimental/simd>
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

  /**
     * Position flexibility. Specifies, whether the particle is fixed in position,
     * e.g. being part of a wall like in nano tubes or free moving
     */
  uint8_t fixed{};

  Particle(const std::array<double, DIMENSIONS> &position,
           const std::array<double, DIMENSIONS> &velocity,
           const std::array<double, DIMENSIONS> &force,
           const std::array<double, DIMENSIONS> &old_force,
           double mass,
           int type,
           int fixed) : position(position), velocity(velocity), force(force), old_force(old_force), mass(mass),
                        type(type), fixed(fixed) {}

  Particle(
      // for visualization, we need always 3 coordinates
      // -> in case of 2d, we use only the first and the second
      std::array<double, DIMENSIONS> x_arg, std::array<double, DIMENSIONS> v_arg, double m_arg,
      int type = 0, int fixed = 0) : position(x_arg), velocity(v_arg), mass(m_arg), type(type), fixed(fixed) {
    for (size_t i = 0; i < DIMENSIONS; ++i) {
      force[i] = 0.0;
      old_force[i] = 0.0;
    }
  }
};

/***
 * A structure of Arrays for the particles
 *
 * */
template<const size_t DIMENSIONS>
class Particles {
 public:
  std::array<std::vector<double>, DIMENSIONS> positions{};

  std::array<std::vector<double>, DIMENSIONS> velocities{};

  std::array<std::vector<double>, DIMENSIONS> forces{};

  std::array<std::vector<double>, DIMENSIONS> old_forces{};

  std::vector<double> mass{};

  std::vector<long> type{};
  std::vector<uint8_t> fixed{};
  std::vector<uint8_t> active{};
  std::vector<size_t> cell{};
  std::vector<size_t> block{};
  std::vector<size_t> color{};

#if DEBUG
  /**
     * particle ids only used with debugging. Helps to set debug points for suspicious particles
     * */
  std::vector<size_t> ids{};
#endif

  auto insert_particle(Particle<DIMENSIONS> p) {

    for (size_t i = 0; i < DIMENSIONS; ++i) {
      positions[i].emplace_back(p.position[i]);
    }
    for (size_t i = 0; i < DIMENSIONS; ++i) {
      velocities[i].emplace_back(p.velocity[i]);
    }
    for (size_t i = 0; i < DIMENSIONS; ++i) {
      forces[i].emplace_back(p.force[i]);
    }
    for (size_t i = 0; i < DIMENSIONS; ++i) {
      old_forces[i].emplace_back(p.old_force[i]);
    }

    mass.emplace_back(p.mass);
    type.emplace_back(p.type);
    fixed.emplace_back(p.fixed);
    active.emplace_back(true);
    cell.emplace_back(0);
    block.emplace_back(0);
    color.emplace_back(0);
#if DEBUG
    ids.emplace_back(size);
#endif
    size++;
  }

  size_t size{};

  template<typename Callable>
  void sort(Callable c) {
    // recaulcates the cell
    for (size_t i = 0; i < size; ++i) {
      auto [_color, _block, _cell] = c(i);
      color[i] = _color;
      block[i] = _block;
      cell[i] = _cell;
    }
    // This is ugly but i haven't found a better way
    std::apply([&](auto &...ps) {
      return std::apply([&](auto &...vs) {
        return std::apply([&](auto &...fs) {
          return std::apply([&](auto &...ofs) {
            auto zip = ranges::zip_view(
                cell,
                block,
                color,
                ps...,
                vs...,
                fs...,
                ofs...,
                mass, type, active,fixed
#if DEBUG
                ,
                ids
#endif
            );
            ranges::sort(zip,
                         [](auto tuple1, auto tuple2) {
                           // Do not sort by color for better cache behaviour
                           if (std::get<1>(tuple1) != std::get<1>(tuple2)) {
                             return std::get<1>(tuple1) < std::get<1>(tuple2);
                           }
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

    for (size_t d = 0; d < DIMENSIONS; ++d) {
      for (size_t i = 0; i < size; i++) {
        forces[d][i] = 0;
      }
    }
  }
};
