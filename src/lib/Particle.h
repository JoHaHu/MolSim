#pragma once

#include <array>
#include <ranges>
#include <vector>
/**
 * the Particle class representing the simulated particles
 * */
class Particle {
 public:
  /**
   * Position of the particle
   */
  std::array<double, 3> position{};

  /**
   * Velocity of the particle
   */
  std::array<double, 3> velocity{};

  /**
   * Force effective on this particle
   */
  std::array<double, 3> force{};

  /**
   * Force which was effective on this particle
   */
  std::array<double, 3> old_force{};

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
      std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg,
      int type = 0) : position(x_arg), velocity(v_arg), mass(m_arg), type(type)

  {
  }
};

class Particles {
 public:
  std::vector<double> position_x{};
  std::vector<double> position_y{};
  std::vector<double> position_z{};

  std::vector<double> velocity_x{};
  std::vector<double> velocity_y{};
  std::vector<double> velocity_z{};

  std::vector<double> force_x{};
  std::vector<double> force_y{};
  std::vector<double> force_z{};

  std::vector<double> old_force_x{};
  std::vector<double> old_force_y{};
  std::vector<double> old_force_z{};

  std::vector<double> mass{};

  std::vector<int> type{};
  std::vector<size_t> active{};

  void insert_particle(Particle p) {
    position_x.emplace_back(p.position[0]);
    position_y.emplace_back(p.position[1]);
    position_z.emplace_back(p.position[2]);

    velocity_x.emplace_back(p.velocity[0]);
    velocity_y.emplace_back(p.velocity[1]);
    velocity_z.emplace_back(p.velocity[2]);

    force_x.emplace_back(p.force[0]);
    force_y.emplace_back(p.force[1]);
    force_z.emplace_back(p.force[2]);

    old_force_x.emplace_back(p.old_force[0]);
    old_force_y.emplace_back(p.old_force[1]);
    old_force_z.emplace_back(p.old_force[2]);

    mass.emplace_back(p.mass);
    type.emplace_back(p.type);
    active.emplace_back(1);
    size++;
  }
  size_t size{};

  void swap_force() {
    std::swap(old_force_x, force_x);
    std::swap(old_force_y, force_y);
    std::swap(old_force_z, force_z);

    for (size_t i = 0; i < size; i++) {
      force_x[i] = 0;
      force_y[i] = 0;
      force_z[i] = 0;
    }
  }
};
