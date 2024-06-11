#pragma once

#include <array>
#include <experimental/simd>
#include <ranges>
#include <vector>

namespace stdx {
using namespace std::experimental;
using namespace std::experimental::__proposed;
}// namespace stdx

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
using double_v = stdx::native_simd<double>;
using double_mask = stdx::native_simd_mask<double>;
using int_v = stdx::native_simd<int>;
using size_v = stdx::native_simd<size_t>;

struct VectorizedParticle {
  std::array<double_v, 3> position{};
  std::array<double_v, 3> velocity{};
  std::array<double_v, 3> force{};
  std::array<double_v, 3> old_force{};
  double_v mass{};
  int_v type;
  size_v active;

  VectorizedParticle(const std::array<double_v, 3> &position, const std::array<double_v, 3> &velocity, const std::array<double_v, 3> &force, const std::array<double_v, 3> &oldForce, const double_v &mass, const int_v &type, const size_v &active) : position(position), velocity(velocity), force(force), old_force(oldForce), mass(mass), type(type), active(active) {}
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

  auto store_force_single(VectorizedParticle p1, size_t index) {
    force_x[index] = p1.force[0][0];
    force_y[index] = p1.force[1][0];
    force_z[index] = p1.force[2][0];
  }

  auto store_force_vector(VectorizedParticle p, size_t index, double_mask mask) {
    p.force[0].copy_to(&force_x[index], stdx::element_aligned);
    p.force[1].copy_to(&force_y[index], stdx::element_aligned);
    p.force[2].copy_to(&force_z[index], stdx::element_aligned);
  }

  auto load_vectorized_single(size_t index) -> VectorizedParticle {
    return VectorizedParticle({double_v(position_x[index]), double_v(position_y[index]), double_v(position_z[index])},
                              {double_v(velocity_x[index]), double_v(velocity_y[index]), double_v(velocity_z[index])},
                              {double_v(force_x[index]), double_v(force_y[index]), double_v(force_z[index])},
                              {double_v(old_force_x[index]), double_v(old_force_y[index]), double_v(old_force_z[index])},
                              double_v(mass[index]),
                              int_v(type[index]),
                              size_v(active[index]));
  }
  auto load_vectorized(size_t index) -> VectorizedParticle {
    auto pos_x_vector = double_v();
    pos_x_vector.copy_from(&position_x[index], stdx::element_aligned);
    auto pos_y_vector = double_v();
    pos_y_vector.copy_from(&position_y[index], stdx::element_aligned);
    auto pos_z_vector = double_v();
    pos_z_vector.copy_from(&position_z[index], stdx::element_aligned);

    auto velocity_x_vector = double_v();
    velocity_x_vector.copy_from(&velocity_x[index], stdx::element_aligned);
    auto velocity_y_vector = double_v();
    velocity_y_vector.copy_from(&velocity_y[index], stdx::element_aligned);
    auto velocity_z_vector = double_v();
    velocity_z_vector.copy_from(&velocity_z[index], stdx::element_aligned);

    auto force_x_vector = double_v();
    force_x_vector.copy_from(&force_x[index], stdx::element_aligned);
    auto force_y_vector = double_v();
    force_y_vector.copy_from(&force_y[index], stdx::element_aligned);
    auto force_z_vector = double_v();
    force_z_vector.copy_from(&force_z[index], stdx::element_aligned);

    auto old_force_x_vector = double_v();
    old_force_x_vector.copy_from(&old_force_x[index], stdx::element_aligned);
    auto old_force_y_vector = double_v();
    old_force_y_vector.copy_from(&old_force_y[index], stdx::element_aligned);
    auto old_force_z_vector = double_v();
    old_force_z_vector.copy_from(&old_force_z[index], stdx::element_aligned);

    auto mass_vector = double_v();
    mass_vector.copy_from(&mass[index], stdx::element_aligned);

    auto type_vector = int_v();
    type_vector.copy_from(&type[index], stdx::element_aligned);

    auto active_vector = size_v();
    active_vector.copy_from(&active[index], stdx::element_aligned);

    return VectorizedParticle({pos_x_vector, pos_y_vector, pos_z_vector},
                              {velocity_x_vector, velocity_y_vector, velocity_z_vector},
                              {force_x_vector, force_y_vector, force_z_vector},
                              {old_force_x_vector, old_force_y_vector, old_force_z_vector},
                              mass_vector,
                              type_vector,
                              active_vector);
  }

  auto insert_particle(Particle p) -> size_t {
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
    return size - 1;
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
