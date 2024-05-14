/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <array>
#include <string>
#include <atomic>
#include <iostream>

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

  /**
   * Static atomic ID counter shared by all objects that assigns the ID value to a new object when created
   */
  static std::atomic<int> nextID;

  /**
   * ID of the particle to differentiate between particles and compare if a particle is the same
   */
  int id;

  explicit Particle(int type = 0);

  Particle(const Particle &other);

  Particle(const Particle &&other) = delete;

  Particle(
      // for visualization, we need always 3 coordinates
      // -> in case of 2d, we use only the first and the second
      std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg,
      int type = 0);

  ~Particle();

  auto operator==(Particle &other) const -> bool;

  auto operator=(Particle const &other) -> Particle & = default;
  auto operator=(Particle &&other) -> Particle & = default;

  [[nodiscard]] auto to_string() const -> std::string;
};

auto operator<<(std::ostream &stream, Particle &p) -> std::ostream &;
