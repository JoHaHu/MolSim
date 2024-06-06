//
// Created by TimSc on 05.06.2024.
//

#include <array>
#include <iostream>
#include <sstream>
#include <string>

#ifndef PSEMOLDYN_GROUPD_CUBOID_H
#define PSEMOLDYN_GROUPD_CUBOID_H

/**
 * @class Cuboid
 * @brief Can store the information that defines a cuboid, namely the coordinates, its particle counts on all dimensions, and its initial velocity.
 * Also contains a string returning function to enable the printing of cuboid specifications for testing and debugging.
 */
class Cuboid {
 public:
  std::array<double, 3> coordinates;
  std::array<double, 3> particles;
  std::array<double, 3> velocity;

  // Default constructor
  Cuboid() : coordinates({0.0, 0.0, 0.0}), particles({0.0, 0.0, 0.0}), velocity({0.0, 0.0, 0.0}) {}

  // Parameterized constructor
  Cuboid(const std::array<double, 3> &coordinates, const std::array<double, 3> &particles, const std::array<double, 3> &velocity)
      : coordinates(coordinates), particles(particles), velocity(velocity) {}

  // Copy constructor
  Cuboid(const Cuboid &other) = default;

  // Move constructor
  Cuboid(Cuboid &&other) noexcept = default;

  // Copy assignment operator
  Cuboid &operator=(const Cuboid &other) = default;

  // Move assignment operator
  Cuboid &operator=(Cuboid &&other) noexcept = default;

  // Destructor
  ~Cuboid() = default;

  std::string toString() const {
    std::string output;
    output.append("Coordinates: " + arrayToString(coordinates));
    output.append(" | Particle counts: " + arrayToString(particles));
    output.append(" | Velocity: " + arrayToString(velocity));
    return output;
  }

  template<std::size_t N>
  std::string arrayToString(const std::array<double, N> &array) const {
    std::ostringstream oss;
    oss << "[";
    for (std::size_t i = 0; i < array.size(); ++i) {
      oss << array[i];
      if (i < array.size() - 1) {
        oss << ", ";
      }
    }
    oss << "]";
    return oss.str();
  }
};

#endif// PSEMOLDYN_GROUPD_CUBOID_H
