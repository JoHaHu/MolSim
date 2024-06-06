//
// Created by TimSc on 05.06.2024.
//

#include <array>
#include <sstream>
#include <string>

#ifndef PSEMOLDYN_GROUPD_DISC_H
#define PSEMOLDYN_GROUPD_DISC_H

/**
 * @class Disc
 * @brief Can store the information that defines a disc for the Lennard Jones force simulation, namely the coordinates, its velocity, and radius.
 * Also contains a string returning function to enable the printing of disc specifications for testing and debugging.
 */
class Disc {
 public:
  std::array<double, 3> coordinates;
  std::array<double, 3> velocity;
  int radius;

  // Default constructor
  Disc() : coordinates({0.0, 0.0, 0.0}), velocity({0.0, 0.0, 0.0}), radius(0) {}

  // Parameterized constructor
  Disc(const std::array<double, 3> &coordinates, const std::array<double, 3> &velocity, int radius)
      : coordinates(coordinates), velocity(velocity), radius(radius) {}

  // Copy constructor
  Disc(const Disc &other) = default;

  // Move constructor
  Disc(Disc &&other) noexcept = default;

  // Copy assignment operator
  Disc &operator=(const Disc &other) = default;

  // Move assignment operator
  Disc &operator=(Disc &&other) noexcept = default;

  // Destructor
  ~Disc() = default;

  std::string toString() const {
    std::string output;
    output.append("Coordinates: " + arrayToString(coordinates));
    output.append(" | Velocity: " + arrayToString(velocity));
    output.append(" | Radius: " + std::to_string(radius));
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

#endif// PSEMOLDYN_GROUPD_DISC_H
