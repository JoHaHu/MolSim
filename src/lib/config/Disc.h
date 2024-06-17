#pragma once
#include <array>
#include <sstream>
#include <string>


/**
 * @class Disc
 * @brief Can store the information that defines a disc, namely the coordinates, its velocity, and radius.
 * Also contains a string returning function to enable the printing of disc specifications for testing and debugging.
 */
class Disc {
 public:
  std::vector<double> coordinates;
  std::vector<double> velocity;
  int radius;

  // Default constructor
  Disc() : coordinates({0.0, 0.0, 0.0}), velocity({0.0, 0.0, 0.0}), radius(0) {}

  // Parameterized constructor
  Disc(const std::vector<double> &coordinates, const std::vector<double> &velocity, int radius)
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

  std::string arrayToString(const std::vector<double> &array) const {
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
