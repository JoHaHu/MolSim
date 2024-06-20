#pragma once

#include <array>
#include <iostream>
#include <sstream>
#include <string>

/*
 * @class DoubleHelix can store the information that defines a Double Helix, namely the coordinates, its initial velocity, its radius, pitch and height.
 * Also contains a string returning function to enable the printing of cuboid specifications for testing and debugging.
 *
 * */

class DoubleHelix {
 public:
  std::vector<double> coordinates;
  std::vector<double> velocity;
  double radius;
  double pitch;
  double height;

  DoubleHelix(std::vector<double> &coordinates, std::vector<double> &velocity, double radius, double pitch, double height)
      : coordinates(coordinates),
        velocity(velocity),
        radius(radius),
        pitch(pitch),
        height(height) {}

  std::string toString() {
    std::string output;
    output.append("Coordinates: " + arrayToString(coordinates));
    output.append(" | Velocity: " + arrayToString(velocity));
    output.append(" | Radius: " + std::to_string(radius));
    output.append(" | Pitch: " + std::to_string(pitch));
    output.append(" | Height: " + std::to_string(height));
    return output;
  }

  auto arrayToString(const std::vector<double> &array) -> std::string {
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
