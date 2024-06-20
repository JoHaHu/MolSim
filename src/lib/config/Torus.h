#pragma once
#include <array>
#include <iostream>
#include <sstream>
#include <string>

/**
* @class can store the information that defines a torus, namely the coordinates, its initial velocity and major and minor radius.
   * Also contains a string returning function to enable the printing of cuboid specifications for testing and debugging
   *
   * */

class Torus {
 public:
  std::vector<double> coordinates;
  std::vector<double> velocity;
  double major_radius;
  double minor_radius;

  Torus(std::vector<double> &coordinates, std::vector<double> &velocity, double major_radius, double minor_radius)
      : coordinates(coordinates),
        velocity(velocity),
        major_radius(major_radius),
        minor_radius(minor_radius) {}

  std::string toString() {
    std::string output;
    output.append("Coordinates: " + arrayToString(coordinates));
    output.append(" | Velocity: " + arrayToString(velocity));
    output.append(" | Major radius: " + std::to_string(major_radius));
    output.append(" | Minor radius: " + std::to_string(minor_radius));
    return output;
  }

  auto arrayToString(std::vector<double> &array) -> std::string {
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
