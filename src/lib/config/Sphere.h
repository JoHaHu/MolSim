#pragma once
#include <array>
#include <iostream>
#include <sstream>
#include <string>

/**
* @class Sphere can store the information that defines a sphere, namely the coordinates, its initial velocity and radius.
   * Also contains a string returning function to enable the printing of cuboid specifications for testing and debugging
   *
   * */

class Sphere {
 public:
  std::vector<double> coordinates;
  std::vector<double> velocity;
  int radius;
  double mesh_width;

  Sphere(const std::vector<double> &coordinates, const std::vector<double> &velocity, int radius, double mesh_width)
      : coordinates(coordinates), velocity(velocity), radius(radius), mesh_width(mesh_width) {}

  std::string toString() {
    std::string output;
    output.append("Coordinates: " + arrayToString(coordinates));
    output.append(" | Velocity: " + arrayToString(velocity));
    output.append(" | Radius: " + std::to_string(radius));
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