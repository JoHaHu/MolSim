//
// Created by TimSc on 06.06.2024.
//

#include <array>
#include <iostream>
#include <sstream>
#include <string>

#ifndef PSEMOLDYN_GROUPD_DOUBLEHELIX_H
#define PSEMOLDYN_GROUPD_DOUBLEHELIX_H

/**
   * @class can store the information that defines a torus, namely the coordinates, its initial velocity and major and minor radius
   * also contains a string returning function to enable the printing of cuboid specifications for testing and debugging
   *
   * */

class DoubleHelix {
 public:
  std::array<double, 3> coordinates;
  std::array<double, 3> velocity;
  double radius;
  double pitch;
  double height;

  DoubleHelix(const std::array<double, 3> &coordinates, const std::array<double, 3> &velocity, double radius, double pitch, double height)
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

  template<std::size_t N>
  auto arrayToString(const std::array<double, N> &array) -> std::string {
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

#endif//PSEMOLDYN_GROUPD_DOUBLEHELIX_H
