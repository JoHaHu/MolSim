//
// Created by TimSc on 05.06.2024.
//

#include <array>
#include <sstream>
#include <string>

#ifndef PSEMOLDYN_GROUPD_CELESTIALBODY_H
#define PSEMOLDYN_GROUPD_CELESTIALBODY_H

/**
   * @class can store the information that defines a celestial body for the gravitational simulation, namely the coordinates, its velocity and mass
   * also contains a string returning function to enable the printing of celestial body specifications for testing and debugging
   *
   * */
class CelestialBody {

 public:
  CelestialBody(const std::array<double, 3> &coordinates, const std::array<double, 3> &velocity, double mass) : coordinates(coordinates), velocity(velocity), mass(mass) {}

  std::array<double, 3> coordinates;
  std::array<double, 3> velocity;
  double mass;

  std::string toString() {
    std::string output;
    output.append("Coordinates: " + arrayToString(coordinates));
    output.append(" | Velocity: " + arrayToString(velocity));
    output.append(" | Mass: " + std::to_string(mass));

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

#endif//PSEMOLDYN_GROUPD_CELESTIALBODY_H
