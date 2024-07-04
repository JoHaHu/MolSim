
#pragma once

#include <array>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

/**
 * @class Cuboid
 * @brief Can store the information that defines a cuboid, namely the coordinates, its particle counts on all dimensions, and its initial velocity.
 * Also contains a string returning function to enable the printing of cuboid specifications for testing and debugging.
 */
class Cuboid {
public:
    std::vector<double> coordinates{};
    std::vector<double> particles{};
    std::vector<double> velocity{};
    double mass;
    double spacing;
    int type;
    bool fixed;

    // Parameterized constructor
    Cuboid(std::vector<double> &coordinates,
           std::vector<double> &particles,
           std::vector<double> &velocity,
           double mass, double spacing, int type)
            : coordinates(coordinates), particles(particles), velocity(velocity), mass(mass), spacing(spacing),
              type(type), fixed(false) {};

    // Constructor for fixed cuboids (additional parameter)
    Cuboid(std::vector<double> &coordinates,
           std::vector<double> &particles,
           std::vector<double> &velocity,
           double mass, double spacing, int type, bool fixed)
            : coordinates(coordinates), particles(particles), velocity(velocity), mass(mass), spacing(spacing),
              type(type), fixed(fixed) {};

    std::string toString() const {
        std::string output;
        output.append("Coordinates: " + arrayToString(coordinates));
        output.append(" | Particle counts: " + arrayToString(particles));
        output.append(" | Velocity: " + arrayToString(velocity));
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
