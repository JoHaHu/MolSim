/*
 * Particle.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "Particle.h"

#include "utils/ArrayUtils.h"
#include <iostream>
#include <sstream>

#include "spdlog/spdlog.h"
#include "utils/LoggerManager.h"

// initialising the static atomic counter for the IDs starting with 1
std::atomic<int> Particle::nextID{0};

Particle::Particle(int type_arg) : type(type_arg), id(nextID++) {
  spdlog::trace("Particle generated!");
  force = {0., 0., 0.};
  old_force = {0., 0., 0.};
  spdlog::debug("Particle initialized with type {}", type);
}

Particle::Particle(const Particle &other) : position(other.position), velocity(other.velocity), force(other.force), old_force(other.old_force), mass(other.mass), type(other.type), id(nextID++) {

  spdlog::trace("Particle generated by copy!");
  spdlog::debug("Copied particle: position = ({}, {}, {}), velocity = ({}, {}, {}), force = ({}, {}, {}), old force = ({}, {}, {}), mass = {}, type = {}",
                position[0], position[1], position[2], velocity[0], velocity[1], velocity[2], force[0], force[1], force[2], old_force[0], old_force[1], old_force[2], mass, type);
}

//! A constructor for a particle object.
/*!
  The constructor takes four parameters that set the properties of the particle object
  \param x_arg array of integers that specify the coordinates of the particle
  \param v_arg array of integers that specify the velocity of each axis
  \param m_arg double value that specifies the mass of the particle
  \param type_arg integer value that specifies the type of particle
  \return the particle object
*/
Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg)
    : position(x_arg), velocity(v_arg), force({0., 0., 0.}), old_force({0., 0., 0.}), mass(m_arg), type(type_arg), id(nextID++) {
  spdlog::trace("Particle generated!");
  spdlog::debug("Particle initialized: position = ({}, {}, {}), velocity = ({}, {}, {}), mass = {}, type = {}",
                position[0], position[1], position[2], velocity[0], velocity[1], velocity[2], mass, type);
}

//! A destructor for a particle object.
/*!
  The destructor deletes a particle object
*/
Particle::~Particle() {
  spdlog::trace("Particle destructed!");
}

auto Particle::to_string() const -> std::string {
  std::stringstream stream;
  stream << "Particle: X:" << position << " v: " << velocity << " f: " << force << " old_force: " << old_force << " type: " << type;
  return stream.str();
}

auto Particle::operator==(Particle &other) const -> bool {
  return (position == other.position) and (velocity == other.velocity) and (force == other.force) and (type == other.type) and (mass == other.mass)
      and (old_force == other.old_force);
}

auto operator<<(std::ostream &stream, Particle &p) -> std::ostream & {
  stream << p.to_string();
  return stream;
}
