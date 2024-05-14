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

Particle::Particle(int type_arg) {
    type = type_arg;
    spdlog::trace("Particle generated!");
    f = {0., 0., 0.};
    old_f = {0., 0., 0.};
}

Particle::Particle(const Particle &other) {
    x = other.x;
    v = other.v;
    f = other.f;
    old_f = other.old_f;
    m = other.m;
    type = other.type;

    spdlog::trace("Particle generated by copy!");
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
    : x(x_arg), v(v_arg), f({0., 0., 0.}), old_f({0., 0., 0.}), m(m_arg), type(type_arg) {
    spdlog::trace("Particle generated!");
}

//! A destructor for a particle object.
/*!
  The destructor deletes a particle object
*/
Particle::~Particle() {
    spdlog::trace("Particle destructed!");
}


std::string Particle::toString() const {
    std::stringstream stream;
    stream << "Particle: X:" << x << " v: " << v << " f: " << f << " old_f: " << old_f << " type: " << type;
    return stream.str();
}

bool Particle::operator==(Particle &other) const {
    return (x == other.x) and (v == other.v) and (f == other.f) and (type == other.type) and (m == other.m)
           and (old_f == other.old_f);
}

std::ostream &operator<<(std::ostream &stream, Particle &p) {
    stream << p.toString();
    return stream;
}
