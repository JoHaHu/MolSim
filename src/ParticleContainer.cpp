#include "ParticleContainer.h"
#include <iostream>

std::pair<Particle &, Particle &> ParticleContainer::Iterator::operator*() {
  return {*i1, *i2};
}

ParticleContainer::Iterator &ParticleContainer::Iterator::operator++() {

  if (++i2 <= end && i1 <= end) {
    return *this;
  } else if (++i1 <= end) {
    i2 = std::next(i1);
    return *this;
  } else {
    i1 = end;
    i2 = end;
    return *this;
  }
}

ParticleContainer::Iterator ParticleContainer::Iterator::operator++(int) {
  Iterator tmp = *this;
  ++(*this);
  return tmp;
}

bool ParticleContainer::Iterator::operator==(Iterator &i) {
  return i1 == i.i1 && i2 == i.i2 && end == i.end;
}
ParticleContainer::Iterator::Iterator(std::vector<Particle>::iterator i1,
                                      std::vector<Particle>::iterator i2,
                                      std::vector<Particle>::iterator end) : i1(i1), i2(i2), end(end) {

}

ParticleContainer::Iterator ParticleContainer::begin() {
  return {
      particles.begin(), std::next(particles.begin()), particles.end()
  };
}

ParticleContainer::Iterator ParticleContainer::end() {
  return {
      particles.end(), particles.end(), particles.end()
  };
}

ParticleContainer::ParticleContainer(int capacity) {
  particles = std::vector<Particle>(capacity);
}
