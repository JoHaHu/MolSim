#include "ParticleContainer.h"
#include <iostream>
#include <utility>

std::pair<Particle &, Particle &> ParticleContainer::PairIterator::operator*() {
  return {*i1, *i2};
}

ParticleContainer::PairIterator &ParticleContainer::PairIterator::operator++() {

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

ParticleContainer::PairIterator ParticleContainer::PairIterator::operator++(int) {
  PairIterator tmp = *this;
  ++(*this);
  return tmp;
}

bool ParticleContainer::PairIterator::operator==(PairIterator &i) {
  return i1 == i.i1 && i2 == i.i2 && end == i.end;
}
ParticleContainer::PairIterator::PairIterator(std::vector<Particle>::iterator i1,
                                              std::vector<Particle>::iterator i2,
                                              std::vector<Particle>::iterator end) : i1(i1), i2(i2), end(end) {

}

ParticleContainer::PairIterator ParticleContainer::begin_pair() {
  return {
      particles.begin(), std::next(particles.begin()), particles.end()
  };
}

ParticleContainer::PairIterator ParticleContainer::end_pair() {
  return {
      particles.end(), particles.end(), particles.end()
  };
}

ParticleContainer::ParticleContainer(int capacity) {
  particles = std::vector<Particle>(capacity);
}
ParticleContainer::Iterator ParticleContainer::begin() {
  return particles.begin();
}
ParticleContainer::Iterator ParticleContainer::end() {
  return particles.end();
}
unsigned long ParticleContainer::size() {
  return particles.size();
}
ParticleContainer::ParticleContainer(std::vector<Particle> particles) : particles(std::move(particles)) {

}
