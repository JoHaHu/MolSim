#include "ParticleContainer.h"
#include <iostream>
#include <utility>
#include "spdlog/spdlog.h"

/**
 * Dereferences the iterator to access the current pair of particles.
 * @return A pair consisting of references to two adjacent particles.
 */
auto ParticleContainer::PairIterator::operator*() -> std::pair<Particle &, Particle &> {
  return {*i1, *i2};
}

/**
 * Prefix increment operator. Advances the iterator to the next pair of particles.
 * @return Reference to the updated iterator.
 */
auto ParticleContainer::PairIterator::operator++() -> ParticleContainer::PairIterator & {
  if (i2 != end) ++i2;
  if (i2 == end && i1 != end) {
    ++i1;
    if (i1 != end) i2 = std::next(i1);
  }
  if (i1 == end) i2 = end;
  return *this;
}

/**
 * Postfix increment operator. Advances the iterator to the next pair of particles.
 * @return The iterator before it was incremented.
 */
auto ParticleContainer::PairIterator::operator++(int) -> ParticleContainer::PairIterator {
  PairIterator tmp = *this;
  ++(*this);
  return tmp;
}

/**
 * Checks if this iterator is equal to another iterator.
 */
auto ParticleContainer::PairIterator::operator==(PairIterator &i) -> bool {
  return i1 == i.i1 && i2 == i.i2 && end == i.end;
}

ParticleContainer::PairIterator::PairIterator(std::vector<Particle>::iterator i1,
                                              std::vector<Particle>::iterator i2,
                                              std::vector<Particle>::iterator end)
    : i1(i1), i2(i2), end(end) {
  spdlog::trace("PairIterator initialized.");
}

auto ParticleContainer::begin_pair() -> ParticleContainer::PairIterator {
  spdlog::trace("Creating begin_pair iterator.");
  return {particles.begin(), std::next(particles.begin()), particles.end()};
}

auto ParticleContainer::end_pair() -> ParticleContainer::PairIterator {
  spdlog::trace("Creating end_pair iterator.");
  return {particles.end(), particles.end(), particles.end()};
}

ParticleContainer::ParticleContainer(int capacity) {
  particles = std::vector<Particle>();
  particles.reserve(capacity);
  spdlog::info("ParticleContainer initialized with capacity: {}", capacity);
}

auto ParticleContainer::begin() -> ParticleContainer::Iterator {
  return particles.begin();
}

auto ParticleContainer::end() -> ParticleContainer::Iterator {
  return particles.end();
}

auto ParticleContainer::size() -> unsigned long {
  return particles.size();
}

ParticleContainer::ParticleContainer(std::vector<Particle> &particles) : particles(std::move(particles)) {
  spdlog::debug("ParticleContainer initialized with {} particles.", this->particles.size());
}
