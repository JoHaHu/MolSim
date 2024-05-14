#include "ParticleContainer.h"
#include "spdlog/spdlog.h"
#include <iostream>
#include <utility>

/**
 * Dereferences the iterator to access the current pair of particles.
 * @return A pair of references to two adjacent particles.
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
    if (i1 != end) {
      i2 = std::next(i1);
      if (i2 == end) {
        i1 = end;
      }
    } else {
      i2 = end;
    }
  }
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
 * @param i The other iterator to compare with.
 * @return True if the iterators are equal, otherwise false.
 */
auto ParticleContainer::PairIterator::operator==(PairIterator &i) -> bool {
  return i1 == i.i1 && i2 == i.i2 && end == i.end;
}

/**
 * Constructor for PairIterator.
 * @param i1 Iterator to the first particle in the pair.
 * @param i2 Iterator to the second particle in the pair.
 * @param end Iterator to the end of the container.
 */
ParticleContainer::PairIterator::PairIterator(std::vector<Particle>::iterator i1,
                                              std::vector<Particle>::iterator i2,
                                              std::vector<Particle>::iterator end)
    : i1(i1), i2(i2), end(end) {
  spdlog::trace("PairIterator initialized.");
}

/**
 * Creates an iterator to the beginning pair of particles.
 * @return An iterator to the first pair of particles.
 */
auto ParticleContainer::begin_pair() -> ParticleContainer::PairIterator {
  spdlog::trace("Creating begin_pair iterator.");
  return {particles.begin(), std::next(particles.begin()), particles.end()};
}

/**
 * Creates an iterator to the end pair of particles.
 * @return An iterator to the end of the container.
 */
auto ParticleContainer::end_pair() -> ParticleContainer::PairIterator {
  spdlog::trace("Creating end_pair iterator.");
  return {particles.end(), particles.end(), particles.end()};
}

/**
 * Constructs a ParticleContainer with a specified initial capacity.
 * @param capacity The initial capacity of the container.
 */
ParticleContainer::ParticleContainer(int capacity) {
  particles = std::vector<Particle>();
  particles.reserve(capacity);
  spdlog::info("ParticleContainer initialized with capacity: {}", capacity);
}

/**
 * Returns an iterator to the beginning of the container.
 * @return An iterator to the first particle.
 */
auto ParticleContainer::begin() -> ParticleContainer::Iterator {
  return particles.begin();
}

/**
 * Returns an iterator to the end of the container.
 * @return An iterator past the last particle.
 */
auto ParticleContainer::end() -> ParticleContainer::Iterator {
  return particles.end();
}

/**
 * Returns the number of particles in the container.
 * @return The size of the container.
 */
auto ParticleContainer::size() -> unsigned long {
  return particles.size();
}

/**
 * Constructs a ParticleContainer with an initial set of particles.
 * @param particles A vector of particles to initialize the container with.
 */
ParticleContainer::ParticleContainer(std::vector<Particle> &particles) : particles(std::move(particles)) {
  spdlog::debug("ParticleContainer initialized with {} particles.", this->particles.size());
}
