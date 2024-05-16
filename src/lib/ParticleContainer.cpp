#include "ParticleContainer.h"
#include "spdlog/spdlog.h"
#include <utility>

/**
 * @brief Dereferences the iterator to access the current pair of particles.
 * @return A pair consisting of references to two adjacent particles.
 */
auto ParticleContainer::PairIterator::operator*() -> std::pair<Particle &, Particle &> {
  return {*i1, *i2};
}

/**
 * @brief Prefix increment operator. Advances the iterator to the next pair of particles.
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
 * @brief Postfix increment operator. Advances the iterator to the next pair of particles.
 * @return The iterator before it was incremented.
 */
auto ParticleContainer::PairIterator::operator++(int) -> ParticleContainer::PairIterator {
  PairIterator tmp = *this;
  ++(*this);
  return tmp;
}

/**
 * @brief Checks if this iterator is equal to another iterator.
 * @param i The iterator to compare with.
 * @return True if the iterators are equal, false otherwise.
 */
auto ParticleContainer::PairIterator::operator==(PairIterator &i) -> bool {
  return i1 == i.i1 && i2 == i.i2 && end == i.end;
}

/**
 * @brief Constructs a PairIterator.
 * @param i1 Iterator to the first particle.
 * @param i2 Iterator to the second particle.
 * @param end Iterator to the end of the particle container.
 */
ParticleContainer::PairIterator::PairIterator(std::vector<Particle>::iterator i1,
                                              std::vector<Particle>::iterator i2,
                                              std::vector<Particle>::iterator end)
    : i1(i1), i2(i2), end(end) {
  spdlog::trace("PairIterator initialized.");
}

/**
 * @brief Creates an iterator to the beginning pair of particles.
 * @return A PairIterator to the beginning pair of particles.
 */
auto ParticleContainer::begin_pair() -> ParticleContainer::PairIterator {
  spdlog::trace("Creating begin_pair iterator.");
  return {particles.begin(), std::next(particles.begin()), particles.end()};
}

/**
 * @brief Creates an iterator to the end pair of particles.
 * @return A PairIterator to the end pair of particles.
 */
auto ParticleContainer::end_pair() -> ParticleContainer::PairIterator {
  spdlog::trace("Creating end_pair iterator.");
  return {particles.end(), particles.end(), particles.end()};
}

/**
 * @brief Constructs a ParticleContainer with a specified capacity.
 * @param capacity The initial capacity of the container.
 */
ParticleContainer::ParticleContainer(int capacity) {
  particles = std::vector<Particle>();
  particles.reserve(capacity);
  spdlog::debug("ParticleContainer initialized with capacity: {}", capacity);
}

/**
 * @brief Returns an iterator to the beginning of the particles.
 * @return An iterator to the beginning of the particles.
 */
auto ParticleContainer::begin() -> ParticleContainer::Iterator {
  return particles.begin();
}

/**
 * @brief Returns an iterator to the end of the particles.
 * @return An iterator to the end of the particles.
 */
auto ParticleContainer::end() -> ParticleContainer::Iterator {
  return particles.end();
}

/**
 * @brief Returns the number of particles in the container.
 * @return The number of particles.
 */
auto ParticleContainer::size() -> unsigned long {
  return particles.size();
}

/**
 * @brief Constructs a ParticleContainer with a given vector of particles.
 * @param particles A vector of particles to initialize the container with.
 */
ParticleContainer::ParticleContainer(std::vector<Particle> particles) : particles(std::move(particles)) {
  spdlog::debug("ParticleContainer initialized with {} particles.", this->particles.size());
}