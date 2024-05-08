#include "ParticleContainer.h"
#include <iostream>
#include <utility>

/**
 * Dereferences the iterator to access the current pair of particles.
 * @return A pair consisting of references to two adjacent particles.
 */
std::pair<Particle &, Particle &> ParticleContainer::PairIterator::operator*() {
    return {*i1, *i2};
}

/**
 * Prefix increment operator. Advances the iterator to the next pair of particles.
 * @return Reference to the updated iterator.
 */
ParticleContainer::PairIterator &ParticleContainer::PairIterator::operator++() {
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
ParticleContainer::PairIterator ParticleContainer::PairIterator::operator++(int) {
    PairIterator tmp = *this;
    ++(*this);
    return tmp;
}

/**
 * Checks if this iterator is equal to another iterator.
 */
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
