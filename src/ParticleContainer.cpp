#include "ParticleContainer.h"
#include <iostream>

void ParticleContainer::addParticle(const Particle& particle) {
    particles.push_back(particle);
}

void ParticleContainer::printPairs() const {
    for (auto i = particles.begin(); i != particles.end(); ++i) {
        for (auto j = std::next(i); j != particles.end(); ++j) {
            std::cout << "Particle Pair: (" << i->toString() << ") and (" << j->toString() << ")\n";
        }
    }
}

ParticleContainer::Iterator ParticleContainer::begin() {
    return particles.begin();
}

ParticleContainer::Iterator ParticleContainer::end() {
    return particles.end();
}

ParticleContainer::ConstIterator ParticleContainer::begin() const {
    return particles.begin();
}

ParticleContainer::ConstIterator ParticleContainer::end() const {
    return particles.end();
}
