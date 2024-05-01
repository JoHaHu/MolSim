#ifndef PARTICLECONTAINER_H
#define PARTICLECONTAINER_H

#include "Particle.h"
#include <vector>

class ParticleContainer {
public:
    void addParticle(const Particle& particle);
    void printPairs() const;

    using Iterator = std::vector<Particle>::iterator;
    using ConstIterator = std::vector<Particle>::const_iterator;

    Iterator begin();
    Iterator end();
    [[nodiscard]] ConstIterator begin() const;
    [[nodiscard]] ConstIterator end() const;

private:
    std::vector<Particle> particles;
};

#endif // PARTICLECONTAINER_H
