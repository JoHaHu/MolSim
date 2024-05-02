#pragma once

#include "Particle.h"
#include <vector>

class ParticleContainer {
public:
 using Iterator = std::vector<Particle>::iterator;

 class PairIterator {
 public:
  /**
   * Constructs a new PairIterator.
   * @param i1 Iterator pointing to the first particle in a pair.
   * @param i2 Iterator pointing to the second particle in a pair.
   * @param end Iterator representing the end of the particle collection.
   */
  PairIterator(
   std::vector<Particle>::iterator i1,
   std::vector<Particle>::iterator i2,
   std::vector<Particle>::iterator end
  );

  using difference_type = int;

  std::pair<Particle &, Particle &> operator*();

  PairIterator &operator++();

  PairIterator operator++(int);

  bool operator==(PairIterator &);

 private:
  std::vector<Particle>::iterator i1;
  std::vector<Particle>::iterator i2;
  std::vector<Particle>::iterator end;
 };

 Iterator begin();

 Iterator end();

 unsigned long size();

 PairIterator begin_pair();

 PairIterator end_pair();

 /**
  * Constructor that reserves space for a specific number of particles.
  * @param capacity The number of particles for which memory should be preallocated.
  */
 explicit ParticleContainer(int capacity);

 /**
  * Constructor that initializes the container with a list of particles.
  * @param particles A vector of particles to initialize the container.
  */
 explicit ParticleContainer(std::vector<Particle> particles);

private:
 std::vector<Particle> particles;
};

//static_assert(std::output_iterator<ParticleContainer::PairIterator, std::pair<Particle &, Particle &>>);
