#pragma once

#include "Particle.h"
#include <vector>

class ParticleContainer {

 public:

  using Iterator = std::vector<Particle>::iterator;

  class PairIterator {
   public:
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

  /***
    * @param capacity Indicates the capacity of the underlying backing storage
    * */
  explicit ParticleContainer(int capacity);
  explicit ParticleContainer(std::vector<Particle> particles);
 private:
  std::vector<Particle> particles;
};

//static_assert(std::output_iterator<ParticleContainer::PairIterator, std::pair<Particle &, Particle &>>);
