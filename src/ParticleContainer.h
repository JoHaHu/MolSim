#include "Particle.h"
#include <vector>

class ParticleContainer {
  /***
   * @param capacity Indicates the capacity of the underlying backing storage
   * */
  explicit ParticleContainer(int capacity);

 public:

  class Iterator {
   public:
    explicit Iterator(std::vector<Particle>::iterator iterator);
    using difference_type = int;
    std::pair<Particle &, Particle &> operator*() const;

    Iterator &operator++() const;
    Iterator operator++(int) const;
   private:
    std::vector<Particle>::iterator i1;
    std::vector<Particle>::iterator i2;

  };

  Iterator begin();
  Iterator end();

 private:
  std::vector<Particle> particles;
};

static_assert(std::output_iterator<ParticleContainer::Iterator, std::pair<Particle &, Particle &>>);
