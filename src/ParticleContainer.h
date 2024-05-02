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
    Iterator(
        std::vector<Particle>::iterator i1,
        std::vector<Particle>::iterator i2,
        std::vector<Particle>::iterator end
    );

    using difference_type = int;
    std::pair<Particle &, Particle &> operator*();

    Iterator &operator++();
    Iterator operator++(int);
    bool operator==(Iterator &);
   private:
    std::vector<Particle>::iterator i1;
    std::vector<Particle>::iterator i2;
    std::vector<Particle>::iterator end;
  };

  Iterator begin();
  Iterator end();

 private:
  std::vector<Particle> particles;
};

static_assert(std::output_iterator<ParticleContainer::Iterator, std::pair<Particle &, Particle &>>);
