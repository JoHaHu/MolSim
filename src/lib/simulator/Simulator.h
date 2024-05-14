#pragma once

#include "lib/Particle.h"
#include "lib/ParticleContainer.h"
#include "lib/utils/ArrayUtils.h"
#include <cmath>

//! Function for force calculation
/*!
  calculates the force for all particles, takes no arguments and has no return value
  */
auto calculateF(Particle const &p1, Particle const &p2) -> std::array<double, 3>;

//! Function for position calculation
/*!
  calculates the position for all particles, takes no arguments and has no return value
*/
auto calculateX(Particle const &p1, Particle const &p2) -> std::array<double, 3>;

//! Function for velocity calculation
/*!
  calculates the velocity for all particles, takes no arguments and has no return value
*/
auto calculateV(Particle const &p1, Particle const &p2) -> std::array<double, 3>;

//! Function for plotting particles
/*!
  plots the particles of the particle array, takes an integer value and has no return value
  \param iteration an integer argument that sets the number of iterations
*/
auto plotParticles(ParticleContainer &pc, int iteration) -> void;

template<typename T>
concept Plotter = requires(ParticleContainer &pc, int iteration) {
  { plotParticles(pc, iteration) };
};

template<typename T>
concept Physics = requires(Particle const &p1, Particle const &p2) {
  { calculateF(p1, p2) } -> std::convertible_to<std::array<double, 3>>;
};

template<Physics PY, Plotter PL>
class Simulator {

 private:
  ParticleContainer particles;
  PL plotter;
  PY physics;
  double start_time = 0;
  double end_time = 0;
  double delta_t = 0.014;

 public:
  /***
   * \param particles
   * \param plotter
   * \param physics
   * \param start_time  initialisation of the start time of the simulation with 0
   * \param end_time  initialisation of the end time of the simulation (default 1000)
   * \param delta_t  initialisation of time delta (defaul 0.014)
   * */
  explicit Simulator(
      ParticleContainer particles,
      PL const &plotter,
      PY const &physics,
      auto start_time,
      auto end_time,
      auto delta_t);

  auto run() -> void;
  auto calculateX() -> void;
  auto calculateV() -> void;
  auto calculateF() -> void;
};

template<Physics PY, Plotter PL>
Simulator<PY, PL>::Simulator(
    ParticleContainer particles, PL const &plotter, PY const &physics, auto start_time, auto end_time, auto delta_t)
    : particles(std::move(particles)), plotter(std::move(plotter)), physics(std::move(physics)), start_time(start_time), end_time(end_time), delta_t(delta_t) {
}

template<Physics PY, Plotter PL>
auto Simulator<PY, PL>::run() -> void {
  double current_time = start_time;
  int iteration = 0;

  while (current_time < end_time) {
    // calculate new x
    calculateX();
    // calculate new f
    calculateF();
    // calculate new v
    calculateV();

    iteration++;
    if (iteration % 10 == 0) {
      plotter.plotParticles(particles, iteration);
    }

    std::cout << "Iteration " << iteration << " finished." << std::endl;

    current_time += delta_t;
  }
  std::cout << "Output written. Terminating..." << std::endl;
}
template<Physics PY, Plotter PL>
void Simulator<PY, PL>::calculateX() {
  for (auto &p : particles) {
    p.x = p.x + delta_t * p.v + pow(delta_t, 2) * (1 / (2 * p.m)) * p.old_f;
  }
}

template<Physics PY, Plotter PL>
void Simulator<PY, PL>::calculateV() {
  for (auto &p : particles) {
    p.v = p.v + delta_t * (1 / (2 * p.m)) * (p.old_f + p.f);
  }
}

template<Physics PY, Plotter PL>
void Simulator<PY, PL>::calculateF() {
  for (auto &p : particles) {
    p.old_f = p.f;
    p.f = {0, 0, 0};
  }

  for (auto pair = particles.begin_pair(); pair != particles.end_pair(); pair++) {
    const auto [p1, p2] = *pair;

    const auto f = physics.calculateF(p1, p2);

    p1.f = p1.f + f;
    p2.f = p2.f - f;
  }
}