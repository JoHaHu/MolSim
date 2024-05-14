#pragma once

#include "lib/Particle.h"
#include "lib/ParticleContainer.h"
#include "lib/utils/ArrayUtils.h"
#include <cmath>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>

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
   * \param start_time
   * \param end_time
   * \param delta_t
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
  spdlog::info("Running simulation...");
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
      spdlog::debug("Iteration {} plotted.", iteration);
    }

    spdlog::debug("Iteration {} finished.", iteration);

    current_time += delta_t;
  }
  spdlog::info("Output written. Terminating...");
}

template<Physics PY, Plotter PL>
void Simulator<PY, PL>::calculateX() {
  spdlog::debug("Updating positions for {} particles.", particles.size());
  for (auto &p : particles) {
    p.x = p.x + delta_t * p.v + pow(delta_t, 2) * (1 / (2 * p.m)) * p.old_f;
    spdlog::trace("Particle position updated: ({}, {}, {})", p.x[0], p.x[1], p.x[2]);
  }
}

template<Physics PY, Plotter PL>
void Simulator<PY, PL>::calculateV() {
  spdlog::debug("Updating velocities for {} particles.", particles.size());
  for (auto &p : particles) {
    p.v = p.v + delta_t * (1 / (2 * p.m)) * (p.old_f + p.f);
    spdlog::trace("Particle velocity updated: ({}, {}, {})", p.v[0], p.v[1], p.v[2]);
  }
}

template<Physics PY, Plotter PL>
void Simulator<PY, PL>::calculateF() {
  spdlog::debug("Starting force calculation for {} particles.", particles.size());
  for (auto &p : particles) {
    p.old_f = p.f;
    p.f = {0, 0, 0};
  }

  for (auto pair = particles.begin_pair(); pair != particles.end_pair(); pair++) {
    const auto [p1, p2] = *pair;

    const auto f = physics.calculateF(p1, p2);

    p1.f = p1.f + f;
    p2.f = p2.f - f;

    spdlog::trace("Force updated for particle pair: ({}, {}, {}) - ({}, {}, {})", p1.f[0], p1.f[1], p1.f[2], p2.f[0], p2.f[1], p2.f[2]);
  }
  spdlog::trace("Force calculation completed.");
}
