#pragma once

#include "ParticleContainer.h"
#include "config/Config.h"
#include "simulator/physics/ForceModel.h"
#include <memory>
#include <optional>

#include "config/Config.h"

namespace simulator::io {

/*!
 * A recursive descent parser to load scenario from file and populate ParticleContainer
 */
class ParticleGenerator {
 private:
  std::shared_ptr<config::Config> config;

  using cuboid_t = std::tuple<std::array<double, 3>, std::array<double, 3>, std::array<int, 3>, double, double, double>;

 public:
  explicit ParticleGenerator(const std::shared_ptr<config::Config> &config);
  /**
   * Load particles based on a input file and returns a particle container and the used force model
   * */
  auto load_particles() -> std::vector<Particle>;

  auto generate_cuboids(const std::vector<Cuboid> &, auto seed, std::vector<Particle> &particles);

  std::vector<Particle> generate_disk_particles(double centerX, double centerY, double initialVx, double initialVy,
                                                int radiusMolecules, double meshwidth, unsigned int seed);

  std::vector<Particle> generate_sphere_particles(double centerX, double centerY, double centerZ, double initialVx,
                                                  double initialVy, double initialVz, int radiusMolecules,
                                                  double meshwidth, unsigned int seed);

  std::vector<Particle> generate_torus_particles(double centerX, double centerY, double centerZ, double initialVx,
                                                 double initialVy, double initialVz, double majorRadius,
                                                 double minorRadius, double meshwidth, unsigned int seed);

  std::vector<Particle> generate_double_helix_particles(double centerX, double centerY, double centerZ,
                                                        double initialVx,
                                                        double initialVy, double initialVz, double helixRadius,
                                                        double helixPitch, double helixHeight, double meshwidth,
                                                        unsigned int seed);
};

}// namespace simulator::io
