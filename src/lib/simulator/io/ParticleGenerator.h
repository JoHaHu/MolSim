#pragma once

#include "config/Config.h"
#include "range/v3/view/enumerate.hpp"
#include "range/v3/view/iota.hpp"
#include "simulator/physics/ForceModel.h"
#include <memory>
#include <optional>

#include <fstream>
#include <sstream>

#include "Particle.h"
#include "spdlog/spdlog.h"
#include "utils/ArrayUtils.h"
namespace simulator::io {

/*!
 * A recursive descent parser to load scenario from file and populate ParticleContainer
 */

template<const size_t DIMENSIONS>
class ParticleGenerator {
 private:
  std::shared_ptr<config::Config> config;

  using cuboid_t = std::tuple<std::array<double, DIMENSIONS>, std::array<double, DIMENSIONS>, std::array<int, DIMENSIONS>, double, double, double>;

 public:
  explicit ParticleGenerator(const std::shared_ptr<config::Config> &config) : config(config) {};
  /**
   * Load particles based on a input file and returns a particle container and the used force model
   * */
  auto load_particles() -> std::vector<Particle<DIMENSIONS>> {
    auto particles = std::vector<Particle<DIMENSIONS>>();
    generate_cuboids(config->cuboids, particles);
    generate_disk_particles(config->spheres, particles);
    return particles;
  };

  void _generate_cuboids(size_t depth, std::array<double, DIMENSIONS> &offset, std::vector<Particle<DIMENSIONS>> &particles, const Cuboid &cuboid, int type) {
    if (depth == 0) {
      std::array<double, DIMENSIONS> new_coords;
      for (size_t i = 0; i < DIMENSIONS; ++i) {
        new_coords[i] = cuboid.coordinates[i] + cuboid.spacing * static_cast<double>(offset[i]);
      }

      std::array<double, DIMENSIONS> velocity;
      for (size_t i = 0; i < DIMENSIONS; ++i) {
        velocity[i] = cuboid.velocity[i];
      }
      int fixed = cuboid.fixed;
      const auto particle = Particle<DIMENSIONS>(
          new_coords, velocity, config->mass[type],
          type, fixed);
      particles.push_back(particle);
    } else {
      for (size_t z : ranges::iota_view(0, cuboid.particles[depth - 1])) {
        offset[depth - 1] = z;
        _generate_cuboids(depth - 1, offset, particles, cuboid, type);
      }
    }
  }

  /**
    * @brief Generates particles from a list of cuboids using the given seed.
    *
    * Iterates through each cuboid and generates particles based on its dimensions and properties.
    * Logs the progress at various levels.
    *
    * @param cuboids A vector of cuboid_t structures containing cuboid properties.
    * @param seed The seed used for random number generation.
    * @return std::vector<Particle> The generated particles.
    */
  auto generate_cuboids(std::vector<Cuboid> &cuboids, std::vector<Particle<DIMENSIONS>> &particles) {
    for (auto cuboid : cuboids) {
      for (double pos : cuboid.coordinates) {
        if (std::isnan(pos)) {
          SPDLOG_ERROR("Invalid cuboid data: NaN value encountered in position.");
          continue;// Skip this cuboid
        }
      }
      for (double vel : cuboid.velocity) {
        if (std::isnan(vel)) {
          SPDLOG_ERROR("Invalid cuboid data: NaN value encountered in velocity.");
          continue;// Skip this cuboid
        }
      }
      std::array<double, DIMENSIONS> temp;
      _generate_cuboids(DIMENSIONS, temp, particles, cuboid, cuboid.type);
    }

    particles.shrink_to_fit();
    return particles;
  };
  /**
            * @brief Generates particles arranged in a disk configuration.
            *
            * This method generates particles within a disk of a specified radius. The disk
            * is centered at (centerX, centerY) and each particle is given an initial velocity
            * (initialVx, initialVy). The density of the particles is controlled by the meshwidth,
            * which determines the spacing between adjacent particles. The number of particles
            * along the radius is specified by radiusMolecules.
            *
            * @param centerX The x-coordinate of the center of the disk.
            * @param centerY The y-coordinate of the center of the disk.
            * @param initialVx The initial velocity in the x-direction for all particles.
            * @param initialVy The initial velocity in the y-direction for all particles.
            * @param radiusMolecules The radius of the disk in terms of the number of molecules.
            * @param meshwidth The reflecting_distance between adjacent particles.
            * @return std::vector<Particle> A vector containing the generated particles.
            */

  auto generate_disk_particles(std::vector<Sphere> &spheres, std::vector<Particle<DIMENSIONS>> &particles) {

    for (auto &sphere : spheres) {

      double radius = sphere.radius * sphere.mesh_width;

      for (int i = -sphere.radius; i <= sphere.radius; ++i) {
        for (int j = -sphere.radius; j <= sphere.radius; ++j) {
          double x = i * sphere.mesh_width;
          double y = j * sphere.mesh_width;
          if (x * x + y * y <= radius * radius) {

            std::array<double, DIMENSIONS> position;
            if constexpr (DIMENSIONS == 2) {
              position = {sphere.coordinates[0] + x, sphere.coordinates[1] + y};
            } else {
              position = {sphere.coordinates[0] + x, sphere.coordinates[1] + y, 0.0};
            }

            std::array<double, DIMENSIONS> velocity;

            for (size_t k = 0; k < DIMENSIONS; ++k) {
              velocity[k] = 0.0;
            }

            particles.emplace_back(position, velocity, 1.0, 0);
          }
        }
      }
    }

    particles.shrink_to_fit();
    return particles;
  };
  /**
   * @brief Generates particles arranged in a sphere configuration.
   *
   * This method generates particles within a sphere of a specified radius. The sphere
   * is centered at (centerX, centerY, centerZ) and each particle is given an initial velocity
   * (initialVx, initialVy, initialVz). The density of the particles is controlled by the meshwidth,
   * which determines the spacing between adjacent particles. The number of particles
   * along the radius is specified by radiusMolecules.
   *
   * @param centerX The x-coordinate of the center of the sphere.
   * @param centerY The y-coordinate of the center of the sphere.
   * @param centerZ The z-coordinate of the center of the sphere.
   * @param initialVx The initial velocity in the x-direction for all particles.
   * @param initialVy The initial velocity in the y-direction for all particles.
   * @param initialVz The initial velocity in the z-direction for all particles.
   * @param radiusMolecules The radius of the sphere in terms of the number of molecules.
   * @param meshwidth The reflecting_distance between adjacent particles.
   * @return std::vector<Particle> A vector containing the generated particles.
   */
  auto generate_sphere_particles(double centerX, double centerY, double centerZ, double initialVx,
                                 double initialVy, double initialVz, int radiusMolecules,
                                 double meshwidth) -> std::vector<Particle<DIMENSIONS>> {

    std::vector<Particle<DIMENSIONS>> particles;

    if constexpr (DIMENSIONS >= 3) {
      double radius = radiusMolecules * meshwidth;

      for (int i = -radiusMolecules; i <= radiusMolecules; ++i) {
        for (int j = -radiusMolecules; j <= radiusMolecules; ++j) {
          for (int k = -radiusMolecules; k <= radiusMolecules; ++k) {
            double x = i * meshwidth;
            double y = j * meshwidth;
            double z = k * meshwidth;
            if (x * x + y * y + z * z <= radius * radius) {
              std::array<double, 3> position = {centerX + x, centerY + y, centerZ + z};
              std::array<double, 3> velocity = {initialVx, initialVy, initialVz};
              particles.emplace_back(position, velocity, 1.0, 0);
            }
          }
        }
      }
    } else {
      SPDLOG_WARN("Spheres need 3 dimensions");
    }
    particles.shrink_to_fit();
    return particles;
  };

  /**
 * @brief Generates particles arranged in a torus configuration.
 *
 * This method generates particles within a torus of specified major and minor radii. The torus
 * is centered at (centerX, centerY, centerZ) and each particle is given an initial velocity
 * (initialVx, initialVy, initialVz). The density of the particles is controlled by the meshwidth,
 * which determines the spacing between adjacent particles.
 *
 * @param centerX The x-coordinate of the center of the torus.
 * @param centerY The y-coordinate of the center of the torus.
 * @param centerZ The z-coordinate of the center of the torus.
 * @param initialVx The initial velocity in the x-direction for all particles.
 * @param initialVy The initial velocity in the y-direction for all particles.
 * @param initialVz The initial velocity in the z-direction for all particles.
 * @param majorRadius The major radius of the torus (reflecting_distance from the center of the tube to the center of the torus).
 * @param minorRadius The minor radius of the torus (radius of the tube).
 * @param meshwidth The reflecting_distance between adjacent particles.
 * @return std::vector<Particle> A vector containing the generated particles.
 */
  auto generate_torus_particles(double centerX, double centerY, double centerZ, double initialVx,
                                double initialVy, double initialVz, double majorRadius,
                                double minorRadius, double meshwidth) -> std::vector<Particle<DIMENSIONS>> {
    std::vector<Particle<DIMENSIONS>> particles;

    if constexpr (DIMENSIONS == 3) {

      // Loop over the angular coordinates
      for (double theta = 0; theta < 2 * M_PI; theta += 4 * meshwidth / majorRadius) {
        for (double phi = 0; phi < 2 * M_PI; phi += meshwidth / minorRadius) {
          double x = (majorRadius + minorRadius * cos(phi)) * cos(theta);
          double y = (majorRadius + minorRadius * cos(phi)) * sin(theta);
          double z = minorRadius * sin(phi);

          std::array<double, 3> position = {centerX + x, centerY + y, centerZ + z};
          std::array<double, 3> velocity{};
          velocity[0] = initialVx;
          velocity[1] = initialVy;
          velocity[2] = initialVz;

          particles.emplace_back(position, velocity, 1.0, 0);
        }
      }
      particles.shrink_to_fit();
    } else {
      SPDLOG_WARN("Torus needs 3 Dimensions");
    }
    return particles;
  };

  /**
 * @brief Generates particles arranged in a double helix configuration.
 *
 * This method generates particles within a double helix with specified radius, pitch, and height. The double helix
 * is centered at (centerX, centerY, centerZ) and each particle is given an initial velocity
 * (initialVx, initialVy, initialVz). The density of the particles is controlled by the meshwidth,
 * which determines the spacing between adjacent particles along the helix.
 *
 * @param centerX The x-coordinate of the center of the double helix.
 * @param centerY The y-coordinate of the center of the double helix.
 * @param centerZ The z-coordinate of the center of the double helix.
 * @param initialVx The initial velocity in the x-direction for all particles.
 * @param initialVy The initial velocity in the y-direction for all particles.
 * @param initialVz The initial velocity in the z-direction for all particles.
 * @param helixRadius The radius of the helix (reflecting_distance from the center of the helix to the path of the particles).
 * @param helixPitch The pitch of the helix (reflecting_distance between consecutive turns of the helix along the z-axis).
 * @param helixHeight The total height of the helix (reflecting_distance covered along the z-axis).
 * @param meshwidth The reflecting_distance between adjacent particles along the helix path.
 * @param seed The seed for random number generation, used in velocity distribution.
 * @return std::vector<Particle> A vector containing the generated particles.
 */
  auto generate_double_helix_particles(double centerX, double centerY, double centerZ,
                                       double initialVx,
                                       double initialVy, double initialVz, double helixRadius,
                                       double helixPitch, double helixHeight, double meshwidth) -> std::vector<Particle<DIMENSIONS>> {
    if (helixRadius <= 0 || helixPitch <= 0 || helixHeight <= 0 || meshwidth <= 0) {
      SPDLOG_ERROR("Invalid parameters for double helix generation.");
      return {};
    }

    std::vector<Particle<DIMENSIONS>> particles;

    if constexpr (DIMENSIONS == 3) {

      double numTurns = helixHeight / helixPitch;
      double thetaIncrement = 12 * meshwidth / helixRadius;

      // First helix
      for (double theta = 0; theta < numTurns * 2 * M_PI; theta += thetaIncrement) {
        double x = helixRadius * cos(theta);
        double y = helixRadius * sin(theta);
        double z = helixPitch * theta / (2 * M_PI);

        std::array<double, 3> position = {centerX + x, centerY + y, centerZ + z};
        std::array<double, 3> velocity{};
        velocity[0] += initialVx;
        velocity[1] += initialVy;
        velocity[2] += initialVz + 1;

        particles.emplace_back(position, velocity, 1.0, particles.size());
      }

      // Second helix
      for (double theta = 0; theta < numTurns * 2 * M_PI; theta += thetaIncrement) {
        double x = helixRadius * cos(theta + M_PI);
        double y = helixRadius * sin(theta + M_PI);
        double z = helixPitch * theta / (2 * M_PI);

        std::array<double, 3> position = {centerX + x, centerY + y, centerZ + z};
        std::array<double, 3> velocity{};
        velocity[0] += initialVx;
        velocity[1] += initialVy;
        velocity[2] += initialVz - 1;

        particles.emplace_back(position, velocity, 1.0, particles.size());
      }

      particles.shrink_to_fit();
    } else {
      SPDLOG_WARN("Helix need 3 Dimensions");
    }
    return particles;
  };
};

}// namespace simulator::io
