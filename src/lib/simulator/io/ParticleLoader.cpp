#include <fstream>
#include <ranges>
#include <sstream>

#include "Particle.h"
#include "ParticleLoader.h"
#include "config/Config.h"
#include "spdlog/spdlog.h"
#include "utils/ArrayUtils.h"
#include "utils/LoggerManager.h"
#include "utils/MaxwellBoltzmannDistribution.h"

namespace simulator::io {
static constexpr const double brownian_motion = 0.1;

ParticleLoader::ParticleLoader(const std::shared_ptr<config::Config> &config) : config(config) {
}

/**
    * @brief Loads particles from the input file.
    *
    * Reads the input file specified in the config, recognizes comments, headers, and force models,
    * and generates particles accordingly. Logs the process at various levels.
    *
    * @return std::tuple<ParticleContainer, simulator::physics::ForceModel> The loaded particles and the force model.
    */
auto ParticleLoader::load_particles() -> std::tuple<ParticleContainer, ForceModel> {
  spdlog::info("Reading file {}", config->input_filename);
  auto input = std::ifstream(config->input_filename);
  auto input_buf = std::istreambuf_iterator<char>(input);

  while (recognize_comment(input_buf)) {
    spdlog::trace("Recognized a comment in the input file.");
  }

  const auto model = config->simulation_type;
  spdlog::debug("Recognized header: Model = {}", static_cast<int>(model));

  ParticleContainer particles(0);
  switch (model) {
    case ForceModel::LennardJones: {
      spdlog::info("Processing Lennard-Jones model.");
      const auto cuboids = *parse_cuboids(input_buf);
      spdlog::debug("Parsed {} cuboids.", cuboids.size());
      const auto seed = config->seed;
      const auto p = generate_cuboids(cuboids, seed);
      spdlog::debug("Seed: ", seed);
      spdlog::debug("Generated {} particles for Lennard-Jones model.", p.size());
      particles = ParticleContainer(p);
      break;
    }
    case ForceModel::Gravity: {
      spdlog::info("Processing Gravity model.");
      const auto p = *parse_gravity(input_buf);
      spdlog::debug("Parsed {} particles for Gravity model.", p.size());
      particles = ParticleContainer(p);
      break;
    }
  }
  spdlog::debug("File {} read successfully.", config->input_filename);
  return {particles, model};
}

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
            * @param meshwidth The distance between adjacent particles.
            * @return std::vector<Particle> A vector containing the generated particles.
            */
std::vector<Particle> ParticleLoader::generate_disk_particles(double centerX, double centerY, double initialVx,
                                                              double initialVy, int radiusMolecules,
                                                              double meshwidth, unsigned int seed) {
  std::vector<Particle> particles;
  double radius = radiusMolecules * meshwidth;

  for (int i = -radiusMolecules; i <= radiusMolecules; ++i) {
    for (int j = -radiusMolecules; j <= radiusMolecules; ++j) {
      double x = i * meshwidth;
      double y = j * meshwidth;
      if (x * x + y * y <= radius * radius) {
        std::array<double, 3> position = {centerX + x, centerY + y, 0.0};
        std::array<double, 3> velocity = {0.0, 0.0, 0.0};

        std::array<double, 3> velocity2D = maxwellBoltzmannDistributedVelocity(brownian_motion, 2, seed);
        velocity[0] = velocity2D[0] + initialVx;
        velocity[1] = velocity2D[1] + initialVy;
        velocity[2] = 0.0;

        particles.emplace_back(position, velocity, 1.0, 0);
      }
    }
  }

  particles.shrink_to_fit();
  return particles;
}

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
   * @param meshwidth The distance between adjacent particles.
   * @return std::vector<Particle> A vector containing the generated particles.
   */
std::vector<Particle> ParticleLoader::generate_sphere_particles(double centerX, double centerY, double centerZ,
                                                                double initialVx, double initialVy,
                                                                double initialVz,
                                                                int radiusMolecules, double meshwidth,
                                                                unsigned int seed) {
  std::vector<Particle> particles;
  double radius = radiusMolecules * meshwidth;

  for (int i = -radiusMolecules; i <= radiusMolecules; ++i) {
    for (int j = -radiusMolecules; j <= radiusMolecules; ++j) {
      for (int k = -radiusMolecules; k <= radiusMolecules; ++k) {
        double x = i * meshwidth;
        double y = j * meshwidth;
        double z = k * meshwidth;
        if (x * x + y * y + z * z <= radius * radius) {
          std::array<double, 3> position = {centerX + x, centerY + y, centerZ + z};
          std::array<double, 3> velocity = maxwellBoltzmannDistributedVelocity(brownian_motion, 3, seed);
          velocity[0] += initialVx;
          velocity[1] += initialVy;
          velocity[2] += initialVz;

          particles.emplace_back(position, velocity, 1.0, 0);
        }
      }
    }
  }

  particles.shrink_to_fit();
  return particles;
}

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
 * @param majorRadius The major radius of the torus (distance from the center of the tube to the center of the torus).
 * @param minorRadius The minor radius of the torus (radius of the tube).
 * @param meshwidth The distance between adjacent particles.
 * @return std::vector<Particle> A vector containing the generated particles.
 */
std::vector<Particle> ParticleLoader::generate_torus_particles(double centerX, double centerY, double centerZ,
                                                               double initialVx, double initialVy, double initialVz,
                                                               double majorRadius, double minorRadius,
                                                               double meshwidth,
                                                               unsigned int seed) {
  std::vector<Particle> particles;

  // Loop over the angular coordinates
  for (double theta = 0; theta < 2 * M_PI; theta += 4 * meshwidth / majorRadius) {
    for (double phi = 0; phi < 2 * M_PI; phi += meshwidth / minorRadius) {
      double x = (majorRadius + minorRadius * cos(phi)) * cos(theta);
      double y = (majorRadius + minorRadius * cos(phi)) * sin(theta);
      double z = minorRadius * sin(phi);

      std::array<double, 3> position = {centerX + x, centerY + y, centerZ + z};
      std::array<double, 3> velocity = maxwellBoltzmannDistributedVelocity(brownian_motion, 3, seed);
      velocity[0] += initialVx;
      velocity[1] += initialVy;
      velocity[2] += initialVz;

      particles.emplace_back(position, velocity, 1.0, 0);
    }
  }
  particles.shrink_to_fit();
  return particles;
}

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
 * @param helixRadius The radius of the helix (distance from the center of the helix to the path of the particles).
 * @param helixPitch The pitch of the helix (distance between consecutive turns of the helix along the z-axis).
 * @param helixHeight The total height of the helix (distance covered along the z-axis).
 * @param meshwidth The distance between adjacent particles along the helix path.
 * @param seed The seed for random number generation, used in velocity distribution.
 * @return std::vector<Particle> A vector containing the generated particles.
 */
std::vector<Particle> ParticleLoader::generate_double_helix_particles(
    double centerX, double centerY, double centerZ,
    double initialVx, double initialVy, double initialVz,
    double helixRadius, double helixPitch, double helixHeight,
    double meshwidth, unsigned int seed) {

  if (helixRadius <= 0 || helixPitch <= 0 || helixHeight <= 0 || meshwidth <= 0) {
    spdlog::error("Invalid parameters for double helix generation.");
    return {};
  }

  std::vector<Particle> particles;

  double numTurns = helixHeight / helixPitch;
  double thetaIncrement = 12 * meshwidth / helixRadius;

  // First helix
  for (double theta = 0; theta < numTurns * 2 * M_PI; theta += thetaIncrement) {
    double x = helixRadius * cos(theta);
    double y = helixRadius * sin(theta);
    double z = helixPitch * theta / (2 * M_PI);

    std::array<double, 3> position = {centerX + x, centerY + y, centerZ + z};
    std::array<double, 3> velocity = maxwellBoltzmannDistributedVelocity(brownian_motion, 3, seed);
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
    std::array<double, 3> velocity = maxwellBoltzmannDistributedVelocity(brownian_motion, 3, seed);
    velocity[0] += initialVx;
    velocity[1] += initialVy;
    velocity[2] += initialVz - 1;

    particles.emplace_back(position, velocity, 1.0, particles.size());
  }

  particles.shrink_to_fit();
  return particles;
}

auto ParticleLoader::recognize_comment(std::istreambuf_iterator<char> &buf) -> bool {
  if (*buf == '#') {
    while (*(++buf) != '\n') {
    }
    ++buf;
    return true;
  }
  return false;
}

auto ParticleLoader::recognize_whitespace(std::istreambuf_iterator<char> &buf) -> bool {
  if (*buf == ' ') {
    while (*(++buf) == ' ') {
    }
    return true;
  }
  return false;
}

auto ParticleLoader::recognize_int(std::istreambuf_iterator<char> &buf) -> std::optional<int> {
  recognize_whitespace(buf);

  auto tmp = std::ostringstream();
  while (isdigit(*buf) != 0) {
    tmp << *buf;
    buf++;
  }
  return std::stoi(tmp.str());
}

/**
    * @brief Extracts and parses a double value from an input stream buffer.
    *
    * @param buf Iterator into the input stream buffer.
    * @return std::optional<double> The parsed double, or std::nullopt if:
    *        - end of buffer reached prematurely
    *        - parsed value is NaN or infinite
    *        - parsing error occurs
    */
auto ParticleLoader::recognize_double(std::istreambuf_iterator<char> &buf) -> std::optional<double> {
  recognize_whitespace(buf);
  if (buf == std::istreambuf_iterator<char>()) return {};

  std::string doubleStr;
  while (buf != std::istreambuf_iterator<char>() && (std::isdigit(*buf) || *buf == '.' || *buf == '-' || *buf == 'e' || *buf == 'E')) {
    doubleStr += *buf++;
  }

  try {
    const double parsedValue = std::stod(doubleStr);
    if (!std::isfinite(parsedValue)) {
      // Combined NaN and infinity check
      spdlog::error("Parsed value not finite: '{}'", doubleStr);
      return {};
    }
    return parsedValue;
  } catch (const std::invalid_argument &) {
    spdlog::error("Failed to parse double: '{}'", doubleStr);
    return {};
  }
}

/**
    * @brief Recognizes a triplet of doubles from the input buffer.
    *
    * @param buf Input buffer iterator.
    * @return std::optional<std::array<double, 3>> The recognized triplet or an empty optional if not recognized.
    */
auto ParticleLoader::recognize_double_triplet(std::istreambuf_iterator<char> &buf)
    -> std::optional<std::array<double, 3>> {
  spdlog::trace("Starting to recognize a double triplet");
  if (buf == std::istreambuf_iterator<char>()) return {};

  try {
    std::array<double, 3> result;
    for (auto i = 0; i < 3; ++i) {
      recognize_whitespace(buf);
      const auto value = recognize_double(buf);
      if (!value) {
        spdlog::error("Failed to recognize double #{}", i + 1);
        return {};
      }
      result[i] = *value;
    }

    spdlog::debug("Recognized double triplet: ({}, {}, {})", result[0], result[1], result[2]);
    return result;
  } catch (const std::exception &e) {
    spdlog::error("Error parsing double triplet: {}", e.what());
    return {};
  }
}

/**
     * @brief Recognizes a triplet of integers from the input buffer.
     *
     * @param buf Input buffer iterator.
     * @return std::optional<std::array<int, 3>> The recognized triplet or an empty optional if not recognized.
     */
auto ParticleLoader::recognize_dimension_triplet(
    std::istreambuf_iterator<char> &buf) -> std::optional<std::array<int, 3>> {
  spdlog::trace("Starting to recognize a dimension triplet");
  recognize_whitespace(buf);
  const auto first = *recognize_int(buf);
  recognize_whitespace(buf);
  const auto second = *recognize_int(buf);
  recognize_whitespace(buf);
  const auto third = *recognize_int(buf);
  spdlog::debug("Recognized dimension triplet: ({}, {}, {})", first, second, third);
  return {{first, second, third}};
}

/**
     * @brief Recognizes a cuboid structure from the input buffer.
     *
     * @param buf Input buffer iterator.
     * @return std::optional<cuboid_t> The recognized cuboid or an empty optional if not recognized.
     */
auto ParticleLoader::recognize_cuboid(std::istreambuf_iterator<char> &buf) -> std::optional<cuboid_t> {
  spdlog::trace("Starting to recognize a cuboid");
  recognize_whitespace(buf);
  const auto position = *recognize_double_triplet(buf);
  const auto velocity = *recognize_double_triplet(buf);
  const auto dim = *recognize_dimension_triplet(buf);
  const auto h = *recognize_double(buf);
  const auto mass = *recognize_double(buf);
  const auto sigma = *recognize_double(buf);
  recognize_end_of_line(buf);
  spdlog::debug(
      "Recognized cuboid: position = ({}, {}, {}), velocity = ({}, {}, {}), dimensions = ({}, {}, {}), h = {}, mass = {}, sigma = {}",
      position[0], position[1], position[2],
      velocity[0], velocity[1], velocity[2],
      dim[0], dim[1], dim[2], h, mass, sigma);
  return {{position, velocity, dim, h, mass, sigma}};
}

/**
     * @brief Recognizes a planet structure from the input buffer.
     *
     * @param buf Input buffer iterator.
     * @return std::optional<Particle> The recognized planet or an empty optional if not recognized.
     */
auto ParticleLoader::recognize_planet(std::istreambuf_iterator<char> &buf) -> std::optional<Particle> {
  spdlog::trace("Starting to recognize a planet");
  recognize_whitespace(buf);
  const auto position = *recognize_double_triplet(buf);
  const auto velocity = *recognize_double_triplet(buf);
  const auto mass = *recognize_double(buf);

  for (double d : position) {
    if (std::isnan(d)) {
      spdlog::error("Invalid planet data: NaN values encountered.");
      return std::nullopt;
    }
  }
  for (double d : velocity) {
    if (std::isnan(d)) {
      spdlog::error("Invalid planet data: NaN values encountered.");
      return std::nullopt;
    }
  }
  if (std::isnan(mass)) {
    spdlog::error("Invalid planet data: NaN mass encountered.");
    return std::nullopt;
  }

  const auto p = Particle(position, velocity, mass, 0);
  recognize_end_of_line(buf);
  spdlog::debug("Recognized planet: position = ({}, {}, {}), velocity = ({}, {}, {}), mass = {}",
                position[0], position[1], position[2],
                velocity[0], velocity[1], velocity[2], mass);
  return {p};
}

/**
     * @brief Parses particles affected by gravity from the input buffer.
     *
     * @param buf Input buffer iterator.
     * @return std::optional<std::vector<Particle>> The parsed particles or an empty optional if parsing fails.
     */
auto ParticleLoader::parse_gravity(std::istreambuf_iterator<char> &buf) -> std::optional<std::vector<Particle>> {
  spdlog::trace("Starting to parse gravity particles");
  auto particles = std::vector<Particle>();
  while (true) {
    if (recognize_end_of_file(buf)) {
      spdlog::trace("Reached end of file while parsing gravity particles");
      break;
    }
    if (const auto particle = recognize_planet(buf)) {
      particles.push_back(*particle);
    } else {
      spdlog::warn("Failed to recognize a planet while parsing gravity particles");
      break;
    }
  }
  spdlog::debug("Parsed {} gravity particles", particles.size());
  return particles;
}

/**
     * @brief Parses cuboid structures from the input buffer.
     *
     * @param buf Input buffer iterator.
     * @return std::optional<std::vector<cuboid_t>> The parsed cuboids or an empty optional if parsing fails.
     */
auto ParticleLoader::parse_cuboids(std::istreambuf_iterator<char> &buf) -> std::optional<std::vector<cuboid_t>> {
  spdlog::trace("Starting to parse cuboids");
  auto particles = std::vector<cuboid_t>();
  while (true) {
    if (recognize_end_of_file(buf)) {
      spdlog::trace("Reached end of file while parsing cuboids");
      break;
    }
    if (const auto particle = recognize_cuboid(buf)) {
      particles.push_back(*particle);
    } else {
      spdlog::warn("Failed to recognize a cuboid while parsing");
      break;
    }
  }
  spdlog::debug("Parsed {} cuboids", particles.size());
  return particles;
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
auto ParticleLoader::generate_cuboids(const std::vector<cuboid_t> &cuboids, auto seed) -> std::vector<Particle> {
  auto particles = std::vector<Particle>();

  for (auto const [index, cuboid] : std::views::enumerate(cuboids)) {
    const auto [position, velocity, dim, h, m, sigma] = cuboid;

    for (double pos : position) {
      if (std::isnan(pos)) {
        spdlog::error("Invalid cuboid data: NaN value encountered in position.");
        continue;// Skip this cuboid
      }
    }
    for (double vel : velocity) {
      if (std::isnan(vel)) {
        spdlog::error("Invalid cuboid data: NaN value encountered in velocity.");
        continue;// Skip this cuboid
      }
    }
    if (std::isnan(h) || std::isnan(m) || std::isnan(sigma)) {
      spdlog::error("Invalid cuboid data: NaN value encountered in h, m, or sigma.");
      continue;// Skip this cuboid
    }

    for (const auto x : std::views::iota(0, dim[0])) {
      for (const auto y : std::views::iota(0, dim[1])) {
        for (const auto z : std::views::iota(0, dim[2])) {
          const auto particle = Particle(
              position + std::array<double, 3>({h * static_cast<double>(x), h * static_cast<double>(y), h * static_cast<double>(z)}), velocity + maxwellBoltzmannDistributedVelocity(brownian_motion, 2, seed), m,
              static_cast<int>(index));
          particles.push_back(std::move(particle));
        }
      }
    }
  }

  particles.shrink_to_fit();
  return particles;
}

auto ParticleLoader::recognize_end_of_line(std::istreambuf_iterator<char> &buf) -> bool {
  recognize_whitespace(buf);
  bool result;
  while (*buf == '\n' || *buf == '\r') {
    buf++;
    result = true;
  }
  return result;
}

auto ParticleLoader::recognize_end_of_file(std::istreambuf_iterator<char> &buf) -> bool {
  recognize_whitespace(buf);
  if (*buf == '\xff') {
    return true;
  }
  return false;
}
}// namespace simulator::io
