
#include <fstream>
#include <ranges>
#include <sstream>

#include "ParticleLoader.h"
#include "lib/Particle.h"
#include "lib/utils/ArrayUtils.h"
#include "lib/utils/LoggerManager.h"
#include "lib/utils/MaxwellBoltzmannDistribution.h"
#include "spdlog/spdlog.h"

namespace simulator::io {
// TODO evaluate operator>> to primitive types
ParticleLoader::ParticleLoader(const std::shared_ptr<config::Config> &config) : config(config) {}

auto ParticleLoader::load_particles() -> std::tuple<ParticleContainer, simulator::physics::ForceModel> {
  spdlog::info("Start loading particles");
  auto input = std::ifstream(config->input_filename);
  auto input_buf = std::istreambuf_iterator<char>(input);

  while (recognize_comment(input_buf)) {}

  auto [model, length] = *recognize_header(input_buf);

  ParticleContainer particles(0);
  switch (model) {
    case physics::ForceModel::LennardJones: {
      auto cuboids = *parse_cuboids(input_buf);
      auto p = generate_cuboids(cuboids, config->seed);
      particles = ParticleContainer(p);
      break;
    }
    case physics::ForceModel::Gravity: {
      auto p = *parse_gravity(input_buf);
      particles = ParticleContainer(p);
      break;
    }
  }
  spdlog::info("particles loaded");
  return {particles, model};
}
auto ParticleLoader::recognize_force_model(std::istreambuf_iterator<char> &buf) -> std::optional<physics::ForceModel> {
  recognize_whitespace(buf);

  auto tmp = std::ostringstream();
  while ((isalpha(*buf) != 0) || *buf == '-') {
    tmp << *buf;
    buf++;
  }
  if (tmp.str() == "lennard-jones") {
    return {physics::ForceModel::LennardJones};
  }
  if (tmp.str() == "gravity") {
    return {physics::ForceModel::Gravity};
  }
  return {};
}
auto ParticleLoader::recognize_comment(std::istreambuf_iterator<char> &buf) -> bool {
  if (*buf == '#') {
    while (*(++buf) != '\n') {}
    ++buf;
    return true;
  }
  return false;
}
auto ParticleLoader::recognize_whitespace(std::istreambuf_iterator<char> &buf) -> bool {
  if (*buf == ' ') {
    while (*(++buf) == ' ') {}
    return true;
  }
  return false;
}

auto ParticleLoader::recognize_header(std::istreambuf_iterator<char> &buf) -> std::optional<std::tuple<physics::ForceModel, int>> {
  recognize_whitespace(buf);
  auto force = *recognize_force_model(buf);
  recognize_whitespace(buf);
  auto length = *recognize_int(buf);
  recognize_end_of_line(buf);
  return {{force, length}};
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
auto ParticleLoader::recognize_double(std::istreambuf_iterator<char> &buf) -> std::optional<double> {
  recognize_whitespace(buf);
  auto tmp = std::ostringstream();
  while ((isdigit(*buf) != 0) || (*buf == '.') || *buf == '-' || *buf == 'e') {
    tmp << *buf;
    buf++;
  }
  return std::stod(tmp.str());
}
auto ParticleLoader::recognize_double_triplet(std::istreambuf_iterator<char> &buf) -> std::optional<std::array<double, 3>> {
  recognize_whitespace(buf);
  auto first = *recognize_double(buf);
  recognize_whitespace(buf);
  auto second = *recognize_double(buf);
  recognize_whitespace(buf);
  auto third = *recognize_double(buf);

  return {{first, second, third}};
}
auto ParticleLoader::recognize_dimension_triplet(std::istreambuf_iterator<char> &buf) -> std::optional<std::array<int, 3>> {
  recognize_whitespace(buf);
  auto first = *recognize_int(buf);
  recognize_whitespace(buf);
  auto second = *recognize_int(buf);
  recognize_whitespace(buf);
  auto third = *recognize_int(buf);
  return {{first, second, third}};
}
auto ParticleLoader::recognize_cuboid(std::istreambuf_iterator<char> &buf) -> std::optional<cuboid_t> {

  recognize_whitespace(buf);
  auto position = *recognize_double_triplet(buf);
  auto velocity = *recognize_double_triplet(buf);
  auto dim = *recognize_dimension_triplet(buf);
  auto h = *recognize_double(buf);
  auto mass = *recognize_double(buf);
  auto sigma = *recognize_double(buf);
  recognize_end_of_line(buf);
  return {{position, velocity, dim, h, mass, sigma}};
}
auto ParticleLoader::recognize_planet(std::istreambuf_iterator<char> &buf) -> std::optional<Particle> {
  recognize_whitespace(buf);
  auto position = *recognize_double_triplet(buf);
  auto velocity = *recognize_double_triplet(buf);
  auto mass = *recognize_double(buf);
  auto p = Particle(position, velocity, mass, 0);
  recognize_end_of_line(buf);
  return {p};
}
auto ParticleLoader::parse_gravity(std::istreambuf_iterator<char> &buf) -> std::optional<std::vector<Particle>> {

  auto particles = std::vector<Particle>();
  while (true) {
    if (recognize_end_of_file(buf)) {
      break;
    }
    if (auto particle = recognize_planet(buf)) {
      particles.emplace_back(*particle);
    } else {
      break;
    }
  }
  return particles;
}
auto ParticleLoader::parse_cuboids(std::istreambuf_iterator<char> &buf) -> std::optional<std::vector<cuboid_t>> {
  auto particles = std::vector<cuboid_t>();
  while (true) {
    if (recognize_end_of_file(buf)) {
      break;
    }
    if (auto particle = recognize_cuboid(buf)) {
      particles.emplace_back(*particle);
    } else {
      break;
    }
  }
  return particles;
}
static constexpr const double brownian_motion = 0.1;
auto ParticleLoader::generate_cuboids(const std::vector<cuboid_t> &cuboids, auto seed) -> std::vector<Particle> {

  auto particles = std::vector<Particle>();

  for (auto const [index, cuboid] : std::views::enumerate(cuboids)) {
    auto [position, velocity, dim, h, m, sigma] = cuboid;
    for (auto x : std::views::iota(0, dim[0])) {
      for (auto y : std::views::iota(0, dim[1])) {
        for (auto z : std::views::iota(0, dim[2])) {
          auto particle = Particle(position + std::array<double, 3>({h * static_cast<double>(x), h * static_cast<double>(y), h * static_cast<double>(z)}), velocity + maxwellBoltzmannDistributedVelocity(brownian_motion, 2, seed), m, static_cast<int>(index));
          particles.emplace_back(particle);
        }
      }
    }
  }

  particles.shrink_to_fit();
  return particles;
}
auto ParticleLoader::recognize_end_of_line(std::istreambuf_iterator<char> &buf) -> bool {
  recognize_whitespace(buf);
  if (*buf == '\n') {
    buf++;
    return true;
  }
  return false;
}
auto ParticleLoader::recognize_end_of_file(std::istreambuf_iterator<char> &buf) -> bool {
  recognize_whitespace(buf);
  if (*buf == '\xff') {
    return true;
  }
  return false;
}

}// namespace simulator::io