#pragma once

#include "lib/ParticleContainer.h"
#include "lib/config/config.h"
#include "lib/simulator/physics/ForceModel.h"
#include <memory>
#include <optional>

namespace simulator::io {

/*!
 * A recursive descent parser to load scenario from file and populate ParticleContainer
 * */
class ParticleLoader {
 private:
  std::shared_ptr<config::Config> config;

  using cuboid_t = std::tuple<std::array<double, 3>, std::array<double, 3>, std::array<int, 3>, double, double, double>;

  auto generate_cuboids(const std::vector<cuboid_t> &, auto seed) -> std::vector<Particle>;
  static auto recognize_comment(std::istreambuf_iterator<char> &buf) -> bool;
  static auto recognize_whitespace(std::istreambuf_iterator<char> &buf) -> bool;
  static auto recognize_end_of_line(std::istreambuf_iterator<char> &buf) -> bool;
  static auto recognize_header(std::istreambuf_iterator<char> &buf) -> std::optional<std::tuple<physics::ForceModel, int>>;
  static auto recognize_force_model(std::istreambuf_iterator<char> &buf) -> std::optional<physics::ForceModel>;
  static auto recognize_int(std::istreambuf_iterator<char> &buf) -> std::optional<int>;
  static auto recognize_double(std::istreambuf_iterator<char> &buf) -> std::optional<double>;
  static auto recognize_double_triplet(std::istreambuf_iterator<char> &buf) -> std::optional<std::array<double, 3>>;
  static auto recognize_dimension_triplet(std::istreambuf_iterator<char> &buf) -> std::optional<std::array<int, 3>>;
  static auto recognize_cuboid(std::istreambuf_iterator<char> &buf) -> std::optional<cuboid_t>;
  static auto recognize_planet(std::istreambuf_iterator<char> &buf) -> std::optional<Particle>;
  static auto parse_gravity(std::istreambuf_iterator<char> &buf) -> std::optional<std::vector<Particle>>;
  static auto parse_cuboids(std::istreambuf_iterator<char> &buf) -> std::optional<std::vector<cuboid_t>>;
  static auto recognize_end_of_file(std::istreambuf_iterator<char> &buf) -> bool;

 public:
  explicit ParticleLoader(const std::shared_ptr<config::Config> &config);
  auto load_particles() -> std::tuple<ParticleContainer, physics::ForceModel>;
};

}// namespace simulator::io
