#include "lib/simulator/io/FileReader.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>

#include "lib/Particle.h"
#include "lib/utils/LoggerManager.h"
#include "spdlog/spdlog.h"

namespace simulator::io {
void FileReader::read_file(std::vector<Particle> &particles, std::string filename) {
  std::array<double, 3> x{};
  std::array<double, 3> v{};
  double m = NAN;
  int num_particles = 0;

  spdlog::debug("Opening file: {}", filename);

  std::ifstream input_file(filename);
  std::string tmp_string;

  if (input_file.is_open()) {
    getline(input_file, tmp_string);
    spdlog::trace("Read line: {}", tmp_string);

    while (tmp_string.empty() or tmp_string[0] == '#') {
      getline(input_file, tmp_string);
      spdlog::trace("Read line: {}", tmp_string);
    }

    std::istringstream numstream(tmp_string);
    numstream >> num_particles;
    spdlog::debug("Reading {} particles.", num_particles);
    getline(input_file, tmp_string);
    spdlog::trace("Read line: {}", tmp_string);

    for (int i = 0; i < num_particles; i++) {
      std::istringstream datastream(tmp_string);

      for (auto &xj : x) {
        datastream >> xj;
      }
      for (auto &vj : v) {
        datastream >> vj;
      }
      if (datastream.eof()) {
        spdlog::error("Error reading file: EOF reached unexpectedly reading from line {}", i);
        exit(-1);
      }
      datastream >> m;
      particles.emplace_back(x, v, m);
      spdlog::trace("Particle {}: position = ({}, {}, {}), velocity = ({}, {}, {}), mass = {}", i, x[0], x[1], x[2], v[0], v[1], v[2], m);

      getline(input_file, tmp_string);
      spdlog::trace("Read line: {}", tmp_string);
    }
    spdlog::debug("Successfully read {} particles from {}", num_particles, filename);
  } else {
    spdlog::error("Error: could not open file {}", filename);
    exit(-1);
  }
}
}// namespace simulator::io