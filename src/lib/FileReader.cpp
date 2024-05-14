/*
 * FileReader.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "FileReader.h"

#include <cmath>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include "spdlog/spdlog.h"
#include "utils/LoggerManager.h"

FileReader::FileReader() = default;

FileReader::~FileReader() = default;

void FileReader::readFile(std::vector<Particle> &particles, char *filename) {
  std::array<double, 3> x{};
  std::array<double, 3> v{};
  double m = NAN;
  int num_particles = 0;

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

      getline(input_file, tmp_string);
      spdlog::trace("Read line: {}", tmp_string);
    }
  } else {
    spdlog::error("Error: could not open file {}", filename);
    exit(-1);
  }
}
