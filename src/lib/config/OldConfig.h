#pragma once

#include <iostream>
#include <memory>

namespace OldConfig {

/**
 * All supported runtime configurations
 * */
class OldConfig {
 private:
  OldConfig() = default;

 public:
  /**
   * build a config from provided cmdline arguments.
   * */
  auto static parse_config(int argc, char *argv[]) -> std::shared_ptr<OldConfig>;
  /**
   * not used right now
   * */
  double start_time{};
  /**
   * the end time of the simulation
   * */
  double end_time{};
  /**
   * the timestep size
   * */
  double delta_t{};
  /**
   * the output file name
   * */
  std::string output_filename;
  /**
   * the input filename
   * */
  std::string input_filename;
  /**
   * the seed for the rng
   * */
  int seed{};
  /**
   * the number of simulation step after which output files are generated
   * */
  int io_interval{};

  /**
   * prints the help message
   */
  static auto print_usage() -> void;
};

}// namespace config
