#pragma once

#include <iostream>
#include <memory>

namespace config {

class Config {
 private:
  Config() = default;

 public:
  auto static parse_config(int argc, char *argv[]) -> std::shared_ptr<Config>;
  double start_time{};
  double end_time{};
  double delta_t{};
  std::string output_filename;
  std::string input_filename;
  std::string log_level = "info";
  int seed{};
  int io_interval{};

  static auto print_usage() -> void;
};

}// namespace config
