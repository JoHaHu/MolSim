#pragma once

#include "lib/config/config.h"
#include <spdlog/sinks/basic_file_sink.h>

class LoggerManager {
 public:
  LoggerManager() = delete;
  static void setup_logger(std::shared_ptr<config::Config> config, spdlog::level::level_enum level = spdlog::level::from_str("info"));
};
