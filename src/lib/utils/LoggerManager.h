#pragma once

#include "config/config.h"
#include <spdlog/sinks/basic_file_sink.h>

class LoggerManager {
 public:
  LoggerManager() = delete;
  static void setup_logger(std::shared_ptr<config::Config> config);
};
