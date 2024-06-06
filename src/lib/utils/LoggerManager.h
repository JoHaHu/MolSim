#pragma once

#include "config/Config.h"
#include <spdlog/sinks/basic_file_sink.h>

/**
 * a helper class providing a static function to setup the logger
 * */
class LoggerManager {
 public:
  LoggerManager() = delete;
  static void setup_logger(const config::Config &config);
};
