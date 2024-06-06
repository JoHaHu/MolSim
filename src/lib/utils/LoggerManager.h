#pragma once

#include "config/Config.h"
#include <spdlog/sinks/basic_file_sink.h>

#include "config/Config.h"

/**
 * a helper class providing a static function to setup the logger
 * */
class LoggerManager {
 public:
  LoggerManager() = delete;
  // TODO remove the class and only keep a static function
  static void setup_logger(const std::shared_ptr<config::Config> &config);
};
