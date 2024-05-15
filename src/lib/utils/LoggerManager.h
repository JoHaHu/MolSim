#pragma once

#include <spdlog/sinks/basic_file_sink.h>

class LoggerManager {
 public:
  LoggerManager() = delete;
  static void setup_logger(spdlog::level::level_enum level = spdlog::level::from_str("info"));
};
