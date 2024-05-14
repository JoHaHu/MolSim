#pragma once

#include <spdlog/sinks/basic_file_sink.h>

class LoggerManager {
 public:
  static void setupLogger(spdlog::level::level_enum level);
};
