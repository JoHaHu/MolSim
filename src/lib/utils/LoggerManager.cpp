#include "LoggerManager.h"

#include <iostream>

#include "spdlog/cfg/env.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

/**
 * Initializes the default logger with console and file sinks.
 *
 * Sets up a default logger to log to both the console and a file. Loads log
 * levels from environment variables and sets the specified log level.
 *
 * @param level The log level to set for the logger.
 */
void LoggerManager::setupLogger(spdlog::level::level_enum level) {
  try {
    // Create logger
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("logs/log.txt", true);

    std::vector<spdlog::sink_ptr> sinks{console_sink, file_sink};
    auto logger = std::make_shared<spdlog::logger>("multi_sink", sinks.begin(), sinks.end());

    spdlog::set_default_logger(logger);
    spdlog::cfg::load_env_levels();
    spdlog::set_level(level);
    spdlog::set_pattern("[%Y-%m-%d %H:%M:%S] [%^%l%$] %v");

    spdlog::info("Logger initialized with level: {}", spdlog::level::to_string_view(level));
  } catch (const spdlog::spdlog_ex &ex) {
    std::cerr << "Log initialization failed: " << ex.what() << std::endl;
  }
}
