#include "LoggerManager.h"

#include "spdlog/cfg/env.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"
#include <iostream>

#include "config/Config.h"

/**
 * Initializes the default logger with console and file sinks.
 *
 * Sets up a default logger to log to both the console and a file. Loads log
 * levels from environment variables and sets the specified log level.
 *
 * @param level The log level to set for the logger.
 */
void LoggerManager::setup_logger(const config::Config &config) {
  try {
    // Read environment variable
    spdlog::cfg::load_env_levels();
    spdlog::level::level_enum level = spdlog::get_level();

    // Create console sink
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_level(level);
    spdlog::trace("Console sink created");

    // Create file sink
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("logs/log.txt", true);
    file_sink->set_level(level);
    spdlog::trace("File sink created");

    // Combine sinks into one logger
    std::vector<spdlog::sink_ptr> sinks{console_sink, file_sink};
    auto logger = std::make_shared<spdlog::logger>("multi_sink", sinks.begin(), sinks.end());

    // Set default logger
    spdlog::set_default_logger(logger);
    spdlog::set_level(level);

    // Disable log when io disabled
    if (config.output_frequency == 0) {
      spdlog::set_level(spdlog::level::off);
    }
    spdlog::set_pattern("[%^%l%$] %v");
    spdlog::debug("Logger initialized with level: {}", spdlog::level::to_string_view(level));
  } catch (const spdlog::spdlog_ex &ex) {
    std::cerr << "Log initialization failed: " << ex.what() << '\n';
    spdlog::critical("Log initialization failed: {}", ex.what());
  }
}
