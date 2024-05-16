#include "LoggerManager.h"

#include "spdlog/cfg/env.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"
#include <iostream>

/**
 * @brief Retrieves the spdlog logging level from the environment variable SPDLOG_LEVEL.
 *
 * Checks the SPDLOG_LEVEL environment variable and converts it to the corresponding
 * spdlog::level::level_enum value. If the environment variable is not set or has
 * an invalid value, it defaults to spdlog::level::info.
 *
 * @return spdlog::level::level_enum The logging level.
 */
spdlog::level::level_enum get_spdlog_level_from_env() {
  const char* env_level = std::getenv("SPDLOG_LEVEL");
  if (env_level) {
    std::string level_str(env_level);
    if (level_str == "trace") return spdlog::level::trace;
    if (level_str == "debug") return spdlog::level::debug;
    if (level_str == "info") return spdlog::level::info;
    if (level_str == "warn") return spdlog::level::warn;
    if (level_str == "error") return spdlog::level::err;
    if (level_str == "critical") return spdlog::level::critical;
    if (level_str == "off") return spdlog::level::off;
  }
  return spdlog::level::info;  // Default level if SPDLOG_LEVEL is not set or invalid
}

/**
 * Initializes the default logger with console and file sinks.
 *
 * Sets up a default logger to log to both the console and a file. Loads log
 * levels from environment variables and sets the specified log level.
 *
 * @param level The log level to set for the logger.
 */
void LoggerManager::setup_logger(std::shared_ptr<config::Config> config) {
  try {
    // Read environment variable
    spdlog::cfg::load_env_levels();
    spdlog::level::level_enum level = get_spdlog_level_from_env();

    // Create console sink
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_level(level);
    spdlog::debug("Console sink created");

    // Create file sink
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("logs/log.txt", true);
    file_sink->set_level(level);
    spdlog::debug("File sink created");

    // Combine sinks into one logger
    std::vector<spdlog::sink_ptr> sinks{console_sink, file_sink};
    auto logger = std::make_shared<spdlog::logger>("multi_sink", sinks.begin(), sinks.end());

    // Set default logger
    spdlog::set_default_logger(logger);
    spdlog::set_level(level);

    // Disable log when io disabled
    if (config->io_interval == 0) {
      spdlog::set_level(spdlog::level::off);
    }

    spdlog::set_pattern("[%Y-%m-%d %H:%M:%S] [%^%l%$] %v");
    spdlog::info("Logger initialized with level: {}", spdlog::level::to_string_view(level));
  } catch (const spdlog::spdlog_ex &ex) {
    std::cerr << "Log initialization failed: " << ex.what() << '\n';
    spdlog::critical("Log initialization failed: {}", ex.what());
  }
}
