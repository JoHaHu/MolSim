#include "LoggerManager.h"
#include <vector>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/async.h>

/**
 * Sets the logging level and ensures immediate flushing at this level.
 *
 * @param level The log level to be set for the logger.
 */
void LoggerManager::setLogLevel(spdlog::level::level_enum level) {
    auto logger = LoggerManager::getLogger();
    logger->set_level(level);
    logger->flush_on(level);
}

/**
 * Retrieves the singleton logger instance.
 *
 * @return Reference to the shared logger instance.
 */
std::shared_ptr<spdlog::logger> &LoggerManager::getLogger() {
    static std::shared_ptr<spdlog::logger> logger = LoggerManager::createLogger();
    return logger;
}

/**
 * Creates and configures the logger with console and file sinks.
 *
 * @return A shared pointer to the newly created logger.
 */
std::shared_ptr<spdlog::logger> LoggerManager::createLogger() {
    spdlog::init_thread_pool(8192, 1);
    std::vector<spdlog::sink_ptr> sinks;
    sinks.push_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
    sinks.push_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>("log.txt", true));

    auto logger = std::make_shared<spdlog::async_logger>("console", sinks.begin(), sinks.end(),
                                                         spdlog::thread_pool(),
                                                         spdlog::async_overflow_policy::block);
    logger->set_level(spdlog::level::info); // Default log level

    return logger;
}
