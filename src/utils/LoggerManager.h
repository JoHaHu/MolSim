#ifndef LOGGER_MANAGER_H
#define LOGGER_MANAGER_H

#include <memory>
#include <spdlog/spdlog.h>

class LoggerManager {
public:
    static void setLogLevel(spdlog::level::level_enum level);

    static std::shared_ptr<spdlog::logger>& getLogger();

private:
    static std::shared_ptr<spdlog::logger> createLogger();
};

#endif // LOGGER_MANAGER_H
