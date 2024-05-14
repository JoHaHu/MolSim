#include "lib/FileReader.h"
#include "lib/utils/LoggerManager.h"

#include <getopt.h>
#include <iostream>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <vector>

#include "lib/simulator/Simulator.h"
#include "lib/simulator/io/VTKPlotter.h"
#include "lib/simulator/physics/Gravity.h"

void print_usage() {
  std::cout << "Usage: ./MolSim -f filename -e end_time -d delta_t -l log_level (trace, debug, info, warn, error, critical)" << std::endl;
}

/**
 * Parses command-line options using getopt.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @param filename Stores the input filename.
 * @param end_time Stores the simulation end time.
 * @param delta_t Stores the time step size.
 * @param log_level Stores the log level.
 * @return True if parsing succeeds and it's okay to proceed, False if help is displayed or parsing fails.
 */
auto InitializeOptions(int argc, char *argv[], char *&filename, double &end_time, double &delta_t, std::string &log_level) -> bool {
  int opt = -1;
  while ((opt = getopt(argc, argv, "hf:e:d:l:")) != -1) {
    switch (opt) {
      case 'h':
        print_usage();
        return false;
      case 'f':
        filename = optarg;
        break;
      case 'e':
        end_time = std::stod(optarg);
        break;
      case 'd':
        delta_t = std::stod(optarg);
        break;
      case 'l':
        log_level = std::string(optarg);
        break;
      default:
        print_usage();
        return false;
    }
  }

  // Check input validity
  if (!filename || end_time == 0 || delta_t == 0) {
    spdlog::error(
        "Erroneous program call!\nUsage: ./MolSim -f filename -e end_time -d delta_t -l log_level (trace, debug, info, warn, error, critical)");
    return false;
  }

  return true;
}

/**
 * Sets the logger's level based on input.
 *
 * @param level Log level as a string.
 */
void SetupLogger(const std::string &level) {
  spdlog::level::level_enum log_level = spdlog::level::from_str(level);
  LoggerManager::setupLogger(log_level);
}

/**
 * Main entry point of the MolSim particle simulation program.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @return 0 if successful, 1 on error.
 *
 * Orchestrates initialization, input validation, and simulation. Logs final output or errors.
 */
auto main(int argc, char *argv[]) -> int {
  char *filename = nullptr;
  double end_time = 0.0;
  double delta_t = 0.0;
  std::string log_level = "info";

  if (!InitializeOptions(argc, argv, filename, end_time, delta_t, log_level)) {
    return 1;
  }

  SetupLogger(log_level);

  spdlog::debug("Filename: {}", filename);
  spdlog::debug("End Time: {}", end_time);
  spdlog::debug("Delta T: {}", delta_t);
  spdlog::debug("Log Level: {}", log_level);

  std::vector<Particle> particles;
  FileReader fileReader;
  try {
    spdlog::info("Reading file: {}", filename);
    fileReader.readFile(particles, filename);
  } catch (const std::exception &e) {
    spdlog::error("Failed to read file: {}", e.what());
    return 1;
  }

  auto plotter = VTKPlotter();
  auto physics = Gravity();

  auto p = std::vector<Particle>();
  fileReader.readFile(p, filename);
  auto particleContainer = ParticleContainer(p);

  auto simulator = Simulator(particleContainer, plotter, physics, 0, end_time, delta_t);
  simulator.run();

  return 0;
}
