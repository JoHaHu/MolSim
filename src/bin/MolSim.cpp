

#include "lib/FileReader.h"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <vector>

#include "lib/utils/LoggerManager.h"
#include <boost/program_options.hpp>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace po = boost::program_options;

/**
 * Parses command-line options.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @param vm Stores parsed options.
 * @return True if parsing succeeds and it's okay to proceed, False if help is displayed or parsing fails.
 */
bool InitializeOptions(int argc, char *argv[], po::variables_map &vm) {
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "produce help message")("filename,f", po::value<std::string>(), "input filename")("end_time,e", po::value<double>(), "simulation end time")("delta_t,d", po::value<double>(), "time step size")("log_level,l", po::value<std::string>()->default_value("info"),
                                                                                                                                                                                                                               "log level (trace, debug, info, warn, error, critical)");

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  } catch (po::error &e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    std::cerr << desc << std::endl;
    return false;
  }

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return false;
  }

  // Check input validity
  if (!vm.count("filename") || !vm.count("end_time") || !vm.count("delta_t")) {
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
int main(int argc, char *argv[]) {
  po::variables_map vm;
  if (!InitializeOptions(argc, argv, vm)) {
    return 1;
  }

  SetupLogger(vm["log_level"].as<std::string>());
  spdlog::info("Hello from MolSim for PSE!");

  std::string filename = vm["filename"].as<std::string>();
  double end_time = vm["end_time"].as<double>();
  double delta_t = vm["delta_t"].as<double>();
  std::string log_level = vm["log_level"].as<std::string>();

  spdlog::debug("Filename: {}", filename);
  spdlog::debug("End Time: {}", end_time);
  spdlog::debug("Delta T: {}", delta_t);
  spdlog::debug("Log Level: {}", log_level);

  std::vector<Particle> particles;
  FileReader fileReader;
  try {
    spdlog::info("Reading file: " + filename);
    fileReader.readFile(particles, filename);
  } catch (const std::exception &e) {
    spdlog::error("Failed to read file: {}", e.what());
    return 1;
  }

  auto plotter = VTKPlotter();
  auto physics = Gravity();

  auto p = std::vector<Particle>();
  FileReader::readFile(p, filename);
  auto particles = ParticleContainer(p);

  auto simulator = Simulator(particles, plotter, physics, 0, end_time, delta_t);
  simulator.run();

  return 0;
}
