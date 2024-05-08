#include "FileReader.h"
#include "outputWriter/VTKWriter.h"
#include "outputWriter/XYZWriter.h"
#include "utils/ArrayUtils.h"
#include <iostream>
#include <filesystem>
#include <vector>
#include <boost/program_options.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include "ParticleContainer.h"
#include "utils/LoggerManager.h"

constexpr double start_time = 0;
ParticleContainer particles = ParticleContainer(4);
namespace po = boost::program_options;
namespace fs = std::filesystem;

/**** forward declaration of the calculation functions ****/

/**
 * calculate the force for all particles
 */
void calculateF();

/**
 * calculate the position for all particles
 */
void calculateX(double delta_t);

/**
 * calculate the position for all particles
 */
void calculateV(double delta_t);

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration);

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
    desc.add_options()
            ("help,h", "produce help message")
            ("filename,f", po::value<std::string>(), "input filename")
            ("end_time,e", po::value<double>(), "simulation end time")
            ("delta_t,d", po::value<double>(), "time step size")
            ("log_level,l", po::value<std::string>()->default_value("info"),
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
        auto logger = LoggerManager::getLogger();
        logger->error(
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
    LoggerManager::setLogLevel(log_level);
}

/**
 * Runs the simulation for given particles and parameters.
 *
 * @param particles List of particles.
 * @param start_time Start time of the simulation.
 * @param end_time End time of the simulation.
 * @param delta_t Time step size.
 *
 * Updates particle states and logs progress every X iterations.
 */
void RunSimulation(double start_time, double end_time, double delta_t) {
    auto logger = LoggerManager::getLogger();
    logger->info("Running simulation...");
    int iteration = 0;
    double current_time = start_time;

    while (current_time < end_time) {
        calculateX(delta_t);
        calculateF();
        calculateV(delta_t);

        iteration++;
        if (iteration % 10 == 0) {
            plotParticles(iteration);
            logger->debug("Iteration {} plotted.", iteration);
        }

        logger->debug("Iteration {} finished.", iteration);
        current_time += delta_t;
    }
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
    auto logger = LoggerManager::getLogger();
    logger->info("Hello from MolSim for PSE!");

    std::string filename = vm["filename"].as<std::string>();
    double end_time = vm["end_time"].as<double>();
    double delta_t = vm["delta_t"].as<double>();
    std::string log_level = vm["log_level"].as<std::string>();

    logger->debug("Filename: {}", filename);
    logger->debug("End Time: {}", end_time);
    logger->debug("Delta T: {}", delta_t);
    logger->debug("Log Level: {}", log_level);

    std::vector<Particle> particles;
    FileReader fileReader;
    try {
        logger->info("Reading file: " + filename);
        fileReader.readFile(particles, filename);
    } catch (const std::exception &e) {
        logger->error("Failed to read file: {}", e.what());
        return 1;
    }

    ParticleContainer particleContainer(particles);
    RunSimulation(start_time, end_time, delta_t);

    logger->info("Output written. Terminating...");

    return 0;
}

void calculateF() {
    auto logger = LoggerManager::getLogger();
    logger->debug("Starting force calculation for {} particles.", particles.size());
    for (auto &p: particles) {
        p.old_f = p.f;
        p.f = {0, 0, 0};
    }

    try {
        for (auto pair = particles.begin_pair(); pair != particles.end_pair(); pair++) {
            auto [p1, p2] = *pair;

            auto x_diff = p2.x - p1.x;

            auto norm = ArrayUtils::L2Norm(x_diff);
            if (norm == 0) {
                logger->warn("Zero distance between particles encountered, skipping force calculation for a pair.");
                continue;
            }
            auto f = (p1.m * p2.m) / pow(norm, 3) * x_diff;
            p1.f = p1.f + f;
            p2.f = p2.f - f;
            logger->trace("Force calculation completed.");
        }
    } catch (const std::exception &e) {
        logger->error("Error during force calculation: {}", e.what());
    }
}


void calculateX(double delta_t) {
    auto logger = LoggerManager::getLogger();
    logger->debug("Updating positions for {} particles.", particles.size());

    try {
        for (auto &p: particles) {
            p.x = p.x + delta_t * p.v + (pow(delta_t, 2) / (2 * p.m)) * p.old_f;

            logger->trace("Particle updated.");
        }
        logger->trace("Position updates completed.");
    } catch (const std::exception &e) {
        logger->error("Error during position update: {}", e.what());
    }
}

void calculateV(double delta_t) {
    auto logger = LoggerManager::getLogger();
    logger->debug("Updating velocities for {} particles.", particles.size());

    try {
        for (auto &p: particles) {
            p.v = p.v + delta_t * (1 / (2 * p.m)) * (p.old_f + p.f);

            logger->trace("Particle velocity updated.");
        }
        logger->trace("Velocity updates completed.");
    } catch (const std::exception &e) {
        logger->error("Error during velocity update: {}", e.what());
    }
}

void plotParticles(int iteration) {
    auto logger = LoggerManager::getLogger();
    logger->debug("Starting particle plotting for iteration {}", iteration);

    std::string out_name("MD_vtk");

    try {
        outputWriter::XYZWriter writer;
        writer.plotParticles(particles, out_name, iteration);
        logger->debug("XYZ particle plotting completed for iteration {}", iteration);

        outputWriter::VTKWriter vtk_writer;
        if (particles.size() > INT_MAX) {
            logger->error("Particle size > INT_MAX! Abort.");
            throw std::runtime_error("Particle size exceeds maximum allowable value (INT_MAX), cannot proceed!");
        }
        vtk_writer.initializeOutput(static_cast<int>(particles.size()));

        for (auto &p: particles) {
            vtk_writer.plotParticle(p);
        }
        vtk_writer.writeFile(out_name, iteration);
        logger->debug("VTK particle plotting and file writing completed for iteration {}", iteration);
    } catch (const std::exception &e) {
        logger->error("Error during particle plotting for iteration {}: {}", iteration, e.what());
    }
}
