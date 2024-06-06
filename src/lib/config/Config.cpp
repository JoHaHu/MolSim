#include "Config.h"
#include "simulator/io/xml_reader/XMLFileReader.h"
#include "spdlog/spdlog.h"
#include <getopt.h>
#include <utility>

namespace config {

void Config::print_usage() {
  std::cout << "Usage: MolSim [options]\n"
               "--input\t\t-i\t\t specify input XML file\n"
               "--help\t\t-h\t\t open usage page\n";
}

/**
 * @brief Creates an object to store the simulation parameters given by the input XML and parsed by XSD.
 *
 *
 * @param simulationType The enumerated type of simulation (differentiating between Gravity and Lennard Jones).
 * @param bodyType The enumerated type of body (differentiating between cuboids and discs).
 * @param baseName The base name.
 * @param endTime The simulation end time.
 * @param output_frequency The write frequency of output files.
 * @param outputFilename The name of the output file.
 * @param totalBodies For gravitational simulation: the number of celestial bodies in the simulation.
 * @param deltaT The time steps of t.
 * @param inputSigma The sigma constant used in calculation of Lennard-Jones forces.
 * @param inputEpsilon The epsilon constant used in calculation of Lennard-Jones forces.
 * @param massM The mass m of one particle in the Lennard Jones force simulation.
 * @param distanceH The distance between particles (mesh width of the grid) in the Lennard Jones force simulation.
 * @param averageBrownianMotion The average brownian motion velocity used for calculation.
 * @param celestialBodies A vector of celestial bodies that should be simulated.
 * @param cuboidVector A vector of cuboids that should be simulated.
 * @param discVector A vector of cuboids that should be simulated.
 * @param discVector A vector of discs that should be simulated.
 * @param RngSeed  An integer seed that enables replication of the simulation process.
 *
 * @return std::shared_ptr<Config> containing the configured parameters.
 *
 */
auto Config::store_config_values(simulator::physics::ForceModel simulationType, BodyType bodyType, std::string &baseName, double endTime, double outputFrequency,
                                 std::string &outputFilename, int totalBodies, double deltaT, double inputSigma, double inputEpsilon, double massM, double distanceH,
                                 double averageBrownianMotion, std::vector<CelestialBody> celestialBodies, std::vector<Cuboid> cuboidVector,
                                 std::vector<Disc> discVector, int RngSeed) -> std::shared_ptr<Config> {
  if (simulationType == simulator::physics::ForceModel::LennardJones) {
    Config config = Config();
    config.body_type = bodyType;
    config.simulation_type = simulationType;
    config.base_name = baseName;
    config.end_time = endTime;
    config.output_frequency = outputFrequency;
    config.output_filename = outputFilename;
    config.delta_t = deltaT;
    config.sigma = inputSigma;
    config.epsilon = inputEpsilon;
    config.mass_m = massM;
    config.distance_h = distanceH;
    config.brownian_motion = averageBrownianMotion;
    config.seed = RngSeed;
    if (bodyType == BodyType::Cub) {
      config.cuboids = std::move(cuboidVector);
    } else {
      config.discs = std::move(discVector);
    }
    return std::make_shared<Config>(config);
  }
  if (simulationType == simulator::physics::ForceModel::Gravity) {
    Config config = Config();
    config.body_type = bodyType;
    config.simulation_type = simulationType;
    config.base_name = baseName;
    config.end_time = endTime;
    config.output_frequency = outputFrequency;
    config.output_filename = outputFilename;
    config.total_bodies = totalBodies;
    config.celestial_bodies = std::move(celestialBodies);
    config.seed = RngSeed;
    return std::make_shared<Config>(config);
  }
  return nullptr;
}

/**
 * @brief Parses command-line arguments to configure the simulation parameters.
 *
 * This function reads options from the command line and sets the corresponding
 * configuration values. It then parses an XML file to set additional parameters
 * for the simulation. If the input arguments are invalid or if help is requested,
 * it prints usage information and exits the program.
 *
 * @param argc The number of command-line arguments.
 * @param argv The array of command-line arguments.
 * @return std::shared_ptr<Config> containing the configured parameters.
 *
 * Command-line Options:
 * - -i, --input <file> : Specify the input XML file.
 * - -h, --help         : Display the usage page.
 *
 * The function performs the following steps:
 * 1. Parses command-line arguments to obtain the input filename.
 * 2. Prints usage information and exits on invalid input or if help is requested.
 * 3. Reads the XML file specified by the input filename to set additional parameters.
 * 4. Returns a shared pointer to the Config object containing all configuration values.
 */
auto Config::parse_config(int argc, char *argv[]) -> std::shared_ptr<Config> {
  std::string input_file;

  const std::string short_options = "i:h";
  const std::array<option, 2> long_options = {{{"input", required_argument, nullptr, 'i'},
                                               {"help", no_argument, nullptr, 'h'}}};

  while (true) {
    const auto opt = getopt_long(argc, argv, short_options.c_str(), long_options.data(), nullptr);
    if (opt == -1) {
      break;
    }
    switch (opt) {
      case 'i':
        input_file = std::string(optarg);
        break;
      case 'h':
        print_usage();
        exit(1);
      case '?':
        print_usage();
        exit(1);
      case ':':
        print_usage();
        exit(3);
      default:
        print_usage();
        exit(2);
    }
  }

  if (input_file.empty()) {
    spdlog::error("Input file not specified");
    print_usage();
    exit(1);
  }

  auto config = std::make_shared<Config>();
  config->input_filename = input_file;

  // Call XMLFileReader::parseXMLData to set additional parameters
  auto sim_config = XMLFileReader::parseXMLData("../input/" + config->input_filename);

  // Set additional parameters from sim_config
  config->simulation_type = sim_config->simulation_type;
  config->body_type = sim_config->body_type;
  config->base_name = sim_config->base_name;
  config->end_time = sim_config->end_time;
  config->output_frequency = sim_config->output_frequency;
  config->output_filename = sim_config->output_filename;
  config->total_bodies = sim_config->total_bodies;
  config->delta_t = sim_config->delta_t;
  config->sigma = sim_config->sigma;
  config->epsilon = sim_config->epsilon;
  config->mass_m = sim_config->mass_m;
  config->distance_h = sim_config->distance_h;
  config->brownian_motion = sim_config->brownian_motion;
  config->celestial_bodies = std::move(sim_config->celestial_bodies);
  config->cuboids = std::move(sim_config->cuboids);
  config->discs = std::move(sim_config->discs);
  config->seed = sim_config->seed;

  return config;
}

}// namespace config
