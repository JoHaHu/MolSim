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
    SPDLOG_ERROR("Input file not specified");
    print_usage();
    exit(0);
  }

  // Call XMLFileReader::parseXMLData to set additional parameters
  auto config = XMLFileReader::parseXMLData("../input/" + input_file);
  return config;
}

}// namespace config
