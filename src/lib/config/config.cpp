#include "config.h"
#include "spdlog/spdlog.h"
#include <getopt.h>
#include <string>

void config::Config::print_usage() {
  std::cout << "Usage: MolSim [options]\n"
               "--output\t-o\t\t specify output file\n"
               "--input\t\t-i\t\t specify input file\n"
               "--help\t\t-h\t\t open usage page\n"
               "--delta_t\t-d\t\t set delta_t\n"
               "--start_time\t-s\t\t set start_time\n"
               "--seed\t\t-r\t\t set seed for rng\n"
               "--end_time\t-e\t\t set end_time\n"
               "--io_interval\t-l\t\t set plot interval, 0 disables plotting\n";
}

/**
 * @brief Parses command-line arguments to configure the simulation parameters.
 *
 * Reads options from the command line and sets the corresponding configuration values.
 * Prints usage and exits on invalid input.
 *
 * @param argc The number of command-line arguments.
 * @param argv The array of command-line arguments.
 * @return std::shared_ptr<Config> containing the configured parameters.
 *
 * Options:
 * - -o, --output <file>       : Specify the output file.
 * - -i, --input <file>        : Specify the input file.
 * - -h, --help                : Display the usage page.
 * - -d, --delta_t <value>     : Set the delta time (double).
 * - -s, --start_time <value>  : Set the start time (double).
 * - -r, --seed <value>        : Set the seed for RNG (int).
 * - -e, --end_time <value>    : Set the end time (double).
 * - -l, --io_interval <value> : Set the plot interval (int), 0 disables plotting.
 */
auto config::Config::parse_config(int argc, char *argv[]) -> std::shared_ptr<Config> {
  // NOLINT(*-avoid-c-arrays)
  double start_time = 0.0;
  double end_time = 1000.0;// NOLINT(*-avoid-magic-numbers)
  double delta_t = 0.014;  // NOLINT(*-avoid-magic-numbers)
  int seed = 42;           // NOLINT(*-avoid-magic-numbers)
  std::string output_file = "MD_vtk";
  std::string input_file;
  int io_interval = 100;//NOLINT(*-avoid-magic-numbers)
  bool generate_disk = false;
  double disk_center_x = 0;
  double disk_center_y = 0.0;
  double disk_initial_vx = 0.0;
  double disk_initial_vy = 0.0;
  int disk_radius_molecules = 100;
  double disk_meshwidth = 1.0;

  const std::string short_options = "ho:i:d:s:e:r:l:";
  const std::array<option, 16> long_options = {
      {{"output", required_argument, nullptr, 'o'},
       {"input", required_argument, nullptr, 'i'},
       {"help", no_argument, nullptr, '?'},
       {"delta_t", required_argument, nullptr, 'd'},
       {"start_time", required_argument, nullptr, 's'},
       {"seed", required_argument, nullptr, 'r'},
       {"end_time", required_argument, nullptr, 'e'},
       {"io_interval", required_argument, nullptr, 'l'}}};

  while (true) {
    const auto opt = getopt_long(argc, argv, short_options.c_str(), long_options.data(), nullptr);
    if (opt == -1) {
      break;
    }
    switch (opt) {
      case 'o':
        output_file = std::string(optarg);
        break;
      case 'i':
        input_file = std::string(optarg);
        break;
      case 'd':
        delta_t = std::stod(optarg);
        break;
      case 's':
        start_time = std::stod(optarg);
        break;
      case 'e':
        end_time = std::stod(optarg);
        break;
      case 'r':
        seed = std::stoi(optarg);
        break;
      case 'l':
        io_interval = std::stoi(optarg);
        break;
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
  };

  auto config = Config();
  config.output_filename = output_file;
  config.input_filename = input_file;
  config.end_time = end_time;
  config.start_time = start_time;
  config.delta_t = delta_t;
  config.seed = seed;
  config.io_interval = io_interval;

  return std::make_shared<Config>(config);
}
