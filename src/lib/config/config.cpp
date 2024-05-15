
#include "config.h"
#include "spdlog/spdlog.h"
#include <getopt.h>

void print_usage() {
  std::cout << "Usage: MolSim [options]\n"
               "--output\t-o\t\t specify output file\n"
               "--input\t\t-i\t\t specify input file\n"
               "--help\t\t-h\t\t open usage page\n"
               "--delta_t\t-d\t\t set delta_t\n"
               "--start_time\t-s\t\t set start_time\n"
               "--seed\t\t-r\t\t set seed for rng\n"
               "--task\t\t-t\t\t select task\n"
               "--end_time\t-e\t\t set end_time\n";
}

auto config::Config::parse_config(int argc, char *argv[]) -> std::shared_ptr<Config> {// NOLINT(*-avoid-c-arrays)
  double start_time = 0.0;
  double end_time = 1000.0;// NOLINT(*-avoid-magic-numbers)
  double delta_t = 0.014;  // NOLINT(*-avoid-magic-numbers)
  int seed = 42;           // NOLINT(*-avoid-magic-numbers)
  std::string output_file;
  std::string input_file;

  const std::string short_options = "ho:i:d:s:e:r:";
  const std::array<option, 8> long_options = {{
      {"output", required_argument, nullptr, 'o'},
      {"input", required_argument, nullptr, 'i'},
      {"help", no_argument, nullptr, '?'},
      {"delta_t", required_argument, nullptr, 'd'},
      {"start_time", required_argument, nullptr, 's'},
      {"seed", required_argument, nullptr, 'r'},
      {"task", required_argument, nullptr, 't'},
      {"end_time", required_argument, nullptr, 'e'},
  }};

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

  if (config.input_filename.empty()) {
    print_usage();
    exit(1);
  }

  return std::make_shared<Config>(config);
}
