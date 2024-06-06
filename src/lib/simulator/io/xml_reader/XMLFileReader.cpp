#include "XMLFileReader.h"
#include "SimulationInputSchema.hxx"
#include "container/boundary.h"
#include <iostream>
#include <memory>
#include <ranges>
#include <spdlog/spdlog.h>
#include <vector>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/PlatformUtils.hpp>

using namespace xercesc;

/**
 * @brief Parses the data from the given XML file using the XSD library and stores it in a Config object.
 *
 * @param xmlFilePath path of the XML file that contains the simulation parameter values
 *
 * creates a Config object to store the simulation parameters given by the XML input
 * for further passing and use inside the project as well as easier testing
 *
 */
auto XMLFileReader::parseXMLData(const std::string &xmlFilePath) -> std::shared_ptr<config::Config> {
  try {
    XMLPlatformUtils::Initialize();
  } catch (const XMLException &exception) {
    char *message = XMLString::transcode(exception.getMessage());
    SPDLOG_WARN("Error during XML file initialization!\n");
    SPDLOG_INFO(message);
    return nullptr;
  }

  try {
    // Parse the XML file and disable validation (according to worksheet)
    std::unique_ptr<Data> data = Data_(xmlFilePath, 0);
    SPDLOG_INFO("Reading XML input to start simulation.");

    // Log the header input values so the user can confirm their correctness
    SPDLOG_INFO("Passing the following header parameters to simulation: " + data->header().base_name() + " | output frequency: " + std::to_string(data->header().output_frequency()) + " | t_end value: " + std::to_string(data->header().t_end()));

    config::Config config = config::Config();

    // Header parsing
    config.base_name = data->header().base_name();
    config.end_time = data->header().t_end();
    config.output_frequency = data->header().output_frequency();
    config.output_filename = data->header().output_file_name();
    config.seed = data->header().seed();

    /** Particle Loader choice **/

    // Linked Cells
    if (data->linked_cells().present()) {
      config.particle_loader_type = ParticleContainerType::LinkedCells;

      // Parsing array of domain size (volume)
      const auto &domain_size = data->linked_cells()->domain_size();
      std::vector<double> domain_size_vector;
      for (const auto &arr : domain_size.value()) {
        domain_size_vector.push_back(arr);
      }
      std::array<double, 3> arr_domain_size{};
      if (domain_size_vector.size() == arr_domain_size.size()) {
        std::copy(domain_size_vector.begin(), domain_size_vector.end(), arr_domain_size.begin());
      }
      config.domain_size = arr_domain_size;

      // Parsing array of domain size (volume)
      const auto &boundary_conditions = data->linked_cells()->boundary_conditions();
      std::vector<BoundaryCondition> boundary_conditions_vec;
      for (const auto &arr : boundary_conditions.boundary_condition()) {
        BoundaryCondition bc = BoundaryCondition::none;
        if (arr.type() == "outflow") {
          bc = BoundaryCondition::outflow;
        } else if (arr.type() == "reflecting") {
          bc = BoundaryCondition::reflecting;
        }
        boundary_conditions_vec.push_back(bc);
      }
      std::array<BoundaryCondition, 6> arr_boundary_conditions{};
      if (boundary_conditions_vec.size() == arr_boundary_conditions.size()) {
        std::copy(boundary_conditions_vec.begin(), boundary_conditions_vec.end(), arr_boundary_conditions.begin());
      }
      config.boundary_conditions = arr_boundary_conditions;

      // Vector
    } else if (data->vector().present()) {
      config.particle_loader_type = ParticleContainerType::Vector;
    }

    // Gravity Model
    if (data->gravity().present()) {

      config.simulation_type = simulator::physics::ForceModel::Gravity;

      std::vector<CelestialBody> temp_bodies;
      // Iterate through all celestial bodies
      const auto &bodies = data->gravity()->celestial_body();
      int body_count = 0;
      for (const auto &body : bodies) {

        std::vector<double> coordinate{};
        std::vector<double> velocity{};
        double mass = body.mass();

        for (const auto &arr : body.coordinate().value()) {
          coordinate.push_back(arr);
        }

        for (const auto &arr : body.velocity().value()) {
          velocity.push_back(arr);
        }

        // initializing arrays for the config object to be created
        std::array<double, 3> arr_coordinate{};
        std::array<double, 3> arr_velocity{};

        // checking if lengths match before copying to static length arrays of size 3, otherwise give error
        if (coordinate.size() == arr_coordinate.size() && velocity.size() == arr_velocity.size()) {
          std::copy(coordinate.begin(), coordinate.end(), arr_coordinate.begin());
          std::copy(velocity.begin(), velocity.end(), arr_velocity.begin());
        } else {
          SPDLOG_WARN("There was an error while parsing the values for the cuboids.");
          return nullptr;
        }

        auto temp_body = CelestialBody(arr_coordinate, arr_velocity, mass);
        temp_bodies.emplace_back(temp_body);
        body_count++;
      }

      config.total_bodies = body_count;
      config.celestial_bodies = temp_bodies;
    }

    // Lennard-Jones Model
    if (data->lennard_jones().present()) {

      // parsing settings of Lennard Jones force simulation model
      config.simulation_type = simulator::physics::ForceModel::LennardJones;
      config.delta_t = data->lennard_jones()->settings().delta_t();
      config.sigma = data->lennard_jones()->settings().sigma();
      config.epsilon = data->lennard_jones()->settings().epsilon();
      config.mass_m = data->lennard_jones()->settings().mass_m();
      config.distance_h = data->lennard_jones()->settings().distance_h();
      config.brownian_motion = data->lennard_jones()->settings().brown_motion();
      config.cutoff_radius = data->lennard_jones()->settings().cutoff_radius();

      // Cuboids
      if (data->lennard_jones()->cuboids().present()) {

        config.body_type = BodyType::Cub;
        std::vector<Cuboid> temp_cuboids;
        // Iterate through all cuboid elements
        const auto &cuboids = data->lennard_jones()->cuboids().get().Cuboid();
        for (const auto &cuboid : cuboids) {

          std::vector<double> coordinate{};
          std::vector<double> particle_numbers{};
          std::vector<double> velocity{};

          for (const auto &arr : cuboid.coordinate().value()) {
            coordinate.push_back(arr);
          }

          for (const auto &arr : cuboid.particle_counts().value()) {
            particle_numbers.push_back(arr);
          }

          for (const auto &arr : cuboid.velocity().value()) {
            velocity.push_back(arr);
          }

          // initializing arrays for the config object to be created
          std::array<double, 3> arr_coordinate{};
          std::array<double, 3> arr_particle_numbers{};
          std::array<double, 3> arr_velocity{};

          // checking if lengths match before copying to static length arrays of size 3, otherwise give error
          if (coordinate.size() == arr_coordinate.size() && particle_numbers.size() == arr_particle_numbers.size() && velocity.size() == arr_velocity.size()) {
            std::copy(coordinate.begin(), coordinate.end(), arr_coordinate.begin());
            std::copy(particle_numbers.begin(), particle_numbers.end(), arr_particle_numbers.begin());
            std::copy(velocity.begin(), velocity.end(), arr_velocity.begin());
          } else {
            SPDLOG_WARN("There was an error while parsing the values for the cuboids.");
            return nullptr;
          }

          // Create an instance of Cuboid
          // append the cuboid to the vector of cuboids for the config

          auto temp_cuboid = Cuboid(arr_coordinate, arr_particle_numbers, arr_velocity);
          temp_cuboids.emplace_back(temp_cuboid);
        }

        config.cuboids = temp_cuboids;

      } else if (data->lennard_jones()->discs().present()) {

        // Discs
        std::vector<Disc> temp_discs;

        config.body_type = BodyType::Dis;
        // Iterate through all disc elements
        const auto &discs = data->lennard_jones()->discs().get().Disc();
        for (const auto &disc : discs) {

          std::vector<double> coordinate{};
          std::vector<double> velocity{};
          int radius = disc.radius();

          for (const auto &arr : disc.coordinate().value()) {
            coordinate.push_back(arr);
          }

          for (const auto &arr : disc.velocity().value()) {
            velocity.push_back(arr);
          }

          // initializing arrays for the config object to be created
          std::array<double, 3> arr_coordinate{};
          std::array<double, 3> arr_velocity{};

          // checking if lengths match before copying to static length arrays of size 3, otherwise give error
          if (coordinate.size() == arr_coordinate.size() && velocity.size() == arr_velocity.size()) {
            std::copy(coordinate.begin(), coordinate.end(), arr_coordinate.begin());
            std::copy(velocity.begin(), velocity.end(), arr_velocity.begin());
          } else {
            SPDLOG_WARN("There was an error while parsing the values for the cuboids.");
            return nullptr;
          }

          auto temp_disc = Disc(arr_coordinate, arr_velocity, radius);
          temp_discs.emplace_back(temp_disc);
        }

        config.discs = temp_discs;
      }

      // Spheres
      else if (data->lennard_jones()->spheres().present()) {

        config.body_type = BodyType::Sph;
        std::vector<Sphere> temp_spheres;
        // Iterate through all disc elements
        const auto &spheres = data->lennard_jones()->spheres().get().Sphere();
        for (const auto &sphere : spheres) {

          std::vector<double> coordinate{};
          std::vector<double> velocity{};
          int radius = sphere.radius();

          for (const auto &arr : sphere.coordinate().value()) {
            coordinate.push_back(arr);
          }

          for (const auto &arr : sphere.velocity().value()) {
            velocity.push_back(arr);
          }

          // initializing arrays for the config object to be created
          std::array<double, 3> arr_coordinate{};
          std::array<double, 3> arr_velocity{};

          // checking if lengths match before copying to static length arrays of size 3, otherwise give error
          if (coordinate.size() == arr_coordinate.size() && velocity.size() == arr_velocity.size()) {
            std::copy(coordinate.begin(), coordinate.end(), arr_coordinate.begin());
            std::copy(velocity.begin(), velocity.end(), arr_velocity.begin());
          } else {
            SPDLOG_WARN("There was an error while parsing the values for the cuboids.");
            return nullptr;
          }

          auto temp_sphere = Sphere(arr_coordinate, arr_velocity, radius);
          temp_spheres.emplace_back(temp_sphere);
        }

        config.spheres = temp_spheres;
      }

      // Tori
      else if (data->lennard_jones()->tori().present()) {

        config.body_type = BodyType::Tor;
        std::vector<Torus> temp_tori;
        // Iterate through all disc elements
        const auto &tori = data->lennard_jones()->tori().get().Torus();
        for (const auto &torus : tori) {

          std::vector<double> coordinate{};
          std::vector<double> velocity{};
          double major_radius = torus.major_radius();
          double minor_radius = torus.minor_radius();

          for (const auto &arr : torus.coordinate().value()) {
            coordinate.push_back(arr);
          }

          for (const auto &arr : torus.velocity().value()) {
            velocity.push_back(arr);
          }

          // initializing arrays for the config object to be created
          std::array<double, 3> arr_coordinate{};
          std::array<double, 3> arr_velocity{};

          // checking if lengths match before copying to static length arrays of size 3, otherwise give error
          if (coordinate.size() == arr_coordinate.size() && velocity.size() == arr_velocity.size()) {
            std::copy(coordinate.begin(), coordinate.end(), arr_coordinate.begin());
            std::copy(velocity.begin(), velocity.end(), arr_velocity.begin());
          } else {
            SPDLOG_WARN("There was an error while parsing the values for the cuboids.");
            return nullptr;
          }

          auto temp_torus = Torus(arr_coordinate, arr_velocity, major_radius, minor_radius);
          temp_tori.emplace_back(temp_torus);
        }
        config.tori = temp_tori;
      }

      // Double Helices
      else if (data->lennard_jones()->double_helices().present()) {

        config.body_type = BodyType::Hel;
        std::vector<DoubleHelix> temp_helices;
        // Iterate through all disc elements
        const auto &helices = data->lennard_jones()->double_helices().get().double_helix();
        for (const auto &helix : helices) {

          std::vector<double> coordinate{};
          std::vector<double> velocity{};
          double radius = helix.radius();
          double pitch = helix.pitch();
          double height = helix.height();

          for (const auto &arr : helix.coordinate().value()) {
            coordinate.push_back(arr);
          }

          for (const auto &arr : helix.velocity().value()) {
            velocity.push_back(arr);
          }

          // initializing arrays for the config object to be created
          std::array<double, 3> arr_coordinate{};
          std::array<double, 3> arr_velocity{};

          // checking if lengths match before copying to static length arrays of size 3, otherwise give error
          if (coordinate.size() == arr_coordinate.size() && velocity.size() == arr_velocity.size()) {
            std::copy(coordinate.begin(), coordinate.end(), arr_coordinate.begin());
            std::copy(velocity.begin(), velocity.end(), arr_velocity.begin());
          } else {
            SPDLOG_WARN("There was an error while parsing the values for the cuboids.");
            return nullptr;
          }
        }
        config.double_helices = temp_helices;
      }
    }
    return std::make_shared<config::Config>(config);

  } catch (const xml_schema::exception &exception) {
    SPDLOG_WARN("An error has occurred! Please look at the exception details here:\n");
    std::cerr << exception << '\n';
  }
  XMLPlatformUtils::Terminate();
  return nullptr;
}