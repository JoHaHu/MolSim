//
// Created by TimSc on 01.06.2024.
//

#include "XMLFileReader.h"
#include "SimulationInputSchema.hxx"
#include <iostream>
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
auto XMLFileReader::parseXMLData(const std::string &xmlFilePath) -> std::shared_ptr<Config> {
  try {
    XMLPlatformUtils::Initialize();
  } catch (const XMLException &exception) {
    char *message = XMLString::transcode(exception.getMessage());
    spdlog::warn("Error during XML file initialization!\n");
    SPDLOG_INFO(message);
    return nullptr;
  }

  try {
    // Parse the XML file and disable validation (according to worksheet)
    std::unique_ptr<Data> data = Data_(xmlFilePath, xml_schema::flags::dont_validate);
    SPDLOG_INFO("Reading XML input to start simulation.");

    // Log the header input values so the user can confirm their correctness
    SPDLOG_INFO("Passing the following header parameters to simulation: " + data->header().base_name() + " | output frequency: " + std::to_string(data->header().output_frequency()) + " | t_end value: " + std::to_string(data->header().t_end()));

    if (data->lennard_jones().present() && data->lennard_jones()->cuboids().present()) {

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
          spdlog::warn("There was an error while parsing the values for the cuboids.");
          return nullptr;
        }

        // Create an instance of Cuboid
        // append the cuboid to the vector of cuboids for the config

        auto temp_cuboid = Cuboid(arr_coordinate, arr_particle_numbers, arr_velocity);
        temp_cuboids.emplace_back(temp_cuboid);
      }

      std::vector<Disc> emptyDiscs;
      std::vector<CelestialBody> emptyCelestials;
      // create sim_config object to store all the parameters
      auto sim_config = Config::store_config_values(SimulationType::LennardJones, BodyType::Cub, data->header().base_name(),
                                                    data->header().t_end(), data->header().output_frequency(), data->header().output_file_name(),
                                                    0, data->lennard_jones()->settings().delta_t(), data->lennard_jones()->settings().sigma(),
                                                    data->lennard_jones()->settings().epsilon(), data->lennard_jones()->settings().mass_m(),
                                                    data->lennard_jones()->settings().distance_h(), data->lennard_jones()->settings().brown_motion(),
                                                    emptyCelestials, temp_cuboids, emptyDiscs, data->header().seed());
      return sim_config;
    }

    else if (data->lennard_jones().present() && data->lennard_jones()->discs().present()) {

      std::vector<Disc> temp_discs;
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
          spdlog::warn("There was an error while parsing the values for the cuboids.");
          return nullptr;
        }

        auto temp_disc = Disc(arr_coordinate, arr_velocity, radius);
        temp_discs.emplace_back(temp_disc);
      }

      std::vector<Cuboid> emptyCuboids;
      std::vector<CelestialBody> emptyCelestials;
      // create sim_config object to store all the parameters
      auto sim_config = Config::store_config_values(SimulationType::LennardJones, BodyType::Dis, data->header().base_name(),
                                                    data->header().t_end(), data->header().output_frequency(), data->header().output_file_name(),
                                                    0, data->lennard_jones()->settings().delta_t(), data->lennard_jones()->settings().sigma(),
                                                    data->lennard_jones()->settings().epsilon(), data->lennard_jones()->settings().mass_m(),
                                                    data->lennard_jones()->settings().distance_h(), data->lennard_jones()->settings().brown_motion(),
                                                    emptyCelestials, emptyCuboids, temp_discs, data->header().seed());
      return sim_config;
    }

    else if (data->gravity().present()) {
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
          spdlog::warn("There was an error while parsing the values for the cuboids.");
          return nullptr;
        }

        auto temp_body = CelestialBody(arr_coordinate, arr_velocity, mass);
        temp_bodies.emplace_back(temp_body);
        body_count++;
      }

      std::vector<Cuboid> emptyCuboids;
      std::vector<Disc> emptyDiscs;
      const double emptyDouble = 0;

      // create sim_config object to store all the parameters
      auto sim_config = Config::store_config_values(SimulationType::Gravity, BodyType::Cub, data->header().base_name(),
                                                    data->header().t_end(), data->header().output_frequency(), data->header().output_file_name(),
                                                    body_count, emptyDouble, emptyDouble, emptyDouble, emptyDouble, emptyDouble,
                                                    emptyDouble, temp_bodies, emptyCuboids, emptyDiscs, data->header().seed());
      return sim_config;
    } else {
      return nullptr;
    }

  } catch (const xml_schema::exception &exception) {
    spdlog::warn("An error has occurred! Please look at the exception details here:\n");
    std::cerr << exception << '\n';
  }
  XMLPlatformUtils::Terminate();
  return nullptr;
}