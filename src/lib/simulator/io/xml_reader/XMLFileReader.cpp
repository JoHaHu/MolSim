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
 * @brief Parses the data from the given XML file using the XSD library.
 *
 * @param xmlFilePath path of the XML file that contains the simulation parameter values
 *
 * creates a SimConfigXML object to store the simulation parameters given by the XML input
 * for further passing and use inside the project
 *
 */

auto XMLFileReader::parseXMLData(const std::string &xmlFilePath) -> std::shared_ptr<SimConfigXML> {
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

    // Log the input values so the user can confirm their correctness
    SPDLOG_INFO("Passing the following parameters to the simulator: ");
    SPDLOG_INFO(data->base_name() + " | output frequency: " + std::to_string(data->output_write_frequency()) + " | t_end value: " + std::to_string(data->t_end()) + " | delta_t value: " + std::to_string(data->delta_t()) + " | epsilon value: " + std::to_string(data->epsilon()) + " | sigma value: " + std::to_string(data->sigma()));

    // vector of cuboids to store all parameters
    auto cuboids = std::vector<Cuboid>{};

    // Iterate through all cuboid elements
    for (const auto &cuboid : data->Cuboid()) {
      std::string logOutput = "Cuboid with specs: coordinate x: [";

      std::vector<int> lower_left_coordinate{};
      std::vector<int> dimensional_particle_numbers{};
      std::vector<int> initial_velocity{};

      auto count = 0;
      for (const auto &arr : cuboid.lower_left_front_coordinate().value()) {
        lower_left_coordinate.push_back(arr);
        logOutput.append(std::to_string(arr));
        if (count++ < 2) {
          logOutput.append(", ");
        }
      }
      logOutput.append("] | particles per dimension: [");

      for (const auto &arr : cuboid.dimensional_particle_numbers().value()) {
        dimensional_particle_numbers.push_back(arr);
        logOutput.append(std::to_string(arr));
        auto const commaConst5 = 5;// because 5 is a magic number
        if (count++ < commaConst5) {
          logOutput.append(", ");
        }
      }
      logOutput.append("] | initial velocity v: [");

      for (const auto &arr : cuboid.initial_velocity().value()) {
        initial_velocity.push_back(arr);
        logOutput.append(std::to_string(arr));
        auto const commaConst8 = 8;// because 8 is a magic number
        if (count++ < commaConst8) {
          logOutput.append(", ");
        }
      }
      logOutput.append("] | mass m of one particle: ");
      logOutput.append(std::to_string(cuboid.distance_h()));
      logOutput.append(" | distance h of particles (mesh width): ");
      logOutput.append(std::to_string(cuboid.mass_m()));

      // initializing arrays for the SimConfigXML object to be created
      std::array<int, 3> arr_lower_left_coordinate{};
      std::array<int, 3> arr_dimensional_particle_numbers{};
      std::array<int, 3> arr_initial_velocity{};

      // checking if lengths match before copying to static length arrays of size 3, otherwise give error
      if (lower_left_coordinate.size() == arr_lower_left_coordinate.size() && dimensional_particle_numbers.size() == arr_dimensional_particle_numbers.size() && initial_velocity.size() == arr_initial_velocity.size()) {
        std::copy(lower_left_coordinate.begin(), lower_left_coordinate.end(), arr_lower_left_coordinate.begin());
        std::copy(dimensional_particle_numbers.begin(), dimensional_particle_numbers.end(), arr_dimensional_particle_numbers.begin());
        std::copy(initial_velocity.begin(), initial_velocity.end(), arr_initial_velocity.begin());
      } else {
        spdlog::warn("There was an error while parsing the values for the cuboids.");
        return nullptr;
      }

      SPDLOG_INFO(logOutput);

      // Create an instance of Cuboid
      // append the cuboid to the vector of cuboids for the config
      auto dist_h = (double) cuboid.distance_h();
      auto mass_m = (double) cuboid.mass_m();
      auto temp_cuboid = Cuboid(arr_lower_left_coordinate, arr_dimensional_particle_numbers, dist_h, mass_m, arr_initial_velocity);
      cuboids.emplace_back(temp_cuboid);
    }

    // create sim_config object to store all the parameters
    auto sim_config = SimConfigXML::store_config_values(data->base_name(), data->output_write_frequency(), data->t_end(), data->delta_t(), data->epsilon(), data->sigma(), data->average_brownian_motion(), cuboids);
    return sim_config;

  } catch (const xml_schema::exception &exception) {
    spdlog::warn("An error has occurred! Please look at the exception details here:\n");
    std::cerr << exception << '\n';
  }
  XMLPlatformUtils::Terminate();
  return nullptr;
}