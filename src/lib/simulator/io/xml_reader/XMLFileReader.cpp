//
// Created by TimSc on 01.06.2024.
//

#include "XMLFileReader.h"
#include "SimulationInputSchema.hxx"
#include <iostream>
#include <spdlog/spdlog.h>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/PlatformUtils.hpp>

using namespace xercesc;

void printXMLData(const std::string &xmlFilePath) {
  try {
    XMLPlatformUtils::Initialize();
  } catch (const XMLException &exception) {
    char *message = XMLString::transcode(exception.getMessage());
    spdlog::warn("Error during XML file initialization!\n");
    SPDLOG_INFO(message);
    return;
  }

  try {
    // Parse the XML file
    std::unique_ptr<Data> data = Data_(xmlFilePath);
    SPDLOG_INFO("Reading XML input to start simulation.");

    // Log the input values so the user can confirm their correctness
    SPDLOG_INFO("Passing the following parameters to the simulator: ");
    SPDLOG_INFO(data->base_name() + ": output frequency: " + std::to_string(data->output_write_frequency()) + " | t_end value: " + std::to_string(data->t_end()) + " | delta_t value: " + std::to_string(data->delta_t()) + " | epsilon value: " + std::to_string(data->epsilon()) + " | sigma value: " + std::to_string(data->sigma()));

    // Iterate through all cuboid elements
    for (const auto &cuboid : data->Cuboid()) {
      std::string logOutput = "Cuboid with specs: coordinate x: [";

      auto count = 0;
      for (const auto &arr : cuboid.lower_left_front_coordinate().value()) {
        logOutput.append(std::to_string(arr));
        if (count++ < 2) {
          logOutput.append(", ");
        }
      }
      logOutput.append("] | particles per dimension: [");

      for (const auto &arr : cuboid.dimensional_particle_numbers().value()) {
        logOutput.append(std::to_string(arr));
        auto const commaConst5 = 5;// because 5 is a magic number
        if (count++ < commaConst5) {
          logOutput.append(", ");
        }
      }
      logOutput.append("] | initial velocity v: [");

      for (const auto &arr : cuboid.initial_velocity().value()) {
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

      SPDLOG_INFO(logOutput);
    }
  } catch (const xml_schema::exception &exception) {
    spdlog::warn("An error has occurred! Please look at the exception details here:\n");
    std::cerr << exception << '\n';
  }
  XMLPlatformUtils::Terminate();
}
