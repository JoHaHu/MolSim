//
// Created by TimSc on 03.06.2024.
//

#include "simulator/io/xml_reader/XMLFileReader.h"
#include "simulator/physics/Gravity.h"
#include "utils/ArrayUtils.h"
#include <array>
#include <gtest/gtest.h>
#include <tuple>

class XMLFileReaderTest : public ::testing::Test {
};

/**
 * Test: input_test_gravity WS1
 *
 * Verifies that the given input in the XML file "eingabe-gravity-ws1.xml" is handled properly for the gravitational simulation
 *
 * Ensures that the object that is created for the simulation has the correct values specified in the XML file.
 * Using expects to be able to differentiate between the values and narrow down specific problems e.g. only occurring to some parameters.
 */
TEST_F(XMLFileReaderTest, input_test_advanced) {

  // Header data
  const std::string base_name = "Collision of two bodies";
  const double end_time = 20;
  const int output_frequency = 100;
  const int seed = 1234;
  const std::string output_file_name = "collision_simulation_cuboids";

  // Cuboid 1
  std::vector<double> coordinates_x1 = {20, 20, 2};
  std::vector<double> dimensions_x1 = {100, 20, 1};
  std::vector<double> velocity_x1 = {0, 0, 0};
  Cuboid cuboid1 = Cuboid(coordinates_x1, dimensions_x1, velocity_x1, 1.0, 1.1225, 0);

  // Cuboid 2
  std::vector<double> coordinates_x2 = {70, 60, 2};
  std::vector<double> dimensions_x2 = {20, 20, 1};
  std::vector<double> velocity_x2 = {0, -10, 0};
  Cuboid cuboid2 = Cuboid(coordinates_x2, dimensions_x2, velocity_x2, 1.0, 1.3458, 0);

  std::vector<Cuboid> cuboids{};
  cuboids.emplace_back(cuboid1);
  cuboids.emplace_back(cuboid2);

  std::vector<double> domain_size = {150, 90, 4};
  std::vector<BoundaryCondition> boundary_conditions = {BoundaryCondition::periodic, BoundaryCondition::reflecting,
                                                        BoundaryCondition::periodic, BoundaryCondition::reflecting,
                                                        BoundaryCondition::periodic, BoundaryCondition::periodic};
  const double mass = 1;
  const double delta_t = 0.0005;
  const std::vector<double> sigmas = {1.0};
  const std::vector<double> epsilons = {5.0};

  auto config = XMLFileReader::parseXMLData("../../input/eingabe-collision-test.xml");

  // confirming correctness of header data
  EXPECT_TRUE(config->base_name == base_name);
  EXPECT_TRUE(config->output_file == output_file_name);
  EXPECT_EQ(config->output_frequency, output_frequency);
  EXPECT_TRUE(config->end_time == end_time);
  EXPECT_EQ(config->seed, seed);
  EXPECT_DOUBLE_EQ(config->delta_t, delta_t);

  EXPECT_TRUE(config->particle_container_type == ParticleContainerType::LinkedCells);
  EXPECT_EQ(config->domain_size, domain_size);
  EXPECT_EQ(config->boundary_conditions, boundary_conditions);

  EXPECT_TRUE(config->simulation_type == simulator::physics::ForceModel::LennardJones);
  EXPECT_EQ(config->sigma, sigmas);
  EXPECT_EQ(config->epsilon, epsilons);

  EXPECT_EQ(config->cuboids.at(0).coordinates, cuboids.at(0).coordinates);
  EXPECT_EQ(config->cuboids.at(0).particles, cuboids.at(0).particles);
  EXPECT_EQ(config->cuboids.at(0).velocity, cuboids.at(0).velocity);
  EXPECT_EQ(config->cuboids.at(0).spacing, cuboids.at(0).spacing);

  EXPECT_EQ(config->cuboids.at(1).coordinates, cuboids.at(1).coordinates);
  EXPECT_EQ(config->cuboids.at(1).particles, cuboids.at(1).particles);
  EXPECT_EQ(config->cuboids.at(1).velocity, cuboids.at(1).velocity);
  EXPECT_EQ(config->cuboids.at(1).spacing, cuboids.at(1).spacing);
};
