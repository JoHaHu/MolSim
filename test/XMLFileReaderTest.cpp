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
TEST_F(XMLFileReaderTest, input_test_gravity) {

  // Header data
  const std::string base_name = "Halley's Comet";
  const double end_time = 5;
  const double output_frequency = 100;
  const int seed = 1234;
  const std::string output_file_name = "sonne_simulation";

  const std::array<double, 3> coordinates_x1 = {0, 0, 0};
  const std::array<double, 3> velocity_x1 = {0, 0, 0};
  const double mass_x1 = 1;

  const std::array<double, 3> coordinates_x2 = {0, 1, 0};
  const std::array<double, 3> velocity_x2 = {-1, 0, 0};
  const double mass_x2 = 3.0e-6;

  const std::array<double, 3> coordinates_x3 = {0, 5.36, 0};
  const std::array<double, 3> velocity_x3 = {-0.425, 0, 0};
  const double mass_x3 = 9.55e-4;

  const std::array<double, 3> coordinates_x4 = {34.75, 0, 0};
  const std::array<double, 3> velocity_x4 = {0, 0.0296, 0};
  const double mass_x4 = 1.0e-14;

  auto config = XMLFileReader::parseXMLData("../../input/eingabe-gravity-ws1.xml");

  // confirming correctness of header data
  EXPECT_TRUE(config->base_name == base_name);
  EXPECT_TRUE(config->end_time == end_time);
  EXPECT_EQ(config->output_frequency, output_frequency);
  EXPECT_EQ(config->seed, seed);
  EXPECT_TRUE(config->output_filename == output_file_name);

  // checking size of the cuboids vector which should in this case contain 2 cuboids
  // -> if not, further testing can be aborted to not cause indexing errors
  ASSERT_EQ(config->total_bodies, 4);

  auto celestial_body = config->celestial_bodies.at(0);
  EXPECT_TRUE(coordinates_x1 == celestial_body.coordinates);
  EXPECT_TRUE(velocity_x1 == celestial_body.velocity);
  EXPECT_DOUBLE_EQ(mass_x1, celestial_body.mass);

  celestial_body = config->celestial_bodies.at(1);
  EXPECT_TRUE(coordinates_x2 == celestial_body.coordinates);
  EXPECT_TRUE(velocity_x2 == celestial_body.velocity);
  EXPECT_DOUBLE_EQ(mass_x2, celestial_body.mass);

  celestial_body = config->celestial_bodies.at(2);
  EXPECT_TRUE(coordinates_x3 == celestial_body.coordinates);
  EXPECT_TRUE(velocity_x3 == celestial_body.velocity);
  EXPECT_DOUBLE_EQ(mass_x3, celestial_body.mass);

  celestial_body = config->celestial_bodies.at(3);
  EXPECT_TRUE(coordinates_x4 == celestial_body.coordinates);
  EXPECT_TRUE(velocity_x4 == celestial_body.velocity);
  EXPECT_DOUBLE_EQ(mass_x4, celestial_body.mass);
};

/**
 * Test: input_test_cuboids WS2
 *
 * Verifies that the given input in the XML file "eingabe-collision-ws2.xml" is handled properly
 *
 * Ensures that the object that is created for the simulation (Config object) has the correct values specified in the XML file.
 * Using expects to be able to differentiate between the values and narrow down specific problems e.g. only occurring to some parameters.
 */
TEST_F(XMLFileReaderTest, input_test_cuboids) {

  // Header data
  const std::string base_name = "Collision of two bodies";
  const double end_time = 5;
  const double output_frequency = 100;
  const int seed = 1234;
  const std::string output_file_name = "collision_simulation";

  // Lennard Jones data
  const double dist_h = 1.1225;
  const double mass_m = 1;
  const double sigma = 1;
  const double epsilon = 5;
  const double delta_t = 0.0002;
  const double brownian_motion = 0.1;

  // values for the first cuboid
  const std::array<double, 3> coordinates_x1 = {0, 0, 0};
  const std::array<double, 3> particles_x1 = {40, 8, 1};
  const std::array<double, 3> velocity_x1 = {0, 0, 0};

  // values for the second cuboid
  const std::array<double, 3> coordinates_x2 = {15, 15, 0};
  const std::array<double, 3> particles_x2 = {8, 8, 1};
  const std::array<double, 3> velocity_x2 = {0, -10, 0};

  auto config = XMLFileReader::parseXMLData("../../input/eingabe-collision-ws2.xml");

  // confirming correctness of header data
  EXPECT_TRUE(config->base_name == base_name);
  EXPECT_TRUE(config->end_time == end_time);
  EXPECT_EQ(config->output_frequency, output_frequency);
  EXPECT_EQ(config->seed, seed);
  EXPECT_TRUE(config->output_filename == output_file_name);

  // confirming correctness of Lennard Jones settings
  EXPECT_EQ(config->sigma, sigma);
  EXPECT_EQ(config->epsilon, epsilon);
  EXPECT_EQ(config->mass_m, mass_m);
  EXPECT_EQ(config->distance_h, dist_h);
  EXPECT_EQ(config->delta_t, delta_t);
  EXPECT_EQ(config->brownian_motion, brownian_motion);

  // checking size of the cuboids vector which should in this case contain 2 cuboids
  // -> if not, further testing can be aborted to not cause indexing errors
  ASSERT_EQ(config->cuboids.size(), 2);

  auto cuboid = config->cuboids.at(0);
  EXPECT_TRUE(coordinates_x1 == cuboid.coordinates);
  EXPECT_TRUE(velocity_x1 == cuboid.velocity);
  EXPECT_TRUE(particles_x1 == cuboid.particles);

  cuboid = config->cuboids.at(1);
  EXPECT_TRUE(coordinates_x2 == cuboid.coordinates);
  EXPECT_TRUE(velocity_x2 == cuboid.velocity);
  EXPECT_TRUE(particles_x2 == cuboid.particles);
};

/**
 * Test: input_test_discs WS3
 *
 * Verifies that the given input in the XML file "eingabe-disc-ws3.xml" is handled properly
 *
 * Ensures that the object that is created for the simulation (Config object) has the correct values specified in the XML file.
 * Using expects to be able to differentiate between the values and narrow down specific problems e.g. only occurring to some parameters.
 */
TEST_F(XMLFileReaderTest, input_test_disc) {

  // Header data
  const std::string base_name = "Simulation of a falling drop - Wall";
  const double end_time = 10;
  const double output_frequency = 100;
  const int seed = 1234;
  const std::string output_file_name = "falling_drop_simulation";

  // Lennard Jones data
  const double dist_h = 1.1225;
  const double mass_m = 1;
  const double sigma = 1;
  const double epsilon = 5;
  const double delta_t = 0.00005;
  const double brownian_motion = 0.1;

  // values for the disc
  const std::array<double, 3> coordinates_x1 = {60, 25, 0};
  const std::array<double, 3> velocity_x1 = {0, -10, 0};
  const int radius = 15;

  auto config = XMLFileReader::parseXMLData("../../input/eingabe-disc-ws3.xml");

  // confirming correctness of header data
  EXPECT_TRUE(config->base_name == base_name);
  EXPECT_TRUE(config->end_time == end_time);
  EXPECT_EQ(config->output_frequency, output_frequency);
  EXPECT_EQ(config->seed, seed);
  EXPECT_TRUE(config->output_filename == output_file_name);

  // confirming correctness of Lennard Jones settings
  EXPECT_EQ(config->sigma, sigma);
  EXPECT_EQ(config->epsilon, epsilon);
  EXPECT_EQ(config->mass_m, mass_m);
  EXPECT_EQ(config->distance_h, dist_h);
  EXPECT_EQ(config->delta_t, delta_t);
  EXPECT_EQ(config->brownian_motion, brownian_motion);

  // checking size of the cuboids vector which should in this case contain 2 cuboids
  // -> if not, further testing can be aborted to not cause indexing errors
  ASSERT_EQ(config->discs.size(), 1);

  auto disc = config->discs.at(0);
  EXPECT_TRUE(coordinates_x1 == disc.coordinates);
  EXPECT_TRUE(velocity_x1 == disc.velocity);
  EXPECT_EQ(radius, disc.radius);
}
