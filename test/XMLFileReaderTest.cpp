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
 public:
  /**
  // CONSTANTS (to be changed when values in XML or input files are changed)
  **/

  const std::string base_name = "ExampleName";
  const int output_frequency = 100;
  const double end_time = 5;
  const double delta_t = 0.0002;
  const double epsilon = 5;
  const double sigma = 1;

  // values for the first cuboid
  const std::array<int, 3> coordinates_x1 = {0, 0, 0};
  const std::array<int, 3> particles_x1 = {40, 8, 1};
  const std::array<int, 3> velocity_x1 = {0, 0, 0};
  double dist_h_x1 = 1.1225;
  double mass_m_x1 = 1;

  // values for the second cuboid
  const std::array<int, 3> coordinates_x2 = {15, 15, 0};
  const std::array<int, 3> particles_x2 = {8, 8, 1};
  const std::array<int, 3> velocity_x2 = {0, -10, 0};
  double dist_h_x2 = 3.3;
  double mass_m_x2 = 4.4;

  // Helper function to print std::array if necessary
  template<typename T, std::size_t N>
  static void printArray(const std::array<T, N> &arr) {
    std::cout << "[";
    for (std::size_t i = 0; i < N; ++i) {
      std::cout << arr[i];
      if (i < N - 1) {
        std::cout << ", ";
      }
    }
    std::cout << "]" << '\n';
  }
};

/**
 * Test: input_test
 *
 * Verifies that the given input in the XML file "eingabe-falling-drop.xml" is handled properly
 *
 * Ensures that the object that is created for the simulation (SimConfigXML object) has the correct values specified in the XML file.
 * Using expects to be able to differentiate between the values and narrow down specific problems e.g. only occurring to some parameters.
 */
TEST_F(XMLFileReaderTest, input_test) {
  auto config = XMLFileReader::parseXMLData("../../input/eingabe-falling-drop.xml");

  // comparing the string and double values of the config object to the ones defined here as constants and specified in the XML
  EXPECT_TRUE(config->base_name == base_name);
  EXPECT_TRUE(config->output_frequency == output_frequency);
  EXPECT_TRUE(config->end_time == end_time);
  EXPECT_TRUE(config->delta_t == delta_t);
  EXPECT_TRUE(config->sigma == sigma);
  EXPECT_TRUE(config->epsilon == epsilon);

  // checking size of the cuboids vector which should in this case contain 2 cuboids
  // -> if not, further testing can be aborted to not cause indexing errors
  ASSERT_EQ(config->cuboids.size(), 2);

  auto cuboid = config->cuboids.at(0);
  EXPECT_TRUE(coordinates_x1 == get<0>(cuboid));
  EXPECT_TRUE(particles_x1 == get<1>(cuboid));
  EXPECT_TRUE(dist_h_x1 == get<2>(cuboid));
  EXPECT_TRUE(mass_m_x1 == get<3>(cuboid));
  EXPECT_TRUE(velocity_x1 == get<4>(cuboid));

  cuboid = config->cuboids.at(1);
  EXPECT_TRUE(coordinates_x2 == get<0>(cuboid));
  EXPECT_TRUE(particles_x2 == get<1>(cuboid));
  EXPECT_TRUE(dist_h_x2 == get<2>(cuboid));
  EXPECT_TRUE(mass_m_x2 == get<3>(cuboid));
  EXPECT_TRUE(velocity_x2 == get<4>(cuboid));

  // to iterate over a tuple and print out the values ==> for more detailed debugging
  //std::apply([](auto &&...args) { ((std::cout << args << '\n'), ...); }, cuboid);
};
