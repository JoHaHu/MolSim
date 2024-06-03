//
// Created by TimSc on 03.06.2024.
//

#include "simulator/io/xml_reader/XMLFileReader.h"
#include "simulator/physics/Gravity.h"
#include "utils/ArrayUtils.h"
#include <array>
#include <gtest/gtest.h>

class XMLFileReaderTest : public ::testing::Test {
 public:

  // CONSTANTS (to be changed when values in XML or input files are changed)
  const std::string base_name = "ExampleName";
  const int output_frequency = 100;
  const double end_time = 5;
  const double delta_t = 0.0002;
  const double epsilon = 5;
  const double sigma = 1;

  const std::array<int, 3> coordinates_x1 = {0, 0 , 0};
  const std::array<int, 3> particles_x1 = {40, 8 , 1};
  const std::array<int, 3> velocity_x1 = {0, 0 , 0};

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
  EXPECT_EQ(config->cuboids.size(), 2);

  int array_count = 0;

  for (const auto &elem : config->cuboids) {
    if (typeid(elem) == typeid(std::array<int, 3>)) {
      std::array<int, 3> array = (const std::array<int, 3> &) elem;
      switch (array_count) {
        case 0:
          EXPECT_TRUE(std::equal(array.begin(), array.end(), coordinates_x1.begin()));
        case 1:
          EXPECT_TRUE(std::equal(array.begin(), array.end(), particles_x1.begin()));
        case 2:
          EXPECT_TRUE(std::equal(array.begin(), array.end(), velocity_x1.begin()));
      }
      array_count++;
    }
    else {

    }
  }
}