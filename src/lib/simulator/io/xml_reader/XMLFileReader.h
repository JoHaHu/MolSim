//
// Created by TimSc on 02.06.2024.
//

#include "../../../config/Config.h"
#include <iostream>
#include <memory>
#include <string>

class XMLFileReader {

 public:
  /**
   * function to parse the data from the given XML using the XSD library and logging the values to the console
   * for the user to confirm. The data is also saved here for further calculation and simulation
   */
  static auto parseXMLData(const std::string &xmlFilePath) -> std::shared_ptr<Config>;
};
