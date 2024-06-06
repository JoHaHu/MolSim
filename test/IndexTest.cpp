#include "container/index.h"
#include <gtest/gtest.h>

class SimpleIndexTest : public ::testing::Test {
};

/*
 * Test that the index in row major ordering behaves okay
 * */
TEST_F(SimpleIndexTest, row_major_ordering_index) {

  auto bounds = std::array<double, 3>({3, 3, 3});
  auto idx = container::index::row_major_index(bounds, 1);
  EXPECT_EQ(idx.boundary, bounds);
  auto dim = std::array<size_t, 3>({3, 3, 3});
  EXPECT_EQ(idx.dimension, dim);
  auto width = std::array<double, 3>({1, 1, 1});
  EXPECT_EQ(idx.width, width);
  auto radius = std::array<size_t, 3>({1, 1, 1});
  EXPECT_EQ(idx.radius, radius);

  for (size_t x = 0; x < 3; ++x) {
    for (size_t y = 0; y < 3; ++y) {
      for (size_t z = 0; z < 3; ++z) {
        auto new_index = idx.dimension_to_index({x, y, z});
        EXPECT_LT(new_index, 27);
      }
    }
  }

  auto position = std::array<double, 3>({1.5, 1.5, 1.5});
  auto cell = idx.dimension_to_index({1, 1, 1});
  EXPECT_EQ(idx.position_to_index(position), cell);
}

class ZOrderIndexTest : public ::testing::Test {
};

//TEST_F(SimpleIndexTest, z_order) {
//}