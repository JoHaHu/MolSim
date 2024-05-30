//#include "Particle.h"
#include "container/index.h"
#include <gtest/gtest.h>

class IndexTest : public ::testing::Test {
};

TEST_F(IndexTest, index_conversions_including_boundaries) {

  auto bounds = std::array<double, 3>({3, 3, 3});

  auto idx = container::index::simple_index(bounds, 1);

  EXPECT_EQ(idx.size(), 27);
  auto test = idx.index_to_dimension(2 + 1 * 5 + 2 * 5 * 5);
  auto res = std::array<size_t, 3>({1, 0, 1});
  EXPECT_EQ(test, res);

  auto test2 = idx.index_to_dimension(3 + 3 * 5 + 3 * 5 * 5);
  auto res2 = std::array<size_t, 3>({2, 2, 2});
  EXPECT_EQ(test2, res2);

  auto test3 = idx.index_to_dimension(1 + 1 * 5 + 1 * 5 * 5);
  auto res3 = std::array<size_t, 3>({0, 0, 0});
  EXPECT_EQ(test3, res3);

  auto test4 = idx.index_to_dimension(1 + 2 * 5 + 3 * 5 * 5);
  auto res4 = std::array<size_t, 3>({0, 1, 2});
  EXPECT_EQ(test4, res4);
  auto test5 = idx.index_to_dimension(3 + 2 * 5 + 1 * 5 * 5);
  auto res5 = std::array<size_t, 3>({2, 1, 0});
  EXPECT_EQ(test5, res5);
}