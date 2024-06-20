#include "container/index.h"
#include <gtest/gtest.h>

class IndexTest : public ::testing::Test {
 public:
  container::index::Index<2> index{};

  container::index::Index<2> mixed_index{};

 protected:
  void SetUp() override {
    index = container::index::Index<2>(
        std::array<double, 2>({3.0, 3.0}),
        {BoundaryCondition::periodic,
         BoundaryCondition::periodic,
         BoundaryCondition::periodic,
         BoundaryCondition::periodic},
        1.0);
    mixed_index = container::index::Index<2>(
        std::array<double, 2>({3.0, 3.0}),
        {BoundaryCondition::periodic,
         BoundaryCondition::reflecting,
         BoundaryCondition::periodic,
         BoundaryCondition::reflecting},
        1.0);
  };
};

TEST_F(IndexTest, test_corrections) {
  EXPECT_EQ(index.dim[0], 3);
  EXPECT_EQ(index.dim[1], 3);

  auto correction = std::array<double, 2>({-3, -3});
  EXPECT_EQ(index.calculate_correction({0, 0}, {-1, -1}), correction);

  auto correction2 = std::array<double, 2>({3, 3});
  EXPECT_EQ(index.calculate_correction({2, 2}, {1, 1}), correction2);
}

TEST_F(IndexTest, test_to_index) {
  EXPECT_EQ(index.dim[0], 3);
  EXPECT_EQ(index.dim[1], 3);

  auto offset_index = 7;
  EXPECT_EQ(index.dimension_to_index({1, 2}), offset_index);

  offset_index = 0;
  EXPECT_EQ(index.dimension_to_index({0, 0}), offset_index);
  offset_index = 8;
  EXPECT_EQ(index.dimension_to_index({2, 2}), offset_index);

  offset_index = 8;
  EXPECT_EQ(index.dimension_to_index({5, 5}), offset_index);
}

TEST_F(IndexTest, test_offset) {
  EXPECT_EQ(index.dim[0], 3);
  EXPECT_EQ(index.dim[1], 3);

  auto offset_index = 7;
  EXPECT_EQ(index.offset({0UL, 0UL}, {1, -1}), offset_index);

  offset_index = 8;
  EXPECT_EQ(index.offset({0UL, 0UL}, {-1, -1}), offset_index);

  offset_index = 3;
  EXPECT_EQ(index.offset({2UL, 2UL}, {1, -1}), offset_index);

  offset_index = 6;
  EXPECT_EQ(index.offset({2UL, 2UL}, {1, 0}), offset_index);

  offset_index = 2;
  EXPECT_EQ(index.offset({2UL, 2UL}, {0, 1}), offset_index);

  offset_index = 0;
  EXPECT_EQ(index.offset({2UL, 2UL}, {1, 1}), offset_index);
}

TEST_F(IndexTest, test_position_to_index) {
  EXPECT_EQ(index.dim[0], 3);
  EXPECT_EQ(index.dim[1], 3);

  auto idx = 0;
  EXPECT_EQ(index.position_to_index({0, 0}), idx);

  idx = 0;
  EXPECT_EQ(index.position_to_index({3, 3}), idx);

  idx = 8;
  EXPECT_EQ(index.position_to_index({2, 2}), idx);

  idx = 0;
  EXPECT_EQ(index.position_to_index({-1, -1}), idx);
  idx = 8;
  EXPECT_EQ(index.position_to_index({-1.5, -1.5}), idx);
}

TEST_F(IndexTest, test_mixed_to_index) {
  EXPECT_EQ(mixed_index.dim[0], 3);
  EXPECT_EQ(mixed_index.dim[1], 3);

  auto offset_index = 7;
  EXPECT_EQ(mixed_index.dimension_to_index({1, 2}), offset_index);

  offset_index = 0;
  EXPECT_EQ(mixed_index.dimension_to_index({0, 0}), offset_index);
  offset_index = 8;
  EXPECT_EQ(mixed_index.dimension_to_index({2, 2}), offset_index);

  offset_index = ULONG_MAX;
  EXPECT_EQ(mixed_index.dimension_to_index({5, 5}), offset_index);
  offset_index = 8;
  EXPECT_EQ(mixed_index.dimension_to_index({5, 2}), offset_index);
}

TEST_F(IndexTest, test_mixed_corrections) {
  EXPECT_EQ(index.dim[0], 3);
  EXPECT_EQ(index.dim[1], 3);

  auto correction = std::array<double, 2>({-3, 0});
  EXPECT_EQ(mixed_index.calculate_correction({0, 0}, {-1, -1}), correction);

  correction = std::array<double, 2>({3, 0});
  EXPECT_EQ(mixed_index.calculate_correction({2, 2}, {1, 1}), correction);

  correction = std::array<double, 2>({0, 0});
  EXPECT_EQ(mixed_index.calculate_correction({1, 1}, {1, 1}), correction);
}