#include "Particle.h"
#include "container/combination.h"
#include <algorithm>
#include <gtest/gtest.h>
#include <ranges>

class CombinationIteratorTest : public ::testing::Test {
};

// Tests that a combination iterator iterates over all unique permutations of a given range
TEST_F(CombinationIteratorTest, correct_order_of_combination) {

  auto particle = std::vector({1, 2, 3, 4});

  auto combs = particle | container::combination;

  auto combs_start = combs.cbegin();
  auto combs_end = combs.cend();

  EXPECT_NE(combs_start, combs_end);
  EXPECT_EQ(*combs_start, std::make_pair(1, 2));
  EXPECT_EQ(combs_end - combs_start, 6);
  combs_start++;

  EXPECT_NE(combs_start, combs_end);
  EXPECT_EQ(*combs_start, std::make_pair(1, 3));
  EXPECT_EQ(combs_end - combs_start, 5);
  combs_start++;

  EXPECT_NE(combs_start, combs_end);
  EXPECT_EQ(*combs_start, std::make_pair(1, 4));
  EXPECT_EQ(combs_end - combs_start, 4);
  combs_start++;

  EXPECT_NE(combs_start, combs_end);
  EXPECT_EQ(*combs_start, std::make_pair(2, 3));
  EXPECT_EQ(combs_end - combs_start, 3);
  combs_start++;

  EXPECT_NE(combs_start, combs_end);
  EXPECT_EQ(*combs_start, std::make_pair(2, 4));
  EXPECT_EQ(combs_end - combs_start, 2);
  combs_start++;

  EXPECT_NE(combs_start, combs_end);
  EXPECT_EQ(*combs_start, std::make_pair(3, 4));
  EXPECT_EQ(combs_end - combs_start, 1);
  combs_start++;

  EXPECT_EQ(combs_start, combs_end);
  EXPECT_EQ(combs_end - combs_start, 0);
}

TEST_F(CombinationIteratorTest, correct_number_of_combination) {
  size_t n = 1024;
  auto vec = std::vector<int>(n);
  auto combs = vec
      | container::combination
      | std::views::transform([](auto) { return 1; });
  EXPECT_EQ(combs.size(), (n * (n - 1)) / 2);
  auto sum = std::ranges::fold_left(combs, 0, std::plus<>());
  EXPECT_EQ(sum, (n * (n - 1)) / 2);
}