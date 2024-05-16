//
// Created by TimSc on 16.05.2024.
//

#include <gtest/gtest.h>

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-avoid-magic-numbers"



class LennardJonesPotentialTest : public ::testing::Test {
 public:

  void SetUp() override {


  }
};

TEST_F(LennardJonesPotentialTest, basic_test) {
  SetUp();

  EXPECT_EQ(1, 1);
}