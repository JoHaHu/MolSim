//#include "Particle.h"
#include "container/index.h"
#include <gtest/gtest.h>

class SimpleIndexTest : public ::testing::Test {
};

TEST_F(SimpleIndexTest, index_conversions_including_boundaries) {

  auto bounds = std::array<double, 3>({3, 3, 3});

  auto idx = container::index::row_major_index(bounds, 1);
}

class ZOrderIndexTest : public ::testing::Test {
};

TEST_F(SimpleIndexTest, z_order) {
}