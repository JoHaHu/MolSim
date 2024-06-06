#include "Particle.h"
#include "container/container.h"
#include <gtest/gtest.h>

class LinkedCellsTest : public ::testing::Test {

 public:
  container::particle_container pc;
  LinkedCellsTest() : pc(container::particle_container(std::vector<Particle>())) {};

 protected:
  void SetUp() override {
    Test::SetUp();

    auto lc = container::linked_cell<container::index::row_major_index>(
        {1.0, 1.0, 1.0},
        1.0,
        {
            container::boundary_condition::outflow,
            container::boundary_condition::outflow,
            container::boundary_condition::outflow,
            container::boundary_condition::outflow,
            container::boundary_condition::outflow,
            container::boundary_condition::outflow,
        },
        1,
        1.0);
    lc.insert(Particle({0.5, 0.5, 0.5}, {0, 0, 0}, 0, 0));
    pc = container::particle_container(std::move(lc));
  }
};

TEST_F(LinkedCellsTest, test_outflow_left) {

  pc.linear([](Particle &p) {
    p.position = {-1.0, 0.5, 0.5};
  });

  pc.boundary([](auto p) {});
  pc.refresh();

  auto counter = 0;
  pc.linear([&counter](auto &p) {
    counter++;
  });
  EXPECT_EQ(counter, 0);
}
TEST_F(LinkedCellsTest, test_outflow_rigth) {

  pc.linear([](Particle &p) {
    p.position = {2.0, 0.5, 0.5};
  });

  pc.boundary([](auto p) {});
  pc.refresh();

  auto counter = 0;
  pc.linear([&counter](auto &p) {
    counter++;
  });
  EXPECT_EQ(counter, 0);
}
TEST_F(LinkedCellsTest, test_outflow_top) {

  pc.linear([](Particle &p) {
    p.position = {0.5, 1.5, 0.5};
  });

  pc.boundary([](auto p) {});
  pc.refresh();

  auto counter = 0;
  pc.linear([&counter](auto &p) {
    counter++;
  });
  EXPECT_EQ(counter, 0);
}
TEST_F(LinkedCellsTest, test_outflow_bottom) {

  pc.linear([](Particle &p) {
    p.position = {0.5, -0.5, 0.5};
  });

  pc.boundary([](auto p) {});
  pc.refresh();

  auto counter = 0;
  pc.linear([&counter](auto &p) {
    counter++;
  });
  EXPECT_EQ(counter, 0);
}

TEST_F(LinkedCellsTest, test_outflow_front) {

  pc.linear([](Particle &p) {
    p.position = {0.5, 0.5, 1.5};
  });

  pc.boundary([](auto p) {});
  pc.refresh();

  auto counter = 0;
  pc.linear([&counter](auto &p) {
    counter++;
  });
  EXPECT_EQ(counter, 0);
}

TEST_F(LinkedCellsTest, test_outflow_back) {

  pc.linear([](Particle &p) {
    p.position = {0.5, 0.5, -1.5};
  });

  pc.boundary([](auto p) {});
  pc.refresh();

  auto counter = 0;
  pc.linear([&counter](auto &p) {
    counter++;
  });
  EXPECT_EQ(counter, 0);
}