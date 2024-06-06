#include "Particle.h"
#include "container/boundary.h"
#include "container/container.h"
#include <gtest/gtest.h>

class LinkedCellsTest : public ::testing::Test {

 public:
  container::particle_container outflow_pc;
  container::particle_container reflecting_pc;
  LinkedCellsTest() : outflow_pc(container::particle_container(std::vector<Particle>())),
                      reflecting_pc(container::particle_container(std::vector<Particle>())){};

 protected:
  void SetUp() override {
    Test::SetUp();

    auto lc = container::linked_cell<container::index::row_major_index>(
        {1.0, 1.0, 1.0},
        1.0,
        {
            BoundaryCondition::outflow,
            BoundaryCondition::outflow,
            BoundaryCondition::outflow,
            BoundaryCondition::outflow,
            BoundaryCondition::outflow,
            BoundaryCondition::outflow,
        },
        1,
        1.0);
    lc.insert(Particle({0.5, 0.5, 0.5}, {0, 0, 0}, 0, 0));
    outflow_pc = container::particle_container(std::move(lc));
    auto lc2 = container::linked_cell<container::index::row_major_index>(
        {1.0, 1.0, 1.0},
        1.0,
        {
            BoundaryCondition::reflecting,
            BoundaryCondition::reflecting,
            BoundaryCondition::reflecting,
            BoundaryCondition::reflecting,
            BoundaryCondition::reflecting,
            BoundaryCondition::reflecting,
        },
        1,
        0.25);
    lc2.insert(Particle({0.5, 0.5, 0.5}, {0, 0, 0}, 0, 0));
    reflecting_pc = container::particle_container(std::move(lc2));
  }
};

TEST_F(LinkedCellsTest, test_outflow_left) {

  outflow_pc.linear([](Particle &p) {
    p.position = {-1.0, 0.5, 0.5};
  });

  outflow_pc.boundary([](auto p) {});
  outflow_pc.refresh();

  auto counter = 0;
  outflow_pc.linear([&counter](auto &p) {
    counter++;
  });
  EXPECT_EQ(counter, 0);
}
TEST_F(LinkedCellsTest, test_outflow_rigth) {

  outflow_pc.linear([](Particle &p) {
    p.position = {2.0, 0.5, 0.5};
  });

  outflow_pc.boundary([](auto p) {});
  outflow_pc.refresh();

  auto counter = 0;
  outflow_pc.linear([&counter](auto &p) {
    counter++;
  });
  EXPECT_EQ(counter, 0);
}
TEST_F(LinkedCellsTest, test_outflow_top) {

  outflow_pc.linear([](Particle &p) {
    p.position = {0.5, 1.5, 0.5};
  });

  outflow_pc.boundary([](auto p) {});
  outflow_pc.refresh();

  auto counter = 0;
  outflow_pc.linear([&counter](auto &p) {
    counter++;
  });
  EXPECT_EQ(counter, 0);
}
TEST_F(LinkedCellsTest, test_outflow_bottom) {

  outflow_pc.linear([](Particle &p) {
    p.position = {0.5, -0.5, 0.5};
  });

  outflow_pc.boundary([](auto p) {});
  outflow_pc.refresh();

  auto counter = 0;
  outflow_pc.linear([&counter](auto &p) {
    counter++;
  });
  EXPECT_EQ(counter, 0);
}

TEST_F(LinkedCellsTest, test_outflow_front) {

  outflow_pc.linear([](Particle &p) {
    p.position = {0.5, 0.5, 1.5};
  });

  outflow_pc.boundary([](auto p) {});
  outflow_pc.refresh();

  auto counter = 0;
  outflow_pc.linear([&counter](auto &p) {
    counter++;
  });
  EXPECT_EQ(counter, 0);
}

TEST_F(LinkedCellsTest, test_outflow_back) {

  outflow_pc.linear([](Particle &p) {
    p.position = {0.5, 0.5, -1.5};
  });

  outflow_pc.boundary([](auto p) {});
  outflow_pc.refresh();

  auto counter = 0;
  outflow_pc.linear([&counter](auto &p) {
    counter++;
  });
  EXPECT_EQ(counter, 0);
}

/*#### REFLECTING CONDITIONS #####*/

TEST_F(LinkedCellsTest, test_reflecting_left) {

  reflecting_pc.linear([](Particle &p) {
    p.position = {0.25, 0.5, 0.5};
  });

  reflecting_pc.boundary([](std::tuple<Particle &, Particle &> p) {
    auto [p1, p2] = p;
    EXPECT_EQ(p2.position[0], -0.25);
  });
  reflecting_pc.refresh();
}

TEST_F(LinkedCellsTest, test_reflecting_right) {

  reflecting_pc.linear([](Particle &p) {
    p.position = {0.75, 0.5, 0.5};
  });

  reflecting_pc.boundary([](std::tuple<Particle &, Particle &> p) {
    auto [p1, p2] = p;
    EXPECT_EQ(p2.position[0], 1.25);
  });
}

TEST_F(LinkedCellsTest, test_reflecting_bottom) {

  reflecting_pc.linear([](Particle &p) {
    p.position = {0.5, 0.25, 0.5};
  });

  reflecting_pc.boundary([](std::tuple<Particle &, Particle &> p) {
    auto [p1, p2] = p;
    EXPECT_EQ(p2.position[1], -0.25);
  });
  reflecting_pc.refresh();
}

TEST_F(LinkedCellsTest, test_reflecting_top) {

  reflecting_pc.linear([](Particle &p) {
    p.position = {0.5, 0.75, 0.5};
  });

  reflecting_pc.boundary([](std::tuple<Particle &, Particle &> p) {
    auto [p1, p2] = p;
    EXPECT_EQ(p2.position[1], 1.25);
  });
}

TEST_F(LinkedCellsTest, test_reflecting_back) {

  reflecting_pc.linear([](Particle &p) {
    p.position = {0.5, 0.5, 0.25};
  });

  reflecting_pc.boundary([](std::tuple<Particle &, Particle &> p) {
    auto [p1, p2] = p;
    EXPECT_EQ(p2.position[2], -0.25);
  });
  reflecting_pc.refresh();
}

TEST_F(LinkedCellsTest, test_reflecting_front) {

  reflecting_pc.linear([](Particle &p) {
    p.position = {0.5, 0.5, 0.75};
  });

  reflecting_pc.boundary([](std::tuple<Particle &, Particle &> p) {
    auto [p1, p2] = p;
    EXPECT_EQ(p2.position[2], 1.25);
  });
}
