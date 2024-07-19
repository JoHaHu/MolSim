#include "Particle.h"
#include "container/LinkedCell.h"
#include "container/ParticleContainer.h"
#include "simulator/physics/LennardJones.h"
#include <gtest/gtest.h>

class LinkedCellsTest : public ::testing::Test {

 public:
  std::unique_ptr<container::Container<3>> outflow_pc;
  std::unique_ptr<container::Container<3>> reflecting_pc;

  LinkedCellsTest() = default;

 protected:
  void SetUp() override {
    Test::SetUp();

    std::vector<double> epsilons = {1.0};
    std::vector<double> sigmas = {1.0};
    std::vector<double> sigma_outflow = {1.0};   // Local variable for sigma
    std::vector<double> sigma_reflecting = {1.0};// Local variable for sigma

    auto lc = std::make_unique<container::LinkedCell<simulator::physics::LennardJonesForce, 3>>(
        simulator::physics::LennardJonesForce(1.0, epsilons, sigmas),
        std::array<double, 3>{1.0, 1.0, 1.0},
        1.0,
        std::array<BoundaryCondition, 6>{
            BoundaryCondition::outflow,
            BoundaryCondition::outflow,
            BoundaryCondition::outflow,
            BoundaryCondition::outflow,
            BoundaryCondition::outflow,
            BoundaryCondition::outflow,
        },
        sigma_outflow);
    lc->insert(Particle<3>({0.5, 0.5, 0.5}, {0, 0, 0}, 1.0, 0));
    outflow_pc = std::move(lc);

    auto lc2 = std::make_unique<container::LinkedCell<simulator::physics::LennardJonesForce, 3>>(
        simulator::physics::LennardJonesForce(1.0, epsilons, sigmas),
        std::array<double, 3>{1.0, 1.0, 1.0},
        1.0,
        std::array<BoundaryCondition, 6>{
            BoundaryCondition::reflecting,
            BoundaryCondition::reflecting,
            BoundaryCondition::reflecting,
            BoundaryCondition::reflecting,
            BoundaryCondition::reflecting,
            BoundaryCondition::reflecting,
        },
        sigma_reflecting);
    lc2->insert(Particle<3>({0.5, 0.5, 0.5}, {0, 0, 0}, 1.0, 0));
    reflecting_pc = std::move(lc2);
  }
};

TEST_F(LinkedCellsTest, test_outflow_left) {
  outflow_pc->linear([](Particles<3> &particles, size_t index) {
    particles.positions[0][index] = -1.0;
  });

  outflow_pc->boundary();
  outflow_pc->refresh();

  auto counter = 0;
  outflow_pc->linear([&counter](Particles<3> &particles, size_t index) {
    counter++;
  });
  EXPECT_EQ(counter, 0);
}

TEST_F(LinkedCellsTest, test_outflow_right) {
  outflow_pc->linear([](Particles<3> &particles, size_t index) {
    particles.positions[0][index] = 2.0;
  });

  outflow_pc->boundary();
  outflow_pc->refresh();

  auto counter = 0;
  outflow_pc->linear([&counter](Particles<3> &particles, size_t index) {
    counter++;
  });
  EXPECT_EQ(counter, 0);
}

TEST_F(LinkedCellsTest, test_outflow_top) {
  outflow_pc->linear([](Particles<3> &particles, size_t index) {
    particles.positions[1][index] = 1.5;
  });

  outflow_pc->boundary();
  outflow_pc->refresh();

  auto counter = 0;
  outflow_pc->linear([&counter](Particles<3> &particles, size_t index) {
    counter++;
  });
  EXPECT_EQ(counter, 0);
}

TEST_F(LinkedCellsTest, test_outflow_bottom) {
  outflow_pc->linear([](Particles<3> &particles, size_t index) {
    particles.positions[1][index] = -0.5;
  });

  outflow_pc->boundary();
  outflow_pc->refresh();

  auto counter = 0;
  outflow_pc->linear([&counter](Particles<3> &particles, size_t index) {
    counter++;
  });
  EXPECT_EQ(counter, 0);
}

TEST_F(LinkedCellsTest, test_outflow_front) {
  outflow_pc->linear([](Particles<3> &particles, size_t index) {
    particles.positions[2][index] = 1.5;
  });

  outflow_pc->boundary();
  outflow_pc->refresh();

  auto counter = 0;
  outflow_pc->linear([&counter](Particles<3> &particles, size_t index) {
    counter++;
  });
  EXPECT_EQ(counter, 0);
}

TEST_F(LinkedCellsTest, test_outflow_back) {
  outflow_pc->linear([](Particles<3> &particles, size_t index) {
    particles.positions[2][index] = -1.5;
  });

  outflow_pc->boundary();
  outflow_pc->refresh();

  auto counter = 0;
  outflow_pc->linear([&counter](Particles<3> &particles, size_t index) {
    counter++;
  });
  EXPECT_EQ(counter, 0);
}

/*#### REFLECTING CONDITIONS #####*/

TEST_F(LinkedCellsTest, test_reflecting_left) {
  reflecting_pc->linear([](Particles<3> &particles, size_t index) {
    particles.positions[0][index] = 0.25;
  });

  reflecting_pc->boundary();
  reflecting_pc->linear([](Particles<3> &particles, size_t index) {
    EXPECT_NEAR(particles.positions[0][index], 0.25, 1e-5);
  });
  reflecting_pc->refresh();
}

TEST_F(LinkedCellsTest, test_reflecting_right) {
  reflecting_pc->linear([](Particles<3> &particles, size_t index) {
    particles.positions[0][index] = 0.75;
  });

  reflecting_pc->boundary();
  reflecting_pc->linear([](Particles<3> &particles, size_t index) {
    EXPECT_NEAR(particles.positions[0][index], 0.75, 1e-5);
  });
}

TEST_F(LinkedCellsTest, test_reflecting_bottom) {
  reflecting_pc->linear([](Particles<3> &particles, size_t index) {
    particles.positions[1][index] = 0.25;
  });

  reflecting_pc->boundary();
  reflecting_pc->linear([](Particles<3> &particles, size_t index) {
    EXPECT_NEAR(particles.positions[1][index], 0.25, 1e-5);
  });
  reflecting_pc->refresh();
}

TEST_F(LinkedCellsTest, test_reflecting_top) {
  reflecting_pc->linear([](Particles<3> &particles, size_t index) {
    particles.positions[1][index] = 0.75;
  });

  reflecting_pc->boundary();
  reflecting_pc->linear([](Particles<3> &particles, size_t index) {
    EXPECT_NEAR(particles.positions[1][index], 0.75, 1e-5);
  });
}

TEST_F(LinkedCellsTest, test_reflecting_back) {
  reflecting_pc->linear([](Particles<3> &particles, size_t index) {
    particles.positions[2][index] = 0.25;
  });

  reflecting_pc->boundary();
  reflecting_pc->linear([](Particles<3> &particles, size_t index) {
    EXPECT_NEAR(particles.positions[2][index], 0.25, 1e-5);
  });
  reflecting_pc->refresh();
}

TEST_F(LinkedCellsTest, test_reflecting_front) {
  reflecting_pc->linear([](Particles<3> &particles, size_t index) {
    particles.positions[2][index] = 0.75;
  });

  reflecting_pc->boundary();
  reflecting_pc->linear([](Particles<3> &particles, size_t index) {
    EXPECT_NEAR(particles.positions[2][index], 0.75, 1e-5);
  });
}
