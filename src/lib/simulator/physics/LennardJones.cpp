#include "LennardJones.h"

#include <array>
#include <cmath>
#include <limits>

#include "lib/utils/ArrayUtils.h"
#include "spdlog/spdlog.h"

namespace simulator::physics {
    /**
   * @brief Calculates the Lennard Jones force between two particles.
   *
   * Computes the Lennard Jones forces based on the positions of two particles.
   *
   * @param particle1 The first particle.
   * @param particle2 The second particle.
   * @return std::array<double, 3> The calculated force vector.
   */
    std::array<double, 3> LennardJones::calculate_force(const Particle &particle1, const Particle &particle2) {
        const auto epsilon = 5;
        const auto sigma = 1;

        spdlog::trace("Entering LennardJones calculate_force");

        spdlog::trace("Particle 1: mass = {}, position = ({}, {}, {})", particle1.mass, particle1.position[0], particle1.position[1], particle1.position[2]);
        spdlog::trace("Particle 2: mass = {}, position = ({}, {}, {})", particle2.mass, particle2.position[0], particle2.position[1], particle2.position[2]);

        const auto x_diff = particle2.position - particle1.position;

        spdlog::trace("Position difference: ({}, {}, {})", x_diff[0], x_diff[1], x_diff[2]);

        auto norm = ArrayUtils::L2Norm(x_diff);
        spdlog::trace("Norm of position difference: {}", norm);

        if (norm == 0) {
            spdlog::warn("Zero distance between particles encountered");
        }

        const auto sigmaOverNorm6 = (sigma/norm) * (sigma/norm) * (sigma/norm) * (sigma/norm) * (sigma/norm) * (sigma/norm);
        const auto sigmaOverNorm12 = sigmaOverNorm6 * sigmaOverNorm6;

        const auto f = 24 * epsilon / norm * norm * (sigmaOverNorm6 - 2 * sigmaOverNorm12) * x_diff;

        spdlog::trace("Calculated force: ({}, {}, {})", f[0], f[1], f[2]);

        spdlog::trace("Exiting LennardJones calculate_force");

        return f;
    }
} // namespace simulator::physics
