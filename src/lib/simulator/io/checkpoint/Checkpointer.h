#pragma once

#include "Checkpoint.hxx"
#include "Particle.h"
#include "container/Container.h"
#include "fstream"
#include "iostream"
#include "spdlog/spdlog.h"
#include "xercesc/dom/DOM.hpp"
#include "xercesc/util/PlatformUtils.hpp"
#include <cstdlib>
#include <string>

using namespace xercesc;

template<const size_t DIMENSIONS>
class Checkpointer {
 public:
  auto load_checkpoint(std::string &filename) -> std::vector<Particle<DIMENSIONS>> {

    try {
      XMLPlatformUtils::Initialize();
    } catch (const XMLException &exception) {
      char *message = XMLString::transcode(exception.getMessage());
      SPDLOG_WARN("Error during XML file initialization!\n");
      SPDLOG_INFO(message);
      exit(1);
    }
    std::unique_ptr<checkpoint> cp = checkpoint_(filename, xml_schema::flags::dont_validate);

    auto particles = std::vector<Particle<DIMENSIONS>>();

    for (auto &p : cp->particle()) {

      std::array<double, DIMENSIONS> force;
      std::array<double, DIMENSIONS> old_force;
      std::array<double, DIMENSIONS> position;
      std::array<double, DIMENSIONS> velocity;

      for (size_t i = 0; i < DIMENSIONS; ++i) {
        switch (i) {
          case 0: {
            force[i] = p.force().x();
            old_force[i] = p.old_force().x();
            position[i] = p.position().x();
            velocity[i] = p.velocity().x();
            break;
          }
          case 1: {
            force[i] = p.force().y();
            old_force[i] = p.old_force().y();
            position[i] = p.position().y();
            velocity[i] = p.velocity().y();
            break;
          }
          case 2: {
            if (p.force().z().present()) {

              force[i] = *p.force().z();
              old_force[i] = *p.old_force().z();
              position[i] = *p.position().z();
              velocity[i] = *p.velocity().z();
            } else {
              force[i] = 0.0;
              force[i] = 0.0;
              force[i] = 0.0;
              force[i] = 0.0;
            }
            break;
          }
          default: {
            SPDLOG_ERROR(" Invalid checkpoint");
            exit(1);
          }
        }
      }
      particles.emplace_back(position, velocity, force, old_force, p.mass(), p.type(), p.fixed());
    }

    return particles;
  }

  auto save_checkpoint(std::string &filename, container::Container<DIMENSIONS> &particles) {
    checkpoint cp = checkpoint();

    particles.linear([&](Particles<DIMENSIONS> &p, size_t index) {
      if constexpr (DIMENSIONS == 2) {
        auto f = particle::force_type(p.forces[0][index], p.forces[1][index]);
        auto of = particle::old_force_type(p.old_forces[0][index], p.old_forces[1][index]);
        auto v = particle::velocity_type(p.velocities[0][index], p.velocities[1][index]);
        auto pos = particle::position_type(p.positions[0][index], p.positions[1][index]);
        auto m = particle::mass_type(p.mass[index]);
        auto type = particle::type_type(p.type[index]);
        auto fixed = particle::fixed_type(p.fixed[index]);
        particle pa = particle(v, f, of, pos, m, type, fixed);
        cp.particle().push_back(pa);
      } else {
        auto f = particle::force_type(p.forces[0][index], p.forces[1][index]);
        f.z(p.forces[2][index]);
        auto of = particle::old_force_type(p.old_forces[0][index], p.old_forces[1][index]);
        of.z(p.old_forces[2][index]);
        auto v = particle::velocity_type(p.velocities[0][index], p.velocities[1][index]);
        v.z(p.velocities[2][index]);
        auto pos = particle::position_type(p.positions[0][index], p.positions[1][index]);
        pos.z(p.positions[2][index]);
        auto m = particle::mass_type(p.mass[index]);
        auto type = particle::type_type(p.type[index]);
        auto fixed = particle::fixed_type(p.fixed[index]);
        particle pa = particle(v, f, of, pos, m, type, fixed);
        cp.particle().push_back(pa);
      }
    });

    std::ofstream file(filename.c_str());

    checkpoint_(file, cp);
  }
};