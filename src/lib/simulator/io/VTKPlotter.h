#pragma once

#include "Plotter.h"
#include "VTKWriter.h"
#include "lib/ParticleContainer.h"
#include "lib/config/config.h"

namespace simulator::io {
class VTKPlotter final : public Plotter {
 private:
  VTKWriter vtk_writer;
  std::shared_ptr<config::Config> config;

 public:
  auto plotParticles(ParticleContainer &particle_container, int iteration) -> void override;
  ~VTKPlotter() override;
  VTKPlotter() = default;

  VTKPlotter(VTKPlotter &plotter) = delete;
  VTKPlotter(VTKPlotter &&plotter) = delete;
  explicit VTKPlotter(const std::shared_ptr<config::Config> &config);

  auto operator=(VTKPlotter &plotter) = delete;
  auto operator=(VTKPlotter &&plotter) -> VTKPlotter & = delete;
};
}// namespace simulator::io