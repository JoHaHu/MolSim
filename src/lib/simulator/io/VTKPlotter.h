#pragma once

#include "Plotter.h"
#include "VTKWriter.h"
#include "lib/ParticleContainer.h"

namespace simulator::io {
class VTKPlotter final : public Plotter {
 private:
  VTKWriter vtk_writer;

 public:
  auto plotParticles(ParticleContainer &particle_container, int iteration) -> void override;
  ~VTKPlotter() override;
  VTKPlotter() = default;

  VTKPlotter(VTKPlotter &plotter) = delete;
  VTKPlotter(VTKPlotter &&plotter) = default;

  auto operator=(VTKPlotter &plotter) = delete;
  auto operator=(VTKPlotter &&plotter) -> VTKPlotter & = default;
};
}// namespace simulator::io