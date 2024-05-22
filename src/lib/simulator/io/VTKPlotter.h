#pragma once

#include "Plotter.h"
#include "VTKWriter.h"
#include "config/config.h"

namespace simulator::io {

/*!
 * plotting vtk files
 */
class VTKPlotter final : public Plotter {
 private:
  VTKWriter vtk_writer;
  std::shared_ptr<config::Config> config;

 public:
  /**
 * @brief Plots particles to a VTK file.
 *
 * Initializes the VTK writer, plots each particle, and writes the file.
 *
 * @param particle_container The container holding the particles to be plotted.
 * @param iteration The current iteration number, used in the output file name.
 */

  auto plotParticles(container::particle_container &particle_container, int iteration) -> void override {
    std::string out_name(config->output_filename);

    assert(particle_container.size() <= INT_MAX);
    vtk_writer.initializeOutput(static_cast<int>(particle_container.size()));

    for (auto &particle : particle_container) {
      vtk_writer.plotParticle(particle);
    }
    vtk_writer.writeFile(out_name, iteration);
  };

  ~VTKPlotter() override = default;
  VTKPlotter() = default;

  VTKPlotter(VTKPlotter &plotter) = delete;
  VTKPlotter(VTKPlotter &&plotter) = delete;
  explicit VTKPlotter(const std::shared_ptr<config::Config> &config){};

  auto operator=(VTKPlotter &plotter) = delete;
  auto operator=(VTKPlotter &&plotter) -> VTKPlotter & = delete;
};
}// namespace simulator::io