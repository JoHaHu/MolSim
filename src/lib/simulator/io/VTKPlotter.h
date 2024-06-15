#pragma once

#include "Plotter.h"
#include "VTKWriter.h"
#include "config/Config.h"
#include <utility>
#include <variant>

namespace simulator::io {

/*!
 * plotting vtk files
 */
template<const size_t DIMENSIONS>
class VTKPlotter final : public Plotter<DIMENSIONS> {
 private:
  VTKWriter<DIMENSIONS> vtk_writer;
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
  auto plotParticles(container::ParticleContainer<DIMENSIONS> &particles, int iteration) -> void override {

    std::string out_name(config->output_filename);
    vtk_writer.initializeOutput(particles.size());
    particles.linear([this](Particles<DIMENSIONS> &particles, size_t index) { vtk_writer.plotParticle(particles, index); });
    vtk_writer.writeFile(out_name, iteration);
  };

  ~VTKPlotter() override = default;
  VTKPlotter() = default;

  VTKPlotter(VTKPlotter &plotter) = delete;
  VTKPlotter(VTKPlotter &&plotter) = delete;
  explicit VTKPlotter(std::shared_ptr<config::Config> config) : config(std::move(config)){};

  auto operator=(VTKPlotter &plotter) = delete;
  auto operator=(VTKPlotter &&plotter) -> VTKPlotter & = delete;
};
}// namespace simulator::io