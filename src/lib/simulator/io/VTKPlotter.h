#pragma once

#include "Plotter.h"
#include "VTKWriter.h"
#include "config/config.h"
#include <utility>
#include <variant>

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
  auto plotParticles(container::particle_container &container, int iteration) -> void override {

    std::string out_name(config->output_filename);

    vtk_writer.initializeOutput(container.size());

    auto linear = container.linear();

    std::visit(overloads{
                   [this](std::ranges::ref_view<std::vector<Particle>> &range) {
                     for (auto &particle : range) {
                       vtk_writer.plotParticle(particle);
                     }
                   },
                   [this](auto &range) {
                     for (const auto &particle : range) {
                       vtk_writer.plotParticle(*particle);
                     }
                   }},
               linear);
    vtk_writer.writeFile(out_name, iteration);
  };

  ~VTKPlotter() override = default;
  VTKPlotter() = default;

  VTKPlotter(VTKPlotter &plotter) = delete;
  VTKPlotter(VTKPlotter &&plotter) = delete;
  explicit VTKPlotter(std::shared_ptr<config::Config> config) : config(std::move(config)) {};

  auto operator=(VTKPlotter &plotter) = delete;
  auto operator=(VTKPlotter &&plotter) -> VTKPlotter & = delete;
};
}// namespace simulator::io