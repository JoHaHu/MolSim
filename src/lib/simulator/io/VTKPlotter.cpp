#include "simulator/io/VTKPlotter.h"
#include "simulator/io/VTKWriter.h"

namespace simulator::io {
/**
 * @brief Plots particles to a VTK file.
 *
 * Initializes the VTK writer, plots each particle, and writes the file.
 *
 * @param particle_container The container holding the particles to be plotted.
 * @param iteration The current iteration number, used in the output file name.
 */
auto VTKPlotter::plotParticles(ParticleContainer &particle_container, int iteration) -> void {

  std::string out_name(config->output_filename);

  assert(particle_container.size() <= INT_MAX);
  vtk_writer.initializeOutput(static_cast<int>(particle_container.size()));

  for (auto &particle : particle_container) {
    vtk_writer.plotParticle(particle);
  }
  vtk_writer.writeFile(out_name, iteration);
}

VTKPlotter::VTKPlotter(const std::shared_ptr<OldConfig::OldConfig> &config) : config(config) {}

VTKPlotter::~VTKPlotter() = default;
}// namespace simulator::io