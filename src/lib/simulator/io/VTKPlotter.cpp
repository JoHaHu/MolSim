#include "lib/simulator/io/VTKPlotter.h"
#include "lib/simulator/io/VTKWriter.h"

namespace simulator::io {
auto VTKPlotter::plotParticles(ParticleContainer &particle_container, int iteration) -> void {

  std::string out_name(config->output_filename);

  assert(particle_container.size() <= INT_MAX);
  vtk_writer.initializeOutput(static_cast<int>(particle_container.size()));

  for (auto &particle : particle_container) {
    vtk_writer.plotParticle(particle);
  }
  vtk_writer.writeFile(out_name, iteration);
}

VTKPlotter::VTKPlotter(const std::shared_ptr<config::Config> &config) : config(config) {}

VTKPlotter::~VTKPlotter() = default;
}// namespace simulator::io