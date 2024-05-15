#include "VTKPlotter.h"
#include "lib/outputWriter/VTKWriter.h"

namespace simulator::io {
auto VTKPlotter::plotParticles(ParticleContainer &particle_container, int iteration) -> void {

  std::string out_name("MD_vtk");

  assert(particle_container.size() <= INT_MAX);
  vtk_writer.initializeOutput(static_cast<int>(particle_container.size()));

  for (auto &particle : particle_container) {
    vtk_writer.plotParticle(particle);
  }
  vtk_writer.writeFile(out_name, iteration);
}

VTKPlotter::~VTKPlotter() = default;
}// namespace simulator::io