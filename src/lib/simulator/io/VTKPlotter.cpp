//
// Created by johannes on 14.05.24.
//

#include "VTKPlotter.h"
#include "lib/outputWriter/VTKWriter.h"
#include "lib/outputWriter/XYZWriter.h"

auto VTKPlotter::plotParticles(ParticleContainer &particles, int iteration) -> void {

  std::string out_name("MD_vtk");

  outputWriter::XYZWriter::plotParticles(particles, out_name, iteration);

  outputWriter::VTKWriter vtk_writer;
  assert(particles.size() <= INT_MAX);
  vtk_writer.initializeOutput(static_cast<int>(particles.size()));

  for (auto &p : particles) {
    vtk_writer.plotParticle(p);
  }
  vtk_writer.writeFile(out_name, iteration);
}
