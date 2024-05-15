#pragma once

#include "lib/Particle.h"
#include "lib/simulator/io/VTKUnstructured.hxx"

#include <list>

namespace simulator::io {

/**
 * This class implements the functionality to generate vtk output from
 * particles.
 */
class VTKWriter {

 public:
  VTKWriter() = default;

  virtual ~VTKWriter() = default;
  /**
   * set up internal data structures and prepare to plot a particle.
   */
  void initializeOutput(int numParticles);

  /**
   * plot type, mass, position, velocity and force of a particle.
   *
   * @note: initializeOutput() must have been called before.
   */
  void plotParticle(Particle &particle);

  /**
   * writes the final output file.
   *
   * @param filename the base name of the file to be written.
   * @param iteration the number of the current iteration,
   *        which is used to generate an unique filename
   */
  void writeFile(const std::string &filename, int iteration);

 private:
  std::unique_ptr<VTKFile_t> vtkFile = std::make_unique<VTKFile_t>(VTKFile_t("UnstructuredGrid"));
};

}// namespace simulator::io