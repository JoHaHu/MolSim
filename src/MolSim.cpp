
#include "FileReader.h"
#include "outputWriter/VTKWriter.h"
#include "outputWriter/XYZWriter.h"
#include "utils/ArrayUtils.h"

#include <iostream>

/**** forward declaration of the calculation functions ****/

/**
 * calculates the force for all particles
 */
void calculateF();

/**
 * calculates the position for all particles
 */
void calculateX();

/**
 * calculates the position for all particles
 */
void calculateV();

//! plotting function that takes an integer value and has no return value
/*!
  \param iteration an integer argument that sets the number of iterations
*/
void plotParticles(int iteration);

constexpr double start_time = 0; /**< initialisation of the start time of the simulation with 0 */
constexpr double end_time = 1000; /**< initialisation of the end time of the simulation with 1000 according to the worksheet task 4.1 */
constexpr double delta_t = 0.014; /**< time delta of 0.014 according to the worksheet task 4.1 */

std::vector<Particle> particles;

//! main function for execution
/*!
  \param argc an integer argument, standard for main function
  \param argsv a char array argument, standard for main function
  \return the return code
*/
int main(int argc, char *argsv[]) {

  std::cout << "Hello from MolSim for PSE!" << std::endl;
  if (argc != 2) {
    std::cout << "Erroneous programme call! " << std::endl;
    std::cout << "./molsym filename" << std::endl;
  }

  FileReader fileReader;
  fileReader.readFile(particles, argsv[1]);

  double current_time = start_time;

  int iteration = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {
    // calculate new x
    calculateX();
    // calculate new f
    calculateF();
    // calculate new v
    calculateV();

    iteration++;
    if (iteration % 10 == 0) {
      plotParticles(iteration);
    }

    std::cout << "Iteration " << iteration << " finished." << std::endl;

    current_time += delta_t;
  }

  std::cout << "output written. Terminating..." << std::endl;
  return 0;
}

void calculateF() {

  for (auto &p : particles) {
    p.old_f = p.f;
    p.f = {0, 0, 0};
  }

  for (long unsigned int i = 0; i < particles.size(); i++) {
    for (long unsigned int j = 0; j < i; j++) {
      auto &p1 = particles[i];
      auto &p2 = particles[j];
      auto x_diff = p2.x - p1.x;

      auto norm = ArrayUtils::L2Norm(x_diff);
      auto f = (p1.m * p2.m) / pow(norm, 3) * x_diff;
      p1.f = p1.f + f;
      p2.f = p2.f - f;
    }
  }
}

void calculateX() {
  for (auto &p : particles) {
    p.x = p.x + delta_t * p.v + pow(delta_t, 2) * (1 / (2 * p.m)) * p.old_f;
  }
}

void calculateV() {
  for (auto &p : particles) {
    p.v = p.v + delta_t * (1 / (2 * p.m)) * (p.old_f + p.f);
  }
}

void plotParticles(int iteration) {

  std::string out_name("MD_vtk");

  outputWriter::XYZWriter writer;
  writer.plotParticles(particles, out_name, iteration);

  outputWriter::VTKWriter vtk_writer;
  assert(particles.size() <= INT_MAX);
  vtk_writer.initializeOutput(static_cast<int>(particles.size()));

  for (auto &p: particles) {
    vtk_writer.plotParticle(p);
  }
  vtk_writer.writeFile(out_name, iteration);
}
