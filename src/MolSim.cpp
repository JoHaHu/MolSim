
#include "FileReader.h"
#include "outputWriter/XYZWriter.h"
#include "utils/ArrayUtils.h"

#include <iostream>
#include <list>

/**** forward declaration of the calculation functions ****/

/**
 * calculate the force for all particles
 */
void calculateF();

/**
 * calculate the position for all particles
 */
void calculateX();

/**
 * calculate the position for all particles
 */
void calculateV();

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration);

constexpr double start_time = 0;
constexpr double end_time = 1000;
constexpr double delta_t = 0.014;

// TODO: what data structure to pick?
std::vector<Particle> particles;

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
    p.setOldF(p.getF());
    p.setF({0., 0., 0.});
  }

  for (long unsigned int i = 1; i <= particles.size(); i++) {
    for (long unsigned int j = 0; j < i; j++) {
      auto &p1 = particles[i];
      auto &p2 = particles[j];
      auto x_diff = p1.getX() - p2.getX();

      auto f = ArrayUtils::elementWiseScalarOp(p1.getM() * p2.getM() / pow(ArrayUtils::L2Norm(x_diff), 3), x_diff, std::multiplies());
      p1.setF(ArrayUtils::elementWisePairOp(p1.getF(), f, std::plus()));
      p2.setF(ArrayUtils::elementWisePairOp(p2.getF(), f, std::plus()));
    }
  }
}

void calculateX() {
  for (auto &p : particles) {
  }
}

void calculateV() {
  for (auto &p : particles) {
    // @TODO: insert calculation of veclocity updates here!
  }
}

void plotParticles(int iteration) {

  std::string out_name("MD_vtk");

  outputWriter::XYZWriter writer;
  writer.plotParticles(particles, out_name, iteration);
}
