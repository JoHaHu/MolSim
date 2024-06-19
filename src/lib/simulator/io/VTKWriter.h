#pragma once

#include "Particle.h"
#include "simulator/io/VTKUnstructured.hxx"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <string>

#include "spdlog/spdlog.h"
namespace simulator::io {

/**
 * This class implements the functionality to generate vtk output from
 * particles.
 */
template<const size_t DIMENSIONS>
class VTKWriter {

 public:
  VTKWriter() = default;

  virtual ~VTKWriter() = default;
  /**
   * set up internal data structures and prepare to plot a particle.
   */
  void initializeOutput(int numParticles) {
    vtkFile = std::make_unique<VTKFile_t>(VTKFile_t("UnstructuredGrid"));

    // per point, we add type, position, velocity and force
    PointData pointData;
    DataArray_t mass(type::Float32, "mass", 1);
    DataArray_t velocity(type::Float32, "velocity", 3);
    DataArray_t forces(type::Float32, "force", 3);
    DataArray_t old_forces(type::Float32, "old_force", 3);
    DataArray_t type(type::Int32, "type", 1);
    pointData.DataArray().push_back(mass);
    pointData.DataArray().push_back(velocity);
    pointData.DataArray().push_back(forces);
    pointData.DataArray().push_back(old_forces);
    pointData.DataArray().push_back(type);
#if DEBUG
    DataArray_t id(type::Int32, "id", 1);
    pointData.DataArray().push_back(id);
#endif
    CellData cellData;// we don't have cell data => leave it empty

    // 3 coordinates
    Points points;
    DataArray_t pointCoordinates(type::Float32, "points", 3);
    points.DataArray().push_back(pointCoordinates);

    Cells cells;// we don't have cells, => leave it empty
    // for some reasons, we have to add a dummy entry for paraview
    DataArray_t cells_data(type::Float32, "types", 0);
    cells.DataArray().push_back(cells_data);

    PieceUnstructuredGrid_t piece(pointData, cellData, points, cells,
                                  numParticles, 0);
    UnstructuredGrid_t unstructuredGrid(piece);
    vtkFile->UnstructuredGrid(unstructuredGrid);
  }

  /**
   * plot type, mass, position, velocity and force of a particle.
   *
   * @note: initializeOutput() must have been called before.
   */
  void plotParticle(Particles<DIMENSIONS> &p, size_t index) {
    if (vtkFile->UnstructuredGrid().present()) {
      SPDLOG_TRACE("UnstructuredGrid is present");
    } else {
      SPDLOG_WARN("ERROR: No UnstructuredGrid present");
    }

    PointData::DataArray_sequence &pointDataSequence =
        vtkFile->UnstructuredGrid()->Piece().PointData().DataArray();
    PointData::DataArray_iterator dataIterator = pointDataSequence.begin();

    dataIterator->push_back(p.mass[index]);
    // cout << "Appended mass data in: " << dataIterator->Name();

    dataIterator++;
    for (int i = 0; i < DIMENSIONS; ++i) {
      dataIterator->push_back(p.velocities[i][index]);
    }
    for (int i = DIMENSIONS; i < 3; ++i) {
      dataIterator->push_back(0);
    }

    // cout << "Appended velocity data in: " << dataIterator->Name();

    dataIterator++;
    for (int i = 0; i < DIMENSIONS; ++i) {
      dataIterator->push_back(p.forces[i][index]);
    }
    for (int i = DIMENSIONS; i < 3; ++i) {
      dataIterator->push_back(0);
    }

    dataIterator++;
    for (int i = 0; i < DIMENSIONS; ++i) {
      dataIterator->push_back(p.old_forces[i][index]);
    }
    for (int i = DIMENSIONS; i < 3; ++i) {
      dataIterator->push_back(0);
    }
    // cout << "Appended force data in: " << dataIterator->Name();

    dataIterator++;
    dataIterator->push_back(p.type[index]);

#if DEBUG
    dataIterator++;
    dataIterator->push_back(p.ids[index]);
#endif

    Points::DataArray_sequence &pointsSequence =
        vtkFile->UnstructuredGrid()->Piece().Points().DataArray();
    Points::DataArray_iterator pointsIterator = pointsSequence.begin();
    for (int i = 0; i < DIMENSIONS; ++i) {
      pointsIterator->push_back(p.positions[i][index]);
    }
    for (int i = DIMENSIONS; i < 3; ++i) {
      pointsIterator->push_back(0);
    }
  }

  /**
   * writes the final output file.
   *
   * @param filename the base name of the file to be written.
   * @param iteration the number of the current iteration,
   *        which is used to generate an unique filename
   */
  void writeFile(const std::string &filename, int iteration) {
    std::stringstream strstr;
    strstr << filename << "_" << std::setfill('0') << std::setw(4) << iteration << ".vtu";
    std::ofstream file(strstr.str().c_str());
    VTKFile(file, *vtkFile);
    // Early release of memory resources
    vtkFile.reset();
  }

 private:
  std::unique_ptr<VTKFile_t> vtkFile = std::make_unique<VTKFile_t>(VTKFile_t("UnstructuredGrid"));
};

}// namespace simulator::io
