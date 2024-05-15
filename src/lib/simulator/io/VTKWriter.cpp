#include "lib/simulator/io/VTKWriter.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "lib/utils/LoggerManager.h"
#include "spdlog/spdlog.h"

namespace simulator::io {

void VTKWriter::initializeOutput(int numParticles) {
  vtkFile = std::make_unique<VTKFile_t>(VTKFile_t("UnstructuredGrid"));

  // per point, we add type, position, velocity and force
  PointData pointData;
  DataArray_t mass(type::Float32, "mass", 1);
  DataArray_t velocity(type::Float32, "velocity", 3);
  DataArray_t forces(type::Float32, "force", 3);
  DataArray_t type(type::Int32, "type", 1);
  pointData.DataArray().push_back(mass);
  pointData.DataArray().push_back(velocity);
  pointData.DataArray().push_back(forces);
  pointData.DataArray().push_back(type);

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

void VTKWriter::writeFile(const std::string &filename, int iteration) {
  std::stringstream strstr;
  strstr << filename << "_" << std::setfill('0') << std::setw(4) << iteration << ".vtu";

  std::ofstream file(strstr.str().c_str());
  VTKFile(file, *vtkFile);
  // Early release of memory resources
  vtkFile.reset();
}

void VTKWriter::plotParticle(Particle &p) {
  if (vtkFile->UnstructuredGrid().present()) {
    spdlog::trace("UnstructuredGrid is present");
  } else {
    spdlog::warn("ERROR: No UnstructuredGrid present");
  }

  PointData::DataArray_sequence &pointDataSequence =
      vtkFile->UnstructuredGrid()->Piece().PointData().DataArray();
  PointData::DataArray_iterator dataIterator = pointDataSequence.begin();

  dataIterator->push_back(p.mass);
  // cout << "Appended mass data in: " << dataIterator->Name();

  dataIterator++;
  dataIterator->push_back(p.velocity[0]);
  dataIterator->push_back(p.velocity[1]);
  dataIterator->push_back(p.velocity[2]);
  // cout << "Appended velocity data in: " << dataIterator->Name();

  dataIterator++;
  dataIterator->push_back(p.old_force[0]);
  dataIterator->push_back(p.old_force[1]);
  dataIterator->push_back(p.old_force[2]);
  // cout << "Appended force data in: " << dataIterator->Name();

  dataIterator++;
  dataIterator->push_back(p.type);

  Points::DataArray_sequence &pointsSequence =
      vtkFile->UnstructuredGrid()->Piece().Points().DataArray();
  Points::DataArray_iterator pointsIterator = pointsSequence.begin();
  pointsIterator->push_back(p.position[0]);
  pointsIterator->push_back(p.position[1]);
  pointsIterator->push_back(p.position[2]);
}

}// namespace simulator::io
