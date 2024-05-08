/*
 * VTKWriter.cpp
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#include "VTKWriter.h"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "utils/LoggerManager.h"

namespace outputWriter {

VTKWriter::VTKWriter() = default;

VTKWriter::~VTKWriter() = default;

void VTKWriter::initializeOutput(int numParticles) {

  vtkFile = new VTKFile_t("UnstructuredGrid");

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

  CellData cellData; // we don't have cell data => leave it empty

  // 3 coordinates
  Points points;
  DataArray_t pointCoordinates(type::Float32, "points", 3);
  points.DataArray().push_back(pointCoordinates);

  Cells cells; // we don't have cells, => leave it empty
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
  delete vtkFile;
}

void VTKWriter::plotParticle(Particle &p) {
  auto logger = LoggerManager::getLogger();
  if (vtkFile->UnstructuredGrid().present()) {
    logger->trace("UnstructuredGrid is present");
  } else {
    logger->warn("ERROR: No UnstructuredGrid present");
  }

  PointData::DataArray_sequence &pointDataSequence =
      vtkFile->UnstructuredGrid()->Piece().PointData().DataArray();
  PointData::DataArray_iterator dataIterator = pointDataSequence.begin();

  dataIterator->push_back(p.m);
  // cout << "Appended mass data in: " << dataIterator->Name();

  dataIterator++;
  dataIterator->push_back(p.v[0]);
  dataIterator->push_back(p.v[1]);
  dataIterator->push_back(p.v[2]);
  // cout << "Appended velocity data in: " << dataIterator->Name();

  dataIterator++;
  dataIterator->push_back(p.old_f[0]);
  dataIterator->push_back(p.old_f[1]);
  dataIterator->push_back(p.old_f[2]);
  // cout << "Appended force data in: " << dataIterator->Name();

  dataIterator++;
  dataIterator->push_back(p.type);

  Points::DataArray_sequence &pointsSequence =
      vtkFile->UnstructuredGrid()->Piece().Points().DataArray();
  Points::DataArray_iterator pointsIterator = pointsSequence.begin();
  pointsIterator->push_back(p.x[0]);
  pointsIterator->push_back(p.x[1]);
  pointsIterator->push_back(p.x[2]);
}

} // namespace outputWriter
