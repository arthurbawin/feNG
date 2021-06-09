#ifndef _FEEXPORTER_
#define _FEEXPORTER_

#include <iostream>

#include "feMesh.h"
#include "feSolution.h"
#include "feNumber.h"

class feExporter{
protected:
  feMesh *_mesh;
  feSolution *_sol;
  feMetaNumber *_metaNumber;
  const std::vector<feSpace*>& _space;
public:
  feExporter(feMesh *mesh, feSolution *sol, feMetaNumber *metaNumber, const std::vector<feSpace*> &space) 
    : _mesh(mesh), _sol(sol), _metaNumber(metaNumber), _space(space) {
  };
  virtual ~feExporter() {}

};

class feExporterVTK : public feExporter{
protected:
public:
  feExporterVTK(std::string vtkFile, feMesh *mesh, feSolution *sol, feMetaNumber *metaNumber,
    const std::vector<feSpace*> &space);
  virtual ~feExporterVTK() {}

  void writeHeader(std::ostream& output);
  void writeNodes(std::ostream& output, feCncGeo *cnc);
  void writeElementsConnectivity(std::ostream& output, feCncGeo *cnc);
  void writeField(std::ostream& output, feCncGeo *cnc, std::string fieldID);
};

#endif