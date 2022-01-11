#ifndef _FEEXPORTER_
#define _FEEXPORTER_

#include <iostream>

#include "feMesh.h"
#include "feSolution.h"
#include "feNumber.h"

/* Supported visualization formats */
typedef enum {
  VTK
} visualizationFormat;

class feExporter {
protected:
  feMesh *_mesh;
  feSolution *_sol;
  feMetaNumber *_metaNumber;
  const std::vector<feSpace *> &_spaces;

public:
  feExporter(feMesh *mesh, feSolution *sol, feMetaNumber *metaNumber,
             const std::vector<feSpace *> &feSpaces)
    : _mesh(mesh), _sol(sol), _metaNumber(metaNumber), _spaces(feSpaces){};
  virtual ~feExporter() {}

  virtual feStatus writeStep(std::string fileName) = 0;
};

class feExporterVTK : public feExporter {
protected:
public:
  feExporterVTK(feMesh *mesh, feSolution *sol, feMetaNumber *metaNumber, const std::vector<feSpace *> &feSpaces)
    : feExporter(mesh, sol, metaNumber, feSpaces){ };
  virtual ~feExporterVTK() {}

  virtual feStatus writeStep(std::string fileName);

private:
  void writeHeader(std::ostream &output);
  void writeNodes(std::ostream &output, feCncGeo *cnc);
  void writeElementsConnectivity(std::ostream &output, feCncGeo *cnc);
  void writeField(std::ostream &output, feCncGeo *cnc, feSpace *intSpace, std::string fieldID);
};

typedef struct feExportData {
  feExporter *exporter;
  int exportEveryNSteps;
  std::string fileNameRoot;
} feExportData;

feStatus createVisualizationExporter(feExporter *&exporter,
                                     visualizationFormat format,
                                     feMetaNumber *metaNumber,
                                     feSolution *solution,
                                     feMesh *mesh,
                                     const std::vector<feSpace*> &feSpaces);

#endif