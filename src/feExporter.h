#ifndef _FEEXPORTER_
#define _FEEXPORTER_

#include <iostream>

#include "feMesh.h"
#include "feSolution.h"
#include "feNumber.h"
#include "feEigenProblem.h"

/* Supported visualization formats */
typedef enum { VTK } visualizationFormat;

class feExporter
{
public:
  feMesh *_mesh;
  feSolution *_sol;
  feMetaNumber *_metaNumber;
  const std::vector<feSpace *> &_spaces;

  feEigenProblem *_eigenProblem;

  bool _addP2Nodes;
  int _writtenNodes;

  std::map<int, int> edgeToMid;

public:
  feExporter(feMesh *mesh, feSolution *sol, feMetaNumber *metaNumber,
             const std::vector<feSpace *> &feSpaces)
    : _mesh(mesh), _sol(sol), _metaNumber(metaNumber), _spaces(feSpaces){};
  virtual ~feExporter() {}

  virtual feStatus writeStep(std::string fileName) = 0;
  virtual feStatus writeEigenvectors(std::string fileName) = 0;
};

class feExporterVTK : public feExporter
{
protected:
public:
  feExporterVTK(feMesh *mesh, feSolution *sol, feMetaNumber *metaNumber,
                const std::vector<feSpace *> &feSpaces)
    : feExporter(mesh, sol, metaNumber, feSpaces){};
  feExporterVTK(feMesh *mesh, feSolution *sol, feEigenProblem *eigenProblem,
                feMetaNumber *metaNumber, const std::vector<feSpace *> &feSpaces);
  virtual ~feExporterVTK() {}

  virtual feStatus writeStep(std::string fileName);
  virtual feStatus writeEigenvectors(std::string fileName);

private:
  void writeHeader(std::ostream &output);
  void writeNodes(std::ostream &output);
  void writeElementsConnectivity(std::ostream &output);
  void writeField(std::ostream &output, feCncGeo *cnc, feSpace *intSpace, std::string fieldID,
                  bool loopOverCnc = false);
  void writeEigenvector(std::ostream &output, feCncGeo *cnc, feSpace *intSpace, std::string fieldID,
                        int eigenPairIndex, size_t nEigenPairs, eigenPair &ep,
                        bool loopOverCnc = false);
};

typedef struct feExportData {
  feExporter *exporter;
  int exportEveryNSteps;
  std::string fileNameRoot;
} feExportData;

feStatus createVisualizationExporter(feExporter *&exporter, visualizationFormat format,
                                     feMetaNumber *metaNumber, feSolution *solution, feMesh *mesh,
                                     const std::vector<feSpace *> &feSpaces);

feStatus createVisualizationExporter(feExporter *&exporter, visualizationFormat format,
                                     feMetaNumber *metaNumber, feSolution *solution,
                                     feEigenProblem *eigenProblem, feMesh *mesh,
                                     const std::vector<feSpace *> &feSpaces);

#endif
