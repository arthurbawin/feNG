#ifndef _FEEXPORTER_
#define _FEEXPORTER_

#include <iostream>

#include "feMesh.h"
#include "feSolution.h"
#include "feNumber.h"
#include "feEigenProblem.h"

/* Supported visualization formats */
typedef enum { VTK } visualizationFormat;

// A hardcoded size for now to avoid mallocs
#define MAX_EXPORTED_FIELDS 8

// A visualization node containing its coordinates
// and its data.
typedef struct vtkNodeStruct {
  double pos[3];
  int globalTag;
  double data[MAX_EXPORTED_FIELDS][3] = 
    {
      {-99., -99., -99.},
      {-99., -99., -99.},
      {-99., -99., -99.},
      {-99., -99., -99.},
      {-99., -99., -99.},
      {-99., -99., -99.},
      {-99., -99., -99.},
      {-99., -99., -99.}
    };
} vtkNode;

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

  int _nNodePerElem;
  int _nEdgePerElem;

  std::map<int, int> _edgeGlobalTag;

  std::set<feCncGeo*> _exportedConnectivities;

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
  std::vector<vtkNode> _vtkNodes;
  bool _recreateVTKNodes = true;

public:
  feExporterVTK(feMesh *mesh, feSolution *sol, feMetaNumber *metaNumber,
                const std::vector<feSpace *> &feSpaces)
    : feExporter(mesh, sol, metaNumber, feSpaces){};
  feExporterVTK(feMesh *mesh, feSolution *sol, feEigenProblem *eigenProblem,
                feMetaNumber *metaNumber, const std::vector<feSpace *> &feSpaces);
  ~feExporterVTK() {}

  feStatus writeStep(std::string fileName);
  feStatus writeEigenvectors(std::string fileName);

private:
  void writeHeader(std::ostream &output);
  
  feStatus createVTKNodes(std::vector<feSpace*> &spacesToExport, std::unordered_map<std::string, int> &fieldTags);
  void writeMesh(std::ostream &output);
  void writeVTKNodes(std::ostream &output, std::unordered_map<std::string, std::pair<int, int>> &fields);

  void writeNodes(std::ostream &output);
  void writeElementsConnectivity(std::ostream &output);
  void writeField(std::ostream &output, feCncGeo *cnc, feSpace *intSpace, std::string fieldID,
                  bool loopOverCnc = false);
  void writeEigenvector(std::ostream &output, feCncGeo *cnc, feSpace *intSpace, std::string fieldID,
                        int eigenPairIndex, size_t nEigenPairs, eigenPair &ep,
                        bool loopOverCnc = false);
};

typedef struct feExportData
{
  feExporter *exporter = nullptr;
  int exportEveryNSteps = 1;
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
