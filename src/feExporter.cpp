#include "feExporter.h"
#include "feVertex.h"

#include <iostream>
#include <fstream>
#include <map>

extern bool FE_VERBOSE;

static int VTK_VERTEX = 1;
static int VTK_LINE = 3;
static int VTK_TRIANGLE = 5;
static int VTK_QUAD = 9;
static int VTK_TETRA = 10;
static int VTK_HEXAHEDRON = 12;
static int VTK_PYRAMID = 14;
static int VTK_QUADRATIC_EDGE = 21;
static int VTK_QUADRATIC_TRIANGLE = 22;
static int VTK_QUADRATIC_QUAD = 23;
static int VTK_QUADRATIC_TETRA = 24;
static int VTK_QUADRATIC_HEXAHEDRON = 25;

// A etoffer au fur et a mesure
std::map<std::string, int> cncToVTKmap = {{"Point0D", VTK_VERTEX},
                                          {"LineP1", VTK_LINE},
                                          {"TriP1", VTK_TRIANGLE},
                                          {"LineP2", VTK_QUADRATIC_EDGE},
                                          {"TriP2", VTK_QUADRATIC_TRIANGLE}};

feStatus createVisualizationExporter(feExporter *&exporter, visualizationFormat format,
                                     feMetaNumber *metaNumber, feSolution *solution, feMesh *mesh,
                                     const std::vector<feSpace *> &feSpaces)
{
  switch(format) {
    case VTK:
      exporter = new feExporterVTK(mesh, solution, metaNumber, feSpaces);
      break;
    default:
      return feErrorMsg(FE_STATUS_ERROR, "Unsupported visualization format.");
  }
  return FE_STATUS_OK;
}

void feExporterVTK::writeHeader(std::ostream &output)
{
  output << "# vtk DataFile Version 4.2\n";
  output << "test output\n";
  output << "ASCII\n";
}

void feExporterVTK::writeNodes(std::ostream &output, feCncGeo *cnc)
{
  // int nNod = cnc->getNbNodes();
  int nNod = _mesh->getNbNodes();
  output << "DATASET UNSTRUCTURED_GRID\n";
  output << "POINTS " << nNod << " double\n";
  Vertex *v;
  for(int i = 0; i < nNod; ++i) {
    // v = _mesh->getVertex(cnc->getNodeConnectivity(i));
    v = _mesh->getVertex(i);
    output << v->x() << " " << v->y() << " " << v->z() << std::endl;
  }
}

void feExporterVTK::writeElementsConnectivity(std::ostream &output, feCncGeo *cnc)
{
  int nElm = cnc->getNbElm();
  int nNodePerElem = cnc->getNbNodePerElem();
  output << "CELLS " << nElm << " " << nElm * (nNodePerElem + 1) << std::endl;
  int nodemax = -1;
  for(int iElm = 0; iElm < nElm; ++iElm) {
    output << nNodePerElem << " ";
    for(int j = 0; j < nNodePerElem; ++j) {
      int node = cnc->getNodeConnectivity(iElm, j);
      nodemax = fmax(nodemax, node);
      output << node << " ";
    }
    output << std::endl;
  }
  output << "CELL_TYPES " << nElm << std::endl;
  int vtkElem = cncToVTKmap[cnc->getForme()];
  for(int iElm = 0; iElm < nElm; ++iElm) {
    output << vtkElem << std::endl;
  }
}

void feExporterVTK::writeField(std::ostream &output, feCncGeo *cnc, feSpace *intSpace,
                               std::string fieldID, bool LoopOverCnc)
{
  std::vector<double> &sol = _sol->getSolutionReference();
  // int nNod = cnc->getNbNodes();
  int nNod = LoopOverCnc ? cnc->getNbNodes() : _mesh->getNbNodes();
  feNumber *n = _metaNumber->getNumbering(fieldID);
  // int nNod = n->getNbNodes();
  output << "SCALARS " << fieldID << " double 1" << std::endl;
  output << "LOOKUP_TABLE default" << std::endl;

  // Write field(s) to a text file
  std::string fileName = "solution" + fieldID + ".txt";
  FILE *f = fopen(fileName.c_str(), "w");
  std::vector<double> V(_mesh->getNbNodes(), 0.);
  for(int iNode = 0; iNode < nNod; ++iNode) {
    int iDOF = n->getDOFNumberAtVertex(iNode);
    if(iDOF >= 0) {
      // There is a degree of freedom at this mesh vertex
      output << sol[iDOF] << std::endl;
    } else {
      /* No dof associated to the mesh vertex.
      Interpolate solution at vertex.
      First locate vertex in elements. */
      Vertex *v = _mesh->getVertex(iNode);
      std::vector<double> x = {v->x(), v->y(), v->z()};
      std::vector<double> r(3, 0.0);
      int elm;
      _mesh->locateVertex(x, elm, r);
      intSpace->initializeAddressingVector(n, elm);
      intSpace->initializeSolution(sol);

      double val;
      if(intSpace->useGlobalFunctions()) {
        std::vector<double> geoCoord = _mesh->getCoord(intSpace->getCncGeoTag(), elm);
        std::vector<double> xc(3, 0.0);
        double rc[3] = {1. / 3., 1. / 3., 1. / 3.};
        intSpace->getCncGeo()->getFeSpace()->interpolateVectorField(geoCoord, rc, xc);

        x[0] -= xc[0];
        x[1] -= xc[1];
        x[2] -= xc[2];

        val = intSpace->interpolateSolution(elm, x);
      } else {
        val = intSpace->interpolateSolution(r.data());
      }

      output << val << std::endl;
    }
    fprintf(f, "%+-16.16e\n", sol[iDOF]);
  }
  fclose(f);
}

feStatus feExporterVTK::writeStep(std::string fileName)
{
  std::filebuf fb;
  if(fb.open(fileName, std::ios::out)) {
    std::ostream output(&fb);
    writeHeader(output);

    // For now, we only print spaces of maximal dimension (non boundary)
    std::vector<feSpace *> spacesToExport;
    std::set<std::string> cncToExport;
    for(feSpace *fS : _spaces) {
      if(fS->getDim() == _mesh->getDim()) {
        spacesToExport.push_back(fS);
        cncToExport.insert(fS->getCncGeoID());
      }
    }

    // For now we only print one domain connectivity, assuming all fields are defined on the same
    // connectivity. To fix this, we have to join the connectivities : check which nodes are shared,
    // and use a global numbering.
    if(cncToExport.size() > 1) {
      feWarning("Multiple domains visualization is not implemented. Only the first domain will be "
                "exported.");
    }

    // Grab the connectivity from any matching fespace
    feCncGeo *cnc;
    for(feSpace *fS : spacesToExport) {
      if(fS->getCncGeoID() == *cncToExport.begin()) {
        // if(fS->getCncGeoID() == "first") {
        cnc = fS->getCncGeo();
        break;
      }
    }

    feInfoCond(FE_VERBOSE > 0, "Exporting geometric connectivity \"%s\" to file \"%s\"",
               cnc->getID().c_str(), fileName.c_str());

    // Write nodes and elements
    writeNodes(output, cnc);
    writeElementsConnectivity(output, cnc);
    // output << "POINT_DATA " << cnc->getNbNodes() << std::endl;
    output << "POINT_DATA " << _mesh->getNbNodes() << std::endl;

    // Write the field associated to each fespace in spacesToExport
    for(feSpace *fS : spacesToExport) {
      if(cnc->getForme() == "TriP1" || cnc->getForme() == "TriP2") {
        writeField(output, cnc, fS, fS->getFieldID(), false);
      } else {
        return feErrorMsg(FE_STATUS_WRITE_ERROR,
                          "So far only P1 and P2 triangles can be exported to a VTK file...");
      }
    }

    fb.close();
  } else {
    return feErrorMsg(FE_STATUS_FILE_CANNOT_BE_OPENED, "Could not open output file.");
  }

  return FE_STATUS_OK;
}
