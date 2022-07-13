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
  int nNod = _mesh->getNbNodes();
  output << "DATASET UNSTRUCTURED_GRID\n";

  if(_addP2Nodes) {
    output << "POINTS " << nNod + _mesh->getNbEdges() << " double\n";
  } else {
    output << "POINTS " << nNod << " double\n";
  }

  Vertex *v;
  for(int i = 0; i < nNod; ++i) {
    // v = _mesh->getVertex(cnc->getNodeConnectivity(i));
    v = _mesh->getVertex(i);
    output << v->x() << " " << v->y() << " " << v->z() << std::endl;
  }

  int cnt = nNod;
  if(_addP2Nodes) {
    /* Additional nodes are added in the order of the edges,
    like in feExporterVTK::writeFields. Both must be modified together. */
    for(auto e : _mesh->_edges) {
      Vertex *v0 = e.getVertex(0);
      Vertex *v1 = e.getVertex(1);
      output << (v0->x() + v1->x()) / 2. << " " << (v0->y() + v1->y()) / 2. << " "
             << (v0->z() + v1->z()) / 2. << std::endl;
      edgeToMid[e.getTag()] = cnt++;
    }
  }
  _writtenNodes = cnt;
}

void feExporterVTK::writeElementsConnectivity(std::ostream &output, feCncGeo *cnc)
{
  int nElm = cnc->getNbElm();
  int nNodePerElem = _addP2Nodes ? 6 : cnc->getNbNodePerElem();
  output << "CELLS " << nElm << " " << nElm * (nNodePerElem + 1) << std::endl;
  for(int iElm = 0; iElm < nElm; ++iElm) {
    output << nNodePerElem << " ";

    if(_addP2Nodes) {
      // Mid-points were added and do not exist in the mesh
      // Vertices
      for(int j = 0; j < 3; ++j) {
        int node = cnc->getNodeConnectivity(iElm, j);
        output << node << " ";
      }
      // Edges
      for(int j = 0; j < 3; ++j) {
        int iEdge = fabs(cnc->getEdgeConnectivity(iElm, j));
        int node = edgeToMid[iEdge];
        output << node << " ";
      }
    } else {
      // Regular case : all nodes exist in the mesh
      for(int j = 0; j < nNodePerElem; ++j) {
        int node = cnc->getNodeConnectivity(iElm, j);
        output << node << " ";
      }
    }

    output << std::endl;
  }
  output << "CELL_TYPES " << nElm << std::endl;
  int vtkElem = cncToVTKmap[cnc->getForme()];
  for(int iElm = 0; iElm < nElm; ++iElm) {
    if(_addP2Nodes)
      output << VTK_QUADRATIC_TRIANGLE << std::endl;
    else
      output << vtkElem << std::endl;
  }
}

void feExporterVTK::writeField(std::ostream &output, feCncGeo *cnc, feSpace *intSpace,
                               std::string fieldID, bool loopOverCnc)
{
  std::vector<double> &sol = _sol->getSolutionReference();
  // int nVertices = cnc->getNbNodes();
  int nVertices = loopOverCnc ? cnc->getNbNodes() : _mesh->getNbNodes();
  feNumber *n = _metaNumber->getNumbering(fieldID);
  // int nVertices = n->getNbNodes();

  output << "POINT_DATA " << _writtenNodes << std::endl;
  output << "SCALARS " << fieldID << " double 1" << std::endl;
  output << "LOOKUP_TABLE default" << std::endl;

  // Write field(s) to a text file
  std::string fileName = "solution" + fieldID + ".txt";
  FILE *f = fopen(fileName.c_str(), "w");

  for(int iVertex = 0; iVertex < nVertices; ++iVertex) {
    int iDOF = n->getDOFNumberAtVertex(iVertex);

    if(iDOF >= 0) {
      // There is a degree of freedom at this mesh vertex
      output << sol[iDOF] << std::endl;
    } else {
      /* No dof associated to the mesh vertex.
      Interpolate solution at vertex. */
      Vertex *v = _mesh->getVertex(iVertex);
      std::vector<double> x = {v->x(), v->y(), v->z()};
      std::vector<double> r(3, 0.0);
      int elm;
      _mesh->locateVertex(x, elm, r);
      intSpace->initializeAddressingVector(n, elm);
      intSpace->initializeSolution(sol);

      double val;
      if(intSpace->useGlobalFunctions()) {
        val = intSpace->interpolateSolution(elm, x);
      } else {
        val = intSpace->interpolateSolution(r.data());
      }

      output << val << std::endl;
    }
    fprintf(f, "%+-16.16e\n", sol[iDOF]);
  }
  fclose(f);

  /* Write the additional P2 nodes. Interpolation is not required if the
  field is quadratic, but it's easier to just interpolate for all fields. */
  if(_addP2Nodes) {
    /* Additional nodes are added in the order of the edges,
    like in feExporterVTK::writeNodes. Both must be modified together. */
    for(auto e : _mesh->_edges) {
      Vertex *v0 = e.getVertex(0);
      Vertex *v1 = e.getVertex(1);
      std::vector<double> x = {(v0->x() + v1->x()) / 2., (v0->y() + v1->y()) / 2.,
                               (v0->z() + v1->z()) / 2.};
      std::vector<double> r(3, 0.0);
      int elm;
      _mesh->locateVertex(x, elm, r);
      intSpace->initializeAddressingVector(n, elm);
      intSpace->initializeSolution(sol);
      double val;
      if(intSpace->useGlobalFunctions()) {
        val = intSpace->interpolateSolution(elm, x);
      } else {
        val = intSpace->interpolateSolution(r.data());
      }
      output << val << std::endl;
    }
  }
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
        cnc = fS->getCncGeo();
        break;
      }
    }

    feInfoCond(FE_VERBOSE > 0, "Exporting geometric connectivity \"%s\" to file \"%s\"",
               cnc->getID().c_str(), fileName.c_str());

    /* Although VTK_HIGHER_ORDER_TRIANGLE and VTK_LAGRANGE_TRIANGLE
    exist in the VTK documentation, I can't find tutorials on how to use them...
    So for now the highest order triangle is the VTK_QUADRATIC_TRIANGLE.
    Quadratic triangles are exported if :
      - the geometry is P2 (curved elements)
      - any field has a P2 or greater interpolant
    In the second case, the mesh is not made of 6-node triangles, so we need to create them
    and evaluate the fields at those mid-edge nodes. */
    int geometryPolynomialDegree;
    if(cnc->getForme() == "TriP1") {
      geometryPolynomialDegree = 1;
    } else if(cnc->getForme() == "TriP2") {
      geometryPolynomialDegree = 2;
    } else {
      return feErrorMsg(FE_STATUS_WRITE_ERROR,
                        "Only P1 and P2 triangles can be exported to a VTK file...");
    }

    int highestFieldPolynomialDegree = 0;
    for(feSpace *fS : spacesToExport)
      highestFieldPolynomialDegree = fmax(highestFieldPolynomialDegree, fS->getPolynomialDegree());

    _addP2Nodes = (geometryPolynomialDegree == 1) && (highestFieldPolynomialDegree >= 2);

    // Write nodes and elements
    writeNodes(output, cnc);
    writeElementsConnectivity(output, cnc);

    // Write the field associated to each fespace in spacesToExport
    for(feSpace *fS : spacesToExport) {
      writeField(output, cnc, fS, fS->getFieldID(), false);
    }

    fb.close();
  } else {
    return feErrorMsg(FE_STATUS_FILE_CANNOT_BE_OPENED, "Could not open output file.");
  }

  return FE_STATUS_OK;
}
