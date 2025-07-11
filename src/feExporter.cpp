#include "feExporter.h"
#include "feVertex.h"
#include "feCncGeo.h"

#include <iostream>
#include <fstream>
#include <map>

#if defined(HAVE_MPI)
  #include "mpi.h"
#endif

#define TOL_ZERO_VTK 1e-20

extern bool FE_VERBOSE;

static int VTK_VERTEX = 1;
static int VTK_LINE = 3;
static int VTK_TRIANGLE = 5;
// static int VTK_QUAD = 9;
static int VTK_TETRA = 10;
// static int VTK_HEXAHEDRON = 12;
// static int VTK_PYRAMID = 14;
static int VTK_QUADRATIC_EDGE = 21;
static int VTK_QUADRATIC_TRIANGLE = 22;
// static int VTK_QUADRATIC_QUAD = 23;
// static int VTK_QUADRATIC_TETRA = 24;
// static int VTK_QUADRATIC_HEXAHEDRON = 25;

std::map<geometricInterpolant, int> cncToVTKmap = {
  {geometricInterpolant::POINTP0, VTK_VERTEX},
  {geometricInterpolant::LINEP1,  VTK_LINE},
  {geometricInterpolant::LINEP2,  VTK_QUADRATIC_EDGE},
  {geometricInterpolant::TRIP1,   VTK_TRIANGLE},
  {geometricInterpolant::TRIP2,   VTK_QUADRATIC_TRIANGLE},
  {geometricInterpolant::TETP1,   VTK_TETRA},
};

static std::map<geometricInterpolant, int> verticesPerElem = {
  {geometricInterpolant::POINTP0, 1},
  {geometricInterpolant::LINEP1,  2},
  {geometricInterpolant::LINEP2,  2},
  {geometricInterpolant::TRIP1,   3},
  {geometricInterpolant::TRIP2,   6},
  {geometricInterpolant::TETP1,   4}
};

static std::map<geometricInterpolant, int> edgePerElem = {
  {geometricInterpolant::POINTP0, 0},
  {geometricInterpolant::LINEP1,  1},
  {geometricInterpolant::LINEP2,  1},
  {geometricInterpolant::TRIP1,   3},
  {geometricInterpolant::TRIP2,   3},
  {geometricInterpolant::TETP1,   6}
};

// Reference coordinates to evaluate the fields
static double referenceCoordinatesTri[3][3] = {
  {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
static double referenceCoordinatesQuadraticTri[6][3] = {
  {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
  {0.5, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.0, 0.5, 0.0}
};

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

feStatus createVisualizationExporter(feExporter *&exporter, visualizationFormat format,
                                     feMetaNumber *metaNumber, feSolution *solution,
                                     feEigenProblem *eigenProblem, feMesh *mesh,
                                     const std::vector<feSpace *> &feSpaces)
{
  switch(format) {
    case VTK:
      exporter = new feExporterVTK(mesh, solution, eigenProblem, metaNumber, feSpaces);
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

void feExporterVTK::writeNodes(std::ostream &output)
{
  int nVertices = _mesh->getNumVertices();
  output << "DATASET UNSTRUCTURED_GRID\n";

  if(_addP2Nodes) {
    output << "POINTS " << nVertices + _mesh->getNumEdges() << " double\n";
  } else {
    output << "POINTS " << nVertices << " double\n";
  }

  Vertex *v;
  for(int i = 0; i < nVertices; ++i) {
    // v = _mesh->getVertex(cnc->getVertexConnectivity(i));
    v = _mesh->getVertex(i);
    output << v->x() << " " << v->y() << " " << v->z() << std::endl;
  }

  int cnt = nVertices;
  if(_addP2Nodes) {
    /* Additional nodes are added in the order of the edges,
    like in feExporterVTK::writeFields. Both must be modified together. */
    for(auto e : _mesh->_edges) {
      Vertex *v0 = e.getVertex(0);
      Vertex *v1 = e.getVertex(1);
      output << (v0->x() + v1->x()) / 2. << " " << (v0->y() + v1->y()) / 2. << " "
             << (v0->z() + v1->z()) / 2. << std::endl;
      _edgeGlobalTag[e.getTag()] = cnt++;
    }
  }
  _writtenNodes = cnt;
}

void feExporterVTK::writeElementsConnectivity(std::ostream &output)
{
  int nElm = 0;
  for(auto *cnc : _mesh->getCncGeo()) {
    if(cnc->getDim() == 2) {
      nElm += cnc->getNumElements();
    }
  }

  int nNodePerElem = _addP2Nodes ? 6 : 3;

  // Write elements connectivities for each cnc wih dim = 2
  output << "CELLS " << nElm << " " << nElm * (nNodePerElem + 1) << std::endl;

  for(auto *cnc : _mesh->getCncGeo()) {
    if(cnc->getDim() == 2) {
      for(int iElm = 0; iElm < cnc->getNumElements(); ++iElm) {
        output << nNodePerElem << " ";

        if(_addP2Nodes) {
          // Mid-points were added and do not exist in the mesh
          // Vertices
          for(int j = 0; j < 3; ++j) {
            int node = cnc->getVertexConnectivity(iElm, j);
            output << node << " ";
          }
          // Edges
          for(int j = 0; j < 3; ++j) {
            int iEdge = fabs(cnc->getEdgeConnectivity(iElm, j));
            int node = _edgeGlobalTag[iEdge];
            output << node << " ";
          }
        } else {
          // Regular case : all nodes exist in the mesh
          for(int j = 0; j < cnc->getNumVerticesPerElem(); ++j) {
            int node = cnc->getVertexConnectivity(iElm, j);
            output << node << " ";
          }
        }

        output << std::endl;
      }
    }
  }

  // Write element types
  output << "CELL_TYPES " << nElm << std::endl;

  for(auto *cnc : _mesh->getCncGeo()) {
    if(cnc->getDim() == 2) {
      int vtkElem = cncToVTKmap[cnc->getInterpolant()];
      for(int iElm = 0; iElm < cnc->getNumElements(); ++iElm) {
        if(_addP2Nodes)
          output << VTK_QUADRATIC_TRIANGLE << std::endl;
        else
          output << vtkElem << std::endl;
      }
    }
  }
}

void feExporterVTK::writeField(std::ostream &output, const feCncGeo *cnc, feSpace *intSpace,
                               std::string fieldID, bool loopOverCnc)
{
  UNUSED(cnc, loopOverCnc);

  std::vector<double> &solVec = _sol->getSolutionReference();
  std::vector<feInt> adr(intSpace->getNumFunctions());
  std::vector<double> sol(intSpace->getNumFunctions());
  // int nVertices = loopOverCnc ? cnc->getNumVertices() : _mesh->getNumVertices();
  int nVertices = _mesh->getNumVertices();
  feNumber *n = _metaNumber->getNumbering(fieldID);
  int nComponents = intSpace->getNumComponents();

  // output << "POINT_DATA " << _writtenNodes << std::endl;
  // output << "SCALARS " << fieldID << " double 1" << std::endl;
  // output << "LOOKUP_TABLE default" << std::endl;

  output << fieldID << " " << nComponents << " " << _writtenNodes << " double" << std::endl;

  int iDOF, elm;
  Vertex *v;
  std::vector<double> x;
  std::vector<double> res(3, 0.);
  double r[3], val;

  // Write field(s) to a text file
  std::string fileName = "solution" + fieldID + ".txt";
  FILE *f = fopen(fileName.c_str(), "w");
  for(int iVertex = 0; iVertex < nVertices; ++iVertex) {
    // Loop over components for vector-valued FE spaces
    for(int iComp = 0; iComp < nComponents; ++iComp) {
      iDOF = n->getDOFNumberAtVertex(iVertex, iComp);

      if(iDOF >= 0) {
        // There is a degree of freedom at this mesh vertex

        // FIXME: This is only valid for Lagrange type elements,
        // where to DOF is the function evaluation. We should interpolate.
        output << solVec[iDOF] << std::endl;

      } else {
        /* No dof associated to the mesh vertex.
        Interpolate solution at vertex. */
        v = _mesh->getVertex(iVertex);
        x = {v->x(), v->y(), v->z()};
        _mesh->locateVertex(x.data(), elm, r);
        intSpace->initializeAddressingVector(elm, adr);

        for(size_t i = 0; i < adr.size(); ++i) {
          sol[i] = solVec[adr[i]];
        }

        if(intSpace->useGlobalFunctions()) {
          val = intSpace->interpolateField(sol, elm, x);
        } else {
          val = intSpace->interpolateField(sol, r);
        }

        output << val << std::endl;
      }
      fprintf(f, "%+-16.16e\n", solVec[iDOF]);
    }
  }
  fclose(f);

  /* Write the additional P2 nodes. Interpolation is not required if the
  field is quadratic, but it's easier to just interpolate for all fields. */
  if(_addP2Nodes) {
    /* Additional nodes are added in the order of the edges,
    like in feExporterVTK::writeNodes. Both must be modified together. */
    Vertex *v0, *v1;
    for(auto e : _mesh->_edges) {
      v0 = e.getVertex(0);
      v1 = e.getVertex(1);
      x = {(v0->x() + v1->x()) / 2., (v0->y() + v1->y()) / 2., (v0->z() + v1->z()) / 2.};
      _mesh->locateVertex(x.data(), elm, r);
      intSpace->initializeAddressingVector(elm, adr);
      for(size_t i = 0; i < adr.size(); ++i) sol[i] = solVec[adr[i]];
      if(intSpace->useGlobalFunctions()) {
        val = intSpace->interpolateField(sol, elm, x);
        output << val << std::endl;
      } else {
        if(nComponents > 1) {
          intSpace->interpolateVectorField(sol, nComponents, r, res);
          for(int i = 0; i < nComponents; ++i) output << res[i] << std::endl;
        } else {
          val = intSpace->interpolateField(sol, r);
          output << val << std::endl;
        }
      }
    }
  }
}

feStatus feExporterVTK::createVTKNodes(std::vector<feSpace *> &spacesToExport,
                                       std::unordered_map<std::string, int> &fieldTags)
{
  // Create the nodes (only once)
  if(_recreateVTKNodes) {

      _recreateVTKNodes = false;

    if(_addP2Nodes)
      _vtkNodes.resize(_mesh->getNumVertices() + _mesh->getNumEdges());
    else
      _vtkNodes.resize(_mesh->getNumVertices());

    // Create nodes with tag and coordinates
    for(auto *cnc : _exportedConnectivities) {

      #if defined(HAVE_OMP)
      #pragma omp parallel for
      #endif
      for(int i = 0; i < cnc->getNumVertices(); ++i) {
        vtkNode node;
        int globalTag = cnc->getUniqueVertexConnectivity(i);
        node.globalTag = globalTag;
        Vertex *v = _mesh->getVertex(globalTag);
        node.pos[0] = v->x();
        node.pos[1] = v->y();
        node.pos[2] = v->z();

        _vtkNodes[globalTag] = node;
      }
    }

    // Add visualization nodes at mid-edges if necessary
    int cnt = _mesh->getNumVertices();
    if(_addP2Nodes) {
      for(auto e : _mesh->_edges) {
        vtkNode node;
        int globalTag = cnt++;
        node.globalTag = globalTag;
        Vertex *v0 = e.getVertex(0);
        Vertex *v1 = e.getVertex(1);
        node.pos[0] = (v0->x() + v1->x()) / 2.;
        node.pos[1] = (v0->y() + v1->y()) / 2.;
        node.pos[2] = (v0->z() + v1->z()) / 2.;

        _edgeGlobalTag[e.getTag()] = globalTag;
        _vtkNodes[globalTag] = node;
      }
    }
  }

  const std::vector<double> &solVec = _sol->getSolution();

  // Set vtkNodes data
  for(auto *space : spacesToExport) {
    std::string field = space->getFieldID();
    int iField = fieldTags[field];
    int nComponents = space->getNumComponents();
    // Write 3 components even for 2D vector fields
    int nComponentsToWrite = (nComponents > 1) ? 3 : 1;
    // int nComponentsToWrite = (nComponents > 1) ? nComponents : 1;
    feNumber *numbering = _metaNumber->getNumbering(field);

    feInfoCond(FE_VERBOSE > VERBOSE_MODERATE, 
      "Exporting field %s on connectivity %s", field.data(), space->getCncGeoID().data());

    const feCncGeo *cnc = space->getCncGeo();

    bool OK = true;
    double failedX[3];

    // Loop over mesh vertices on this connectivity and evaluate field
    #if defined(HAVE_OMP)
    #pragma omp parallel shared(OK, failedX)
    #endif
    {
      int iDOF, elm;
      std::vector<double> x;
      std::vector<double> res(3, 0.);
      double r[3], val;
      std::vector<feInt> adr(space->getNumFunctions());
      std::vector<double> sol(space->getNumFunctions());

      #if defined(HAVE_OMP)
      #pragma omp for
      #endif
      for(int i = 0; i < cnc->getNumVertices(); ++i) {
        int globalTag = cnc->getUniqueVertexConnectivity(i);

        // Loop over components for vector-valued FE spaces
        for(int iComp = 0; iComp < nComponentsToWrite; ++iComp) {

          if(iComp < nComponents) {
            iDOF = numbering->getDOFNumberAtVertex(globalTag, iComp);

            if(iDOF >= 0) {
              // There is a degree of freedom at this mesh vertex

              // FIXME: This is only valid for Lagrange type elements,
              // where to DOF is the function evaluation. We should interpolate.
              val = solVec[iDOF];
              if(fabs(val) < TOL_ZERO_VTK) val = 0.;
              _vtkNodes[globalTag].data[iField][iComp] = val;

            } else {
              /* No dof associated to the mesh vertex,
              e.g. P1 pressure at P2 visualization node.
              Interpolate solution at vertex. */
              vtkNode node = _vtkNodes[globalTag];
              x = {node.pos[0], node.pos[1], node.pos[2]};
              bool wasFound = _mesh->locateVertex(x.data(), elm, r, 1e-6);
              if(!wasFound) {
                feInfo("Could not find point %f - %f - %f when exporting", x[0], x[1], x[2]);
                exit(-1);
              }
              space->initializeAddressingVector(elm, adr);

              for(size_t j = 0; j < adr.size(); ++j) {
                sol[j] = solVec[adr[j]];
              }

              if(space->useGlobalFunctions()) {
                feInfo("FIXME: Choose interpolation function for vector field and global functions.");
                exit(-1);
                val = space->interpolateField(sol, elm, x);
              } else {
                val = space->interpolateVectorFieldComponent(sol, iComp, r);
              }

              if(fabs(val) < TOL_ZERO_VTK) val = 0.;
              _vtkNodes[globalTag].data[iField][iComp] = val;
            }
          } else {
            _vtkNodes[globalTag].data[iField][iComp] = 0.;
          }
        }
      }

      /* Write the additional P2 nodes. Interpolation is not required if the
      field is quadratic, but it's easier to just interpolate for all fields. */
      if(_addP2Nodes) {
        #if defined(HAVE_OMP)
        #pragma omp for
        #endif
        for(int iEdge = 0; iEdge < cnc->getNumEdges(); ++iEdge)
        {
          if(!OK) continue;

          int edgeTag = fabs(cnc->getUniqueEdgeConnectivity(iEdge));
          const Edge *e = _mesh->getEdge(edgeTag);
          int globalTag = _edgeGlobalTag[e->getTag()];
          const Vertex *v0 = e->getVertex(0);
          const Vertex *v1 = e->getVertex(1);
          x = {(v0->x() + v1->x()) / 2., (v0->y() + v1->y()) / 2., (v0->z() + v1->z()) / 2.};

          bool returnLocalElmTag = true;
          std::string targetConnectivity = cnc->getID();
          bool ret =
            _mesh->locateVertex(x.data(), elm, r, 1e-5, returnLocalElmTag, targetConnectivity);

          if(!ret) {
            OK = false;
            #if defined(HAVE_OMP)
            #pragma omp critical
            #endif
            {
              failedX[0] = x[0];
              failedX[1] = x[1];
              failedX[2] = x[2];
            }
            // return feErrorMsg(FE_STATUS_ERROR, "Could not find edge vertex (%+-1.6e, %+-1.6e, %+-1.6e) in the mesh",
            //                   x[0], x[1], x[2]);
          }

          space->initializeAddressingVector(elm, adr);

          for(size_t i = 0; i < adr.size(); ++i) sol[i] = solVec[adr[i]];

          if(space->useGlobalFunctions()) {
            val = space->interpolateField(sol, elm, x);
            // FIXME: vector field with global functions?
            if(fabs(val) < TOL_ZERO_VTK) val = 0.;
            _vtkNodes[globalTag].data[iField][0] = val;
          } else {
            if(nComponents > 1) {
              space->interpolateVectorField(sol, nComponents, r, res);
              for(int i = 0; i < nComponentsToWrite; ++i) {
                _vtkNodes[globalTag].data[iField][i] = res[i];
              }
            } else {
              val = space->interpolateField(sol, r);
              if(fabs(val) < TOL_ZERO_VTK) val = 0.;
              _vtkNodes[globalTag].data[iField][0] = val;
            }
          }
        }
      }
    }

    if(!OK) {
      return feErrorMsg(FE_STATUS_ERROR, "Could not find edge vertex (%+-1.6e, %+-1.6e, %+-1.6e) in the mesh",
                                failedX[0], failedX[1], failedX[2]);
    }
  }

  return FE_STATUS_OK;
}

feStatus feExporterVTK::createDiscontinuousVTKNodes(std::vector<feSpace*> &spacesToExport,
                                                    std::unordered_map<std::string, int> &fieldTags)
{
  if(_recreateVTKNodes) {
    _recreateVTKNodes = false;

    const feCncGeo *cnc = spacesToExport[0]->getCncGeo();
    const std::string &cncName = cnc->getID();

    int numElem = _mesh->getNumInteriorElements();
    int numVerticesPerElem = _mesh->getNumVerticesPerElem(cncName);

    // _nNodePerElem already includes additional P2 nodes if needed
    _vtkNodes.resize(numElem * _nNodePerElem);

    // Loop over mesh elements and create duplicated nodes
    #if defined(HAVE_OMP)
    #pragma omp parallel for
    #endif
    for(int iElm = 0; iElm < cnc->getNumElements(); ++iElm)
    {
      for(int j = 0; j < numVerticesPerElem; ++j)
      {
        vtkNode node;
        int iVertex = _mesh->getVertex(cncName, iElm, j);
        Vertex *v = _mesh->getVertex(iVertex);
        node.pos[0] = v->x();
        node.pos[1] = v->y();
        node.pos[2] = v->z();

        // A unique tag for each duplicated VTK node
        node.globalTag = iElm * _nNodePerElem + j;

        _vtkNodes[node.globalTag] = node;
      }

      if(_addP2Nodes) {
        for(int j = 0; j < numVerticesPerElem; ++j) {
          vtkNode node;
          int iv0 = _mesh->getVertex(cncName, iElm, j);
          int iv1 = _mesh->getVertex(cncName, iElm, (j+1) % 3);
          Vertex *v0 = _mesh->getVertex(iv0);
          Vertex *v1 = _mesh->getVertex(iv1);
          node.pos[0] = (v0->x() + v1->x())/2.;
          node.pos[1] = (v0->y() + v1->y())/2.;
          node.pos[2] = (v0->z() + v1->z())/2.;

          // A unique tag for each duplicated VTK node
          node.globalTag = iElm * _nNodePerElem + 3 + j;

          _vtkNodes[node.globalTag] = node;
        }
      }
    }

    const std::vector<double> &solArray = _sol->getSolution();

    // Set vtkNodes data
    for(auto *space : spacesToExport)
    {
      std::string field = space->getFieldID();
      int iField = fieldTags[field];
      int nComponents = space->getNumComponents();
      // Write 3 components even for 2D vector fields
      int nComponentsToWrite = (nComponents > 1) ? 3 : 1;

      feInfoCond(FE_VERBOSE > VERBOSE_MODERATE, 
        "Exporting field %s on connectivity %s", field.data(), space->getCncGeoID().data());

      // Loop over mesh vertices on this connectivity and evaluate field
      #if defined(HAVE_OMP)
      #pragma omp parallel
      #endif
      {
        std::vector<feInt> adr(space->getNumFunctions());
        std::vector<double> solOnElem(space->getNumFunctions());
        std::vector<double> res(3, 0.);

        ///////////////////////////
        // ElementTransformation T;
        // // const std::vector<int> &edgeconnec = cnc->getEdgeConnectivity();
        // // int sign[4][2] = {{1., -1.},{1., -1.},{1., -1.},{1., -1.}};
        // bool firstForEdge[4] = {false, false, false, false};
        // ///////////////////////////
        // FILE *myFile = fopen("normals.pos", "w");
        // fprintf(myFile, "View \"%s\"{\n", "normals");
        // bool plot = true;
        // ///////////////////////////

        #if defined(HAVE_OMP)
        #pragma omp for
        #endif
        for(int iElm = 0; iElm < cnc->getNumElements(); ++iElm)
        {
          // if(plot) {
          //   if(iElm > 20) plot = false;
          //   // Plot edge normal
          //   for(int j = 0; j < _nNodePerElem; ++j)
          //   {
          //     vtkNode &n0 = _vtkNodes[iElm * 3 + j];
          //     vtkNode &n1 = _vtkNodes[iElm * 3 + ((j + 1) % 3)];
          //     double mid[2] = {(n0.pos[0]+n1.pos[0])/2., (n0.pos[1]+n1.pos[1])/2.};
          //     double t[2] = {n1.pos[0] - n0.pos[0],
          //                    n1.pos[1] - n0.pos[1]};
          //     fprintf(myFile, "VP(%g,%g,%g){%g,%g,%g};\n", n0.pos[0], n0.pos[1], 0.,
          //       t[0], t[1], 0.);
          //     double normal[2] = {t[1], -t[0]};
          //     double norm = sqrt(normal[0]*normal[0]+normal[1]*normal[1]);
          //     double fac = 0.25;
          //     fprintf(myFile, "VP(%g,%g,%g){%g,%g,%g};\n", mid[0], mid[1], 0.,
          //       fac*normal[0]/norm, fac*normal[1]/norm, 0.);
          //   }
          // }

          // cnc->getElementTransformation(iElm, T);
          // #define ntarget 4
          // // int targetEdge[ntarget] = {66, 88, 248, 300};
          // int targetEdge[ntarget] = {248, 300, 315, 330};
          // int sign = 1;
          // int whichEdge, whichDof;
          // bool found = false;
          // for(int it = 0; it < ntarget; ++it) {

          //   for(int iedge = 0; iedge < 3; ++iedge) {
          //     int edge = fabs(cnc->getEdgeConnectivity(iElm, iedge)) - 1;
          //     if(edge == targetEdge[it]) {
          //       // whichEdge = ((iedge+2)%3)+1;
          //       // whichEdge = (iedge+2)%3;
          //       whichEdge = iedge;
          //       whichDof = (iedge+2)%3;
          //       found = true;
          //       if(!firstForEdge[it]){
          //         sign = -1;
          //         firstForEdge[it] = true;
          //       }
          //       break;
          //     }
          //   }
          // }
          // if(!found) {
          //   for(int j = 0; j < _nNodePerElem; ++j)
          //   {
          //     int globalTag = iElm * _nNodePerElem + j;
          //     vtkNode &node = _vtkNodes[globalTag];
          //     for(int iComp = 0; iComp < nComponentsToWrite; ++iComp) {
          //       node.data[iField][iComp] = 0.;
          //     }
          //   }
          //   continue;
          // }

          // {
          //   vtkNode &n0 = _vtkNodes[iElm * 3 + 0];
          //   vtkNode &n1 = _vtkNodes[iElm * 3 + 1];
          //   vtkNode &n2 = _vtkNodes[iElm * 3 + 2];
          //   feInfo("Elm %d : (%f,%f) - (%f,%f) - (%f,%f)", iElm,
          //     n0.pos[0], n0.pos[1],
          //     n1.pos[0], n1.pos[1],
          //     n2.pos[0], n2.pos[1]);
          // }

          // feInfo("sign = %d", sign);
          // feInfo("whichEdge = %d", whichEdge);
          // feInfo("whichDof  = %d", whichDof);
          // feInfo("%+-f,  %+-f", T.dxdr[0], T.dxds[0]);
          // feInfo("%+-f,  %+-f", T.dxdr[1], T.dxds[1]);
          ///////////////////////////

          for(int j = 0; j < _nNodePerElem; ++j)
          {
            int globalTag = iElm * _nNodePerElem + j;
            vtkNode &node = _vtkNodes[globalTag];

            // Evaluate field at visualization node
            space->initializeAddressingVector(iElm, adr);
            for(size_t k = 0; k < adr.size(); ++k) { solOnElem[k] = solArray[adr[k]]; }

            if(_addP2Nodes) {
              space->interpolateVectorField(solOnElem, nComponents, referenceCoordinatesQuadraticTri[j], res);
            } else
            {
              // if(space->getElementType() == elementType::RAVIART_THOMAS)
              // {
              //   space->interpolateVectorFieldRT(T, sign, whichDof, solOnElem, referenceCoordinatesTri[j], res);
              // } else {
                space->interpolateVectorField(solOnElem, nComponents, referenceCoordinatesTri[j], res);
              // }
            }

            for(int iComp = 0; iComp < nComponentsToWrite; ++iComp) {
              node.data[iField][iComp] = res[iComp];
            }
          }
        }
        ///////////////////////////
        // fprintf(myFile, "};\n");
        // fclose(myFile);
        ///////////////////////////
      }
    } // for spacesToExport
  } // if _recreateVTKNodes

  return FE_STATUS_OK;
}

// Write the mesh from VTK nodes
void feExporterVTK::writeMesh(std::ostream &output)
{
  output << "DATASET UNSTRUCTURED_GRID\n";
  output << "POINTS " << _vtkNodes.size() << " double\n";

  for(auto node : _vtkNodes) {
    output << node.pos[0] << " " << node.pos[1] << " " << node.pos[2] << std::endl;
  }

  int nElm = 0;
  for(auto *cnc : _exportedConnectivities) {
    nElm += cnc->getNumElements();
  }

  //
  // Write elements connectivities for each cnc of maximum dimension
  //
  output << "CELLS " << nElm << " " << nElm * (_nNodePerElem + 1) << std::endl;

  for(auto *cnc : _exportedConnectivities)
  {
    int numVerticesPerElem = cnc->getNumVerticesPerElem();

    for(int iElm = 0; iElm < cnc->getNumElements(); ++iElm)
    {
      output << _nNodePerElem << " ";

      if(_addP2Nodes) {
        // Mid-points were added and do not exist in the mesh
        // Vertices
        for(int j = 0; j < 3; ++j) {
          int node = cnc->getVertexConnectivity(iElm, j);
          output << node << " ";
        }
        // Edges
        for(int j = 0; j < 3; ++j) {
          int iEdge = fabs(cnc->getEdgeConnectivity(iElm, j));
          int node = _edgeGlobalTag[iEdge];
          output << node << " ";
        }
      } else {
        // Regular case : all nodes exist in the mesh
        for(int j = 0; j < numVerticesPerElem; ++j) {
          if(_discontinuousMesh) {
            int node = iElm * numVerticesPerElem + j;
            output << node << " ";
          } else {
            int node = cnc->getVertexConnectivity(iElm, j);
            output << node << " ";
          }
        }
      }
      output << std::endl;
    }
  }

  // Write element types
  output << "CELL_TYPES " << nElm << std::endl;

  for(auto *cnc : _exportedConnectivities) {
    int vtkElem = cncToVTKmap[cnc->getInterpolant()];
    for(int iElm = 0; iElm < cnc->getNumElements(); ++iElm) {
      if(_addP2Nodes)
        output << VTK_QUADRATIC_TRIANGLE << std::endl;
      else
        output << vtkElem << std::endl;
    }
  }
}

//
// Write discontinuous mesh.
// The discontinuous VTK nodes were created element by element,
// rather than all P1 vertices followed by additional P2 vertices,
// so the mesh must be written accordingly.
//
void feExporterVTK::writeDiscontinuousMesh(std::ostream &output)
{
  output << "DATASET UNSTRUCTURED_GRID\n";
  output << "POINTS " << _vtkNodes.size() << " double\n";

  for(auto node : _vtkNodes) {
    output << node.pos[0] << " " << node.pos[1] << " " << node.pos[2] << std::endl;
  }

  int nElm = 0;
  for(auto *cnc : _exportedConnectivities) {
    nElm += cnc->getNumElements();
  }

  //
  // Write elements connectivities for each cnc of maximum dimension.
  // The VTK nodes are created in order element by element,
  // so the connectivity is trivially sequential
  //
  output << "CELLS " << nElm << " " << nElm * (_nNodePerElem + 1) << std::endl;

  int cnt = 0;
  for(auto *cnc : _exportedConnectivities) {
    for(int iElm = 0; iElm < cnc->getNumElements(); ++iElm) {
      output << _nNodePerElem << " ";
      for(int j = 0; j < _nNodePerElem; ++j) {
        output << cnt++ << " ";
      }
      output << std::endl;
    }
  }

  // Write element types
  output << "CELL_TYPES " << nElm << std::endl;

  for(auto *cnc : _exportedConnectivities) {
    int vtkElem = cncToVTKmap[cnc->getInterpolant()];
    for(int iElm = 0; iElm < cnc->getNumElements(); ++iElm) {
      if(_addP2Nodes)
        output << VTK_QUADRATIC_TRIANGLE << std::endl;
      else
        output << vtkElem << std::endl;
    }
  }
}

void feExporterVTK::writeVTKNodes(std::ostream &output,
                                  std::unordered_map<std::string, std::pair<int, int> > &fields)
{
  for(const auto &pair : fields) {
    std::string name = pair.first;
    int tag = pair.second.first;
    int numComponents = pair.second.second;

    // if(_discontinuousMesh) {
      // Force vectors to have 3 components to use Paraview vector features
      if(numComponents > 1) numComponents = 3;
    // }

    output << name << " " << numComponents << " " << _vtkNodes.size() << " double" << std::endl;

    for(const auto &node : _vtkNodes) {
      for(int iComp = 0; iComp < numComponents; ++iComp) {
        output << node.data[tag][iComp] << std::endl;
      }
    }
  }
}

//
// Write the solution on the mesh at current time.
// Only implemented for triangles and tetrahedra.
// For now, only a single highest dimensional connectivity is exported, i.e.,
// it does not treat meshes with two domains.
//
feStatus feExporterVTK::writeStep(std::string fileName,
                                  const std::vector<std::string> &requiredFields)
{
#if defined(HAVE_MPI)
  // Write VTK file from master rank only for now
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0) {
#endif
  tic();
  std::filebuf fb;
  if(fb.open(fileName, std::ios::out)) {
    std::ostream output(&fb);
    writeHeader(output);

    // For now, we only print spaces of maximal dimension (non boundary)
    std::vector<feSpace *> spacesToExport;
    std::set<std::string> cncToExport;
    std::set<std::string> fieldsToExport;
    std::unordered_map<std::string, int> fieldTags;

    std::unordered_map<std::string, std::pair<int, int> > fields;

    for(feSpace *fS : _spaces) {

      // If prescribed fields are given, only print those fields
      if(requiredFields.size() > 0) {
        bool found = false;
        for(auto name : requiredFields) {
          if(name == fS->getFieldID()) {
            found = true;
            break;
          }
        }
        // Skip this field if not required
        if(!found) continue;
      }

      if(fS->getDim() == _mesh->getDim()) {
        spacesToExport.push_back(fS);
        cncToExport.insert(fS->getCncGeoID());
        fieldsToExport.insert(fS->getFieldID());
        _exportedConnectivities.insert(fS->getCncGeo());

        std::pair<std::string, std::pair<int, int> > f;
        f.first = fS->getFieldID();
        // f.second.first = fS->getFieldTag(); // Not assigned in feSpace...
        f.second.second = fS->getNumComponents();
        fields.insert(f);
      }
    }

    // Check if any field is discontinuous
    for(feSpace *s : spacesToExport) {
      if(s->isDiscontinuous()) {
        _discontinuousMesh = true;
        break;
      }
    }

    // if(_exportedConnectivities.size() > 1) {
    //   return feErrorMsg(FE_STATUS_ERROR,
    //     "Currently only a single Physical Entity can be exported to a VTK file");
    // }
    const feCncGeo *cnc = *(_exportedConnectivities.begin());

    int cnt = 0;
    for(auto field : fieldsToExport) {
      fieldTags[field] = cnt;
      fields[field].first = cnt++;
    }

    feInfoCond(FE_VERBOSE > 0, "\t\tExporting geometric connectivity \"%s\" to file \"%s\"",
               cnc->getID().c_str(), fileName.c_str());

    geometricInterpolant interpolant = cnc->getInterpolant();

    if(interpolant != geometricInterpolant::TRIP1 &&
       interpolant != geometricInterpolant::TRIP2 &&
       interpolant != geometricInterpolant::TETP1)
    {
      return feErrorMsg(FE_STATUS_ERROR,
                        "Only triangles and tetrahedra can be exported to a VTK file");
    }

    /* Although VTK_HIGHER_ORDER_TRIANGLE and VTK_LAGRANGE_TRIANGLE
    exist in the VTK documentation, I can't find tutorials on how to use them...
    So for now the highest order triangle is the VTK_QUADRATIC_TRIANGLE.
    Quadratic triangles are exported if :
      - the geometry is P2 (curved elements)
      - any field has a P2 or greater interpolant
    In the second case, the mesh is not made of 6-node triangles, so we need to create them
    and evaluate the fields at those mid-edge nodes. */
    int geometryPolynomialDegree = getGeometricInterpolantDegree(interpolant);
    int highestFieldPolynomialDegree = 0;
    for(feSpace *fS : spacesToExport)
      highestFieldPolynomialDegree = fmax(highestFieldPolynomialDegree, fS->getPolynomialDegree());

    _addP2Nodes = (geometryPolynomialDegree == 1) && (highestFieldPolynomialDegree >= 2);

    _nNodePerElem = verticesPerElem.at(interpolant);
    _nEdgePerElem = edgePerElem.at(interpolant);

    if(_addP2Nodes) {
      // Add one visualization vertex per edge
      _nNodePerElem += _nEdgePerElem;
    }

    // Create vtkNodes and write mesh
    if(_discontinuousMesh) {
      feStatus s = createDiscontinuousVTKNodes(spacesToExport, fieldTags);
      if(s != FE_STATUS_OK) { return s; }
      writeDiscontinuousMesh(output);
    } else {
      feStatus s = createVTKNodes(spacesToExport, fieldTags);
      if(s != FE_STATUS_OK) { return s; }
      writeMesh(output);
    }

    // Write the fields
    output << "POINT_DATA " << _vtkNodes.size() << std::endl;
    output << "FIELD FieldData " << fieldsToExport.size() << std::endl;
    writeVTKNodes(output, fields);

    fb.close();
  } else {
    return feErrorMsg(FE_STATUS_FILE_CANNOT_BE_OPENED, "Could not open output file.");
  }

  feInfoCond(FE_VERBOSE > 0, "\t\tWrote VTK file in %f s", toc());

  #if defined(HAVE_MPI)
  }
  #endif

  return FE_STATUS_OK;
}

feExporterVTK::feExporterVTK(feMesh *mesh, feSolution *sol, feEigenProblem *eigenProblem,
                             feMetaNumber *metaNumber, const std::vector<feSpace *> &feSpaces)
  : feExporter(mesh, sol, metaNumber, feSpaces)
{
  _eigenProblem = eigenProblem;
}

void feExporterVTK::writeEigenvector(std::ostream &output, const feCncGeo *cnc, feSpace *intSpace,
                                     std::string fieldID, int eigenPairIndex, size_t nEigenPairs,
                                     eigenPair &ep, bool loopOverCnc)
{
#if defined(HAVE_PETSC)
  std::vector<double> &solVec = _sol->getSolutionReference();
  std::vector<feInt> adr(intSpace->getNumFunctions());
  std::vector<double> sol(intSpace->getNumFunctions());
  int nVertices = loopOverCnc ? cnc->getNumVertices() : _mesh->getNumVertices();
  feNumber *n = _metaNumber->getNumbering(fieldID);

  size_t n_zero = ceil(log10(nEigenPairs + 1));
  std::string notPadded = std::to_string(eigenPairIndex + 1);
  std::string padded = std::string(n_zero - std::min(n_zero, notPadded.length()), '0') + notPadded;

  // Name of the field is the number of the mode followed by the real part of its eigenvalue
  output << "Mode" + padded + "_" + std::to_string(ep.valReal) << " 1 " << _writtenNodes
         << " double" << std::endl;

  int iDOF, elm;
  Vertex *v;
  std::vector<double> x;
  double r[3];

  PetscScalar *vecRealArray;
  VecGetArray(ep.vecReal, &vecRealArray);

  for(int iVertex = 0; iVertex < nVertices; ++iVertex) {
    iDOF = n->getDOFNumberAtVertex(iVertex);

    if(iDOF >= 0) {
      // There is a degree of freedom at this mesh vertex

      if(n->getDOFCodeAtVertex(iVertex) == DOF_ESSENTIAL) {
        // DOF is an essential BC: assign the corresponding BC stored in the feSolution
        output << solVec[iDOF] << std::endl;
      }

      if(n->getDOFCodeAtVertex(iVertex) == DOF_UNKNOWN) {
        // True DOF: read value from the eigenvector (computed at true DOFs only)
        // feInfo("Reading ev %d at DOF %d = %f", eigenPairIndex, iDOF, vecRealArray[iDOF]);
        output << vecRealArray[iDOF] << std::endl;
      }

    } else {
      // feInfo("TEST");
      /* No dof associated to the mesh vertex.
      Interpolate solution at vertex. */
      v = _mesh->getVertex(iVertex);
      x = {v->x(), v->y(), v->z()};
      _mesh->locateVertex(x.data(), elm, r);
      intSpace->initializeAddressingVector(elm, adr);

      // // Initialize solution
      // for(size_t i = 0; i < adr.size(); ++i) {
      //   sol[i] = solVec[adr[i]];
      // }

      // if(intSpace->useGlobalFunctions()) {
      //   val = intSpace->interpolateField(sol, elm, x);
      // } else {
      //   val = intSpace->interpolateField(sol, r);
      // }

      // output << val << std::endl;
    }
  }

  /* Write the additional P2 nodes. Interpolation is not required if the
  field is quadratic, but it's easier to just interpolate for all fields. */
  if(_addP2Nodes) {
    /* Additional nodes are added in the order of the edges,
    like in feExporterVTK::writeNodes. Both must be modified together. */
    Vertex *v0;
    Vertex *v1;
    double val;
    for(auto e : _mesh->_edges) {
      v0 = e.getVertex(0);
      v1 = e.getVertex(1);
      x = {(v0->x() + v1->x()) / 2., (v0->y() + v1->y()) / 2., (v0->z() + v1->z()) / 2.};
      _mesh->locateVertex(x.data(), elm, r);
      intSpace->initializeAddressingVector(elm, adr);
      for(size_t i = 0; i < adr.size(); ++i) sol[i] = vecRealArray[adr[i]];
      if(intSpace->useGlobalFunctions()) {
        val = intSpace->interpolateField(sol, elm, x);
      } else {
        val = intSpace->interpolateField(sol, r);
      }
      output << val << std::endl;
    }
  }

  VecRestoreArray(ep.vecReal, &vecRealArray);

#else
  UNUSED(output, cnc, intSpace, fieldID, eigenPairIndex, nEigenPairs, ep, loopOverCnc);
#endif
}

// For block by block documentation, see the feExporterVTK::writeStep function.
feStatus feExporterVTK::writeEigenvectors(std::string fileName)
{
  std::filebuf fb;
  if(fb.open(fileName, std::ios::out)) {
    std::ostream output(&fb);
    writeHeader(output);

    std::vector<feSpace *> spacesToExport;
    std::set<std::string> cncToExport;
    for(feSpace *fS : _spaces) {
      if(fS->getDim() == _mesh->getDim()) {
        spacesToExport.push_back(fS);
        cncToExport.insert(fS->getCncGeoID());
      }
    }

    if(cncToExport.size() > 1) {
      feWarning("Multiple domains visualization is not implemented. Only the first domain will be "
                "exported.");
    }

    const feCncGeo *cnc = nullptr;
    for(feSpace *fS : spacesToExport) {
      if(fS->getCncGeoID() == *cncToExport.begin()) {
        cnc = fS->getCncGeo();
        break;
      }
    }

    feInfoCond(FE_VERBOSE > 0, "\t\tExporting geometric connectivity \"%s\" to file \"%s\"",
               cnc->getID().c_str(), fileName.c_str());

    int geometryPolynomialDegree;
    if(cnc->getInterpolant() == geometricInterpolant::TRIP1) {
      geometryPolynomialDegree = 1;
    } else if(cnc->getInterpolant() == geometricInterpolant::TRIP2) {
      geometryPolynomialDegree = 2;
    } else {
      return feErrorMsg(FE_STATUS_WRITE_ERROR,
                        "Only P1 and P2 triangles can be exported to a VTK file...");
    }

    int highestFieldPolynomialDegree = 0;
    for(feSpace *fS : spacesToExport)
      highestFieldPolynomialDegree = fmax(highestFieldPolynomialDegree, fS->getPolynomialDegree());

    _addP2Nodes = (geometryPolynomialDegree == 1) && (highestFieldPolynomialDegree >= 2);

    size_t nEigenPairs = _eigenProblem->getNumConvergedPairs();

    writeNodes(output);
    writeElementsConnectivity(output);

    output << "POINT_DATA " << _writtenNodes << std::endl;
    output << "FIELD FieldData " << nEigenPairs << std::endl;

    // Write each eigenvector for the current field
    eigenPair ep;
    if(spacesToExport.size() > 1) {
      feWarning("Currently exporting the eigenvectors for a single field only.");
    }
    for(feSpace *fS : spacesToExport) {
      for(size_t i = 0; i < nEigenPairs; ++i) {
        ep = _eigenProblem->getEigenpair(i);
        writeEigenvector(output, cnc, fS, fS->getFieldID(), i, nEigenPairs, ep, false);
      }
    }

    fb.close();
  } else {
    return feErrorMsg(FE_STATUS_FILE_CANNOT_BE_OPENED, "Could not open output file.");
  }

  return FE_STATUS_OK;
}