#include "feExporter.h"
#include "feVertex.h"
#include "feCncGeo.h"

#include <iostream>
#include <fstream>
#include <map>

extern bool FE_VERBOSE;

static int VTK_VERTEX = 1;
static int VTK_LINE = 3;
static int VTK_TRIANGLE = 5;
// static int VTK_QUAD = 9;
// static int VTK_TETRA = 10;
// static int VTK_HEXAHEDRON = 12;
// static int VTK_PYRAMID = 14;
static int VTK_QUADRATIC_EDGE = 21;
static int VTK_QUADRATIC_TRIANGLE = 22;
// static int VTK_QUADRATIC_QUAD = 23;
// static int VTK_QUADRATIC_TETRA = 24;
// static int VTK_QUADRATIC_HEXAHEDRON = 25;

// A etoffer au fur et a mesure
std::map<geometricInterpolant, int> cncToVTKmap = {
  {geometricInterpolant::POINTP0, VTK_VERTEX},
  {geometricInterpolant::LINEP1, VTK_LINE},
  {geometricInterpolant::TRIP1, VTK_TRIANGLE},
  {geometricInterpolant::LINEP2, VTK_QUADRATIC_EDGE},
  {geometricInterpolant::TRIP2, VTK_QUADRATIC_TRIANGLE}};

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
  for(auto *cnc : _mesh->getCncGeo()){
    if(cnc->getDim() == 2){
      nElm += cnc->getNumElements();
    }
  }

  int nNodePerElem = _addP2Nodes ? 6 : 3;

  // Write elements connectivities for each cnc wih dim = 2
  output << "CELLS " << nElm << " " << nElm * (nNodePerElem + 1) << std::endl;

  for(auto *cnc : _mesh->getCncGeo()){
    if(cnc->getDim() == 2){
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

  for(auto *cnc : _mesh->getCncGeo()){
    if(cnc->getDim() == 2){
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

// void feExporterVTK::writeField(std::ostream &output, feCncGeo *cnc, feSpace *intSpace,
//                                std::string fieldID, bool loopOverCnc)
// {
//   std::vector<double> &solVec = _sol->getSolutionReference();
//   std::vector<feInt> adr(intSpace->getNumFunctions());
//   std::vector<double> sol(intSpace->getNumFunctions());
//   // int nVertices = loopOverCnc ? cnc->getNumVertices() : _mesh->getNumVertices();
//   int nVertices = _mesh->getNumVertices();
//   feNumber *n = _metaNumber->getNumbering(fieldID);
//   int nComponents = intSpace->getNumComponents();

//   // output << "POINT_DATA " << _writtenNodes << std::endl;
//   // output << "SCALARS " << fieldID << " double 1" << std::endl;
//   // output << "LOOKUP_TABLE default" << std::endl;

//   output << fieldID << " " << nComponents << " " << _writtenNodes << " double" << std::endl;

//   int iDOF, elm;
//   Vertex *v;
//   std::vector<double> x;
//   std::vector<double> res(3, 0.);
//   double r[3], val;

//   // Write field(s) to a text file
//   std::string fileName = "solution" + fieldID + ".txt";
//   FILE *f = fopen(fileName.c_str(), "w");
//   for(int iVertex = 0; iVertex < nVertices; ++iVertex) {
//     // Loop over components for vector-valued FE spaces
//     for(int iComp = 0; iComp < nComponents; ++iComp) {
//       iDOF = n->getDOFNumberAtVertex(iVertex, iComp);

//       if(iDOF >= 0) {
//         // There is a degree of freedom at this mesh vertex

//         // FIXME: This is only valid for Lagrange type elements,
//         // where to DOF is the function evaluation. We should interpolate.
//         output << solVec[iDOF] << std::endl;

//       } else {
//         /* No dof associated to the mesh vertex.
//         Interpolate solution at vertex. */
//         v = _mesh->getVertex(iVertex);
//         x = {v->x(), v->y(), v->z()};
//         _mesh->locateVertex(x.data(), elm, r);
//         intSpace->initializeAddressingVector(elm, adr);

//         for(size_t i = 0; i < adr.size(); ++i) {
//           sol[i] = solVec[adr[i]];
//         }

//         if(intSpace->useGlobalFunctions()) {
//           val = intSpace->interpolateField(sol, elm, x);
//         } else {
//           val = intSpace->interpolateField(sol, r);
//         }

//         output << val << std::endl;
//       }
//       fprintf(f, "%+-16.16e\n", solVec[iDOF]);
//     }
//   }
//   fclose(f);

//   /* Write the additional P2 nodes. Interpolation is not required if the
//   field is quadratic, but it's easier to just interpolate for all fields. */
//   if(_addP2Nodes) {
//     /* Additional nodes are added in the order of the edges,
//     like in feExporterVTK::writeNodes. Both must be modified together. */
//     Vertex *v0, *v1;
//     for(auto e : _mesh->_edges) {
//       v0 = e.getVertex(0);
//       v1 = e.getVertex(1);
//       x = {(v0->x() + v1->x()) / 2., (v0->y() + v1->y()) / 2., (v0->z() + v1->z()) / 2.};
//       _mesh->locateVertex(x.data(), elm, r);
//       intSpace->initializeAddressingVector(elm, adr);
//       for(size_t i = 0; i < adr.size(); ++i) sol[i] = solVec[adr[i]];
//       if(intSpace->useGlobalFunctions()) {
//         val = intSpace->interpolateField(sol, elm, x);
//         output << val << std::endl;
//       } else {
//         if(nComponents > 1) {
//           intSpace->interpolateVectorField(sol, nComponents, r, res);
//           for(int i = 0; i < nComponents; ++i) output << res[i] << std::endl;
//         } else {
//           val = intSpace->interpolateField(sol, r);
//           // output << val << std::endl;
//           output << 1.17 << std::endl;
//         }
//       }
//     }
//   }
// }

void feExporterVTK::writeField(std::ostream &output, feCncGeo *cnc, feSpace *intSpace,
                               std::string fieldID, bool loopOverCnc)
{
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
  // std::string fileName = "solution" + fieldID + ".txt";
  // FILE *f = fopen(fileName.c_str(), "w");

  // CHAQUE SOMMET DOIT SAVOIR SUR QUELLE CNC IL EST POUR INTERPOLER

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
      // fprintf(f, "%+-16.16e\n", solVec[iDOF]);
    }
  }
  // fclose(f);

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

feStatus feExporterVTK::createVTKNodes(std::vector<feSpace *> &spacesToExport, std::unordered_map<std::string, int> &fieldTags)
{
  if(_addP2Nodes)
    _vtkNodes.resize(_mesh->getNumVertices() + _mesh->getNumEdges());
  else
    _vtkNodes.resize(_mesh->getNumVertices());

  // Create nodes with tag and coordinates
  for(auto *cnc : _exportedConnectivities) {
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

  int cnt = _mesh->getNumVertices();
  // Add visualization nodes at mid-edges
  if(_addP2Nodes) {
    for(auto e : _mesh->_edges)
    {
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

  int iDOF, elm;
  std::vector<double> x;
  std::vector<double> res(3, 0.);
  double r[3], val;

  std::vector<double> &solVec = _sol->getSolutionReference();

  // Set vtkNodes data
  for(auto *space : spacesToExport) {

    std::string field = space->getFieldID();
    int iField = fieldTags[field];
    int nComponents = space->getNumComponents();
    feNumber *numbering = _metaNumber->getNumbering(field);
    std::vector<feInt> adr(space->getNumFunctions());
    std::vector<double> sol(space->getNumFunctions());

    feInfo("Exporting field %s on connectivity %s", field.data(), space->getCncGeoID().data());

    const feCncGeo *cnc = space->getCncGeo();

    // Loop over mesh vertices on this connectivity
    // and evaluate field
    for(int i = 0; i < cnc->getNumVertices(); ++i) {

      int globalTag = cnc->getUniqueVertexConnectivity(i);

      // Loop over components for vector-valued FE spaces
      for(int iComp = 0; iComp < nComponents; ++iComp)
      {
        iDOF = numbering->getDOFNumberAtVertex(globalTag, iComp);

        if(iDOF >= 0) {
          // There is a degree of freedom at this mesh vertex

          // FIXME: This is only valid for Lagrange type elements,
          // where to DOF is the function evaluation. We should interpolate.
          _vtkNodes[globalTag].data[iField][iComp] = solVec[iDOF];

        } else {
          /* No dof associated to the mesh vertex.
          Interpolate solution at vertex. */
          vtkNode node = _vtkNodes[globalTag];
          x = {node.pos[0], node.pos[1], node.pos[2]};
          _mesh->locateVertex(x.data(), elm, r);
          space->initializeAddressingVector(elm, adr);

          for(size_t i = 0; i < adr.size(); ++i) {
            sol[i] = solVec[adr[i]];
          }

          if(space->useGlobalFunctions()) {
            val = space->interpolateField(sol, elm, x);
          } else {
            val = space->interpolateField(sol, r);
          }

          _vtkNodes[globalTag].data[iField][iComp] = val;
        }
      }
    }

    /* Write the additional P2 nodes. Interpolation is not required if the
    field is quadratic, but it's easier to just interpolate for all fields. */
    if(_addP2Nodes) {
      // for(auto e : _mesh->_edges)
      for(int iEdge = 0; iEdge < cnc->getNumEdges(); ++iEdge)
      {
        int edgeTag = fabs(cnc->getUniqueEdgeConnectivity(iEdge));
        const Edge *e = _mesh->getEdge(edgeTag);
        int globalTag = _edgeGlobalTag[e->getTag()];
        const Vertex *v0 = e->getVertex(0);
        const Vertex *v1 = e->getVertex(1);
        x = { (v0->x() + v1->x()) / 2.,
              (v0->y() + v1->y()) / 2.,
              (v0->z() + v1->z()) / 2. };

        bool returnLocalElmTag = true;
        std::string targetConnectivity = cnc->getID();
        bool ret = _mesh->locateVertex(x.data(), elm, r, 1e-5, returnLocalElmTag, targetConnectivity);

        if(!ret){
          return feErrorMsg(FE_STATUS_ERROR, "Could not find edge vertex (%f,%f,%f) in the mesh", x[0], x[1], x[2]);
        }

        feInfo("%f - %f - %f", r[0], r[1], r[2]);

        space->initializeAddressingVector(elm, adr);

        for(size_t i = 0; i < adr.size(); ++i)
          sol[i] = solVec[adr[i]];

        if(space->useGlobalFunctions())
        {
          val = space->interpolateField(sol, elm, x);
          // FIXME: vector field with global functions?
          _vtkNodes[globalTag].data[iField][0] = val;
        }
        else {
          if(nComponents > 1)
          {
            space->interpolateVectorField(sol, nComponents, r, res);
            for(int i = 0; i < nComponents; ++i){
              _vtkNodes[globalTag].data[iField][i] = res[i];
            }
          } else {
            val = space->interpolateField(sol, r);
            _vtkNodes[globalTag].data[iField][0] = val;
            feInfo("Found midedge (%f,%f) of %d on elm %d on cnc %s : val = %f", x[0], x[1], e->getTag(), elm, cnc->getID().data(), val);
          }
        }
      }
    }
  }
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
  for(auto *cnc :_exportedConnectivities){
    nElm += cnc->getNumElements();
  }

  int nNodePerElem = _addP2Nodes ? 6 : 3;

  // Write elements connectivities for each cnc wih dim = 2
  output << "CELLS " << nElm << " " << nElm * (nNodePerElem + 1) << std::endl;

  for(auto *cnc : _exportedConnectivities){
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

  // Write element types
  output << "CELL_TYPES " << nElm << std::endl;

  for(auto *cnc : _exportedConnectivities){
    int vtkElem = cncToVTKmap[cnc->getInterpolant()];
    for(int iElm = 0; iElm < cnc->getNumElements(); ++iElm) {
      if(_addP2Nodes)
        output << VTK_QUADRATIC_TRIANGLE << std::endl;
      else
        output << vtkElem << std::endl;
    }
  }
}

void feExporterVTK::writeVTKNodes(std::ostream &output, std::unordered_map<std::string, std::pair<int, int>> &fields)
{
  for(auto pair : fields) {

    std::string name = pair.first;
    int tag = pair.second.first;
    int numComponents = pair.second.second;

    output << name << " " << numComponents << " " << _vtkNodes.size() << " double" << std::endl;

    for(auto node : _vtkNodes) {
      for(int iComp = 0; iComp < numComponents; ++iComp) {
        output << node.data[tag][iComp] << std::endl;
      }
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
    std::set<std::string> fieldsToExport;
    std::unordered_map<std::string, int> fieldTags;

    std::unordered_map<std::string, std::pair<int, int>> fields;

    for(feSpace *fS : _spaces) {
      if(fS->getDim() == 2) {
        spacesToExport.push_back(fS);
        cncToExport.insert(fS->getCncGeoID());
        fieldsToExport.insert(fS->getFieldID());
        _exportedConnectivities.insert(fS->getCncGeo());

        std::pair<std::string, std::pair<int, int>> f;
        f.first = fS->getFieldID();
        // f.second.first = fS->getFieldTag(); // Not assigned in feSpace...
        f.second.second = fS->getNumComponents();
        fields.insert(f);
      }
    }

    // Grab the connectivity from any matching fespace
    feCncGeo *cnc;
    for(feSpace *fS : spacesToExport) {
      if(fS->getCncGeoID() == *cncToExport.begin()) {
        cnc = fS->getCncGeo();
        break;
      }
    }

    int cnt = 0;
    for(auto field : fieldsToExport) {
      fieldTags[field] = cnt++;
      fields[field].first = fieldTags[field];
    }

    feInfoCond(FE_VERBOSE > 0, "\t\tExporting geometric connectivity \"%s\" to file \"%s\"",
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

    // ///////////////////////////////////////
    // // Write nodes and elements
    // writeNodes(output);
    // writeElementsConnectivity(output);

    // output << "POINT_DATA " << _writtenNodes << std::endl;
    // output << "FIELD FieldData " << fieldsToExport.size() << std::endl;
    // // Write the field associated to each fespace in spacesToExport
    // for(feSpace *fS : spacesToExport) {
    //   feInfo("Exporting %s - %s", fS->getFieldID().data(), fS->getCncGeoID().data());
    //   writeField(output, cnc, fS, fS->getFieldID(), false);
    //   break;
    // }
    // // writeField(output, cnc, fS, fS->getFieldID(), false);
    // ////////////////////////////////////////


    ////////////////////////////////
    feStatus s = createVTKNodes(spacesToExport, fieldTags);
    if(s != FE_STATUS_OK) {
      return s;
    }

    // for(auto node : _vtkNodes){
    //   feInfo("Created node %d : (%f,%f,%f) - data = [%f,%f,%f], [%f,%f,%f], [%f,%f,%f]",
    //     node.globalTag,
    //     node.pos[0], node.pos[1], node.pos[2],
    //     node.data[0][0], node.data[0][1], node.data[0][2], 
    //     node.data[1][0], node.data[1][1], node.data[1][2], 
    //     node.data[2][0], node.data[2][1], node.data[2][2]
    //   );
    // }

    writeMesh(output);
    output << "POINT_DATA " << _vtkNodes.size() << std::endl;
    output << "FIELD FieldData " << fieldsToExport.size() << std::endl;
    writeVTKNodes(output, fields);
    ////////////////////////////////

    fb.close();
  } else {
    return feErrorMsg(FE_STATUS_FILE_CANNOT_BE_OPENED, "Could not open output file.");
  }

  return FE_STATUS_OK;
}

feExporterVTK::feExporterVTK(feMesh *mesh, feSolution *sol, feEigenProblem *eigenProblem,
                             feMetaNumber *metaNumber, const std::vector<feSpace *> &feSpaces)
  : feExporter(mesh, sol, metaNumber, feSpaces)
{
  _eigenProblem = eigenProblem;
};

void feExporterVTK::writeEigenvector(std::ostream &output, feCncGeo *cnc, feSpace *intSpace,
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
    int elm;
    double val;
    std::vector<double> x;
    double r[3];
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

    feCncGeo *cnc;
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