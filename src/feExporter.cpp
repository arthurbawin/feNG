#include "feExporter.h"
#include "feVertex.h"

#include <iostream>
#include <fstream>
#include <map>

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
std::map<std::string,int> cncToVTKmap = {
  {"Point0D", VTK_VERTEX            },
  {"LineP1" , VTK_LINE              },
  {"TriP1"  , VTK_TRIANGLE          },
  {"LineP2" , VTK_QUADRATIC_EDGE    },
  {"TriP2"  , VTK_QUADRATIC_TRIANGLE}
};

void feExporterVTK::writeHeader(std::ostream& output){
  output << "# vtk DataFile Version 4.2\n";
  output << "test output\n";
  output << "ASCII\n";
}

void feExporterVTK::writeNodes(std::ostream& output, feCncGeo *cnc){
  int nNod = cnc->getNbNodes();
  output << "DATASET UNSTRUCTURED_GRID\n";
  output << "POINTS " << nNod << " double\n";
  Vertex *v;
  for(int i = 0; i < nNod; ++i){
    v = _mesh->getVertex(i);
    output << v->x() << " " << v->y() << " " << v->z() << std::endl;
  }
}

void feExporterVTK::writeElementsConnectivity(std::ostream& output, feCncGeo *cnc){
  int nElm = cnc->getNbElm();
  int nNodePerElem = cnc->getNbNodePerElem();
  output << "CELLS " << nElm << " " << nElm*(nNodePerElem+1) << std::endl;
  for(int iElm = 0; iElm < nElm; ++iElm){
    output << nNodePerElem << " ";
    for(int j = 0; j < nNodePerElem; ++j){
      int node = cnc->getNodeConnectivity(iElm,j);
      output << node << " ";
    }
    output << std::endl;
  }
  output << "CELL_TYPES " << nElm << std::endl;
  int vtkElem = cncToVTKmap[cnc->getForme()];
  for(int iElm = 0; iElm < nElm; ++iElm){
    output << vtkElem << std::endl;
  }
}

void feExporterVTK::writeField(std::ostream& output, feCncGeo *cnc, std::string fieldID){
  std::vector<double>& sol = _sol->getSolutionReference();
  int nNod = cnc->getNbNodes();
  feNumber *n = _metaNumber->getNumbering(fieldID);
  // int nNod = n->getNbNodes();
  output << "SCALARS " << fieldID << " double 1" << std::endl;
  output << "LOOKUP_TABLE default" << std::endl;

  // Write field(s) to a text file
  std::string fileName = "solution" + fieldID + ".txt";
  FILE* f = fopen(fileName.c_str(), "w");

  for(int iNode = 0; iNode < nNod; ++iNode){
    int iDOF = n->getVertexNumber(iNode);
    output << sol[iDOF] << std::endl;
    fprintf(f, "%+-16.16e\n", sol[iDOF]);
  }

  fclose(f);
}

feExporterVTK::feExporterVTK(std::string vtkFile, feMesh *mesh, feSolution *sol, feMetaNumber *metaNumber,
  const std::vector<feSpace*> &space) : feExporter(mesh, sol, metaNumber, space)
{
  std::filebuf fb;
  if(fb.open(vtkFile,std::ios::out)){
    std::ostream output(&fb);
    writeHeader(output);

    // For now, we only print spaces of maximal dimension (non boundary)
    std::vector<feSpace*> spacesToExport;
    std::set<std::string> cncToExport;
    for(feSpace *fS : space){
      if(fS->getDim() == mesh->getDim()){
        spacesToExport.push_back(fS);
        cncToExport.insert(fS->getCncGeoID());
      }
    }

    // For now we only print one domain connectivity, assuming all fields are defined on the same connectivity.
    // To fix this, we have to join the connectivities : check which nodes are shared, and use a global numbering
    if(cncToExport.size() > 1){
      printf("In feExporterVTK : Warning - Multiple domains visualization is not ready yet.\n");
      printf("Only the first domain will be exported.\n");
    }

    // Grab the connectivity from any matching fespace
    feCncGeo *cnc;
    for(feSpace *fS : spacesToExport){
      if(fS->getCncGeoID() == *cncToExport.begin()){
        cnc = fS->getCncGeo();
        break;
      }
    }
    
    // Write nodes and elements
    writeNodes(output, cnc);
    writeElementsConnectivity(output, cnc);
    output << "POINT_DATA " << cnc->getNbNodes() << std::endl;

    // Write the field associated to each fespace in spacesToExport
    for(feSpace *fS : spacesToExport){
      if(cnc->getForme() == "TriP1"){
        writeField(output, cnc, fS->getFieldID());
      } else{
        printf("In feExporterVTK : So far only P1 triangles can be exported to Paraview...\n");
      }
    }

    fb.close();
  }
}
