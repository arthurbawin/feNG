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

void feExporterVTK::writeNodes(std::ostream& output, int iField){
  feNumber *n = _metaNumber->getNumbering(iField);
  int nNod = n->getNbNodes();
  output << "DATASET UNSTRUCTURED_GRID\n";
  output << "POINTS " << nNod << " double\n";
  Vertex *v;
  for(int i = 0; i < nNod; ++i){
    v = _mesh->getVertex(i);
    output << v->x() << " " << v->y() << " " << v->z() << std::endl;
  }
}

void feExporterVTK::writeElementsConnectivity(std::ostream& output, feCncGeo *cnc, int iField){
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
  for(int iElm = 0; iElm < nElm; ++iElm){
    output << cncToVTKmap[cnc->getForme()] << std::endl;
  }
}

void feExporterVTK::writeField(std::ostream& output, feCncGeo *cnc, int iField){
  int nElm = cnc->getNbElm();
  int nNodePerElem = cnc->getNbNodePerElem();
  std::vector<double>& sol = _sol->getSolutionReference();
  // Pas juste : c'est le nombre de noeuds du maillage au complet et pas de cnc...
  feNumber *n = _metaNumber->getNumbering(iField);
  int nNod = n->getNbNodes();
  // output << "POINT_DATA " << nNod << std::endl;
  output << "SCALARS " << _metaNumber->getFieldID(iField) << " double 1" << std::endl;
  output << "LOOKUP_TABLE default" << std::endl;
  for(int iNode = 0; iNode < nNod; ++iNode){
    int iDOF = n->getVertexNumber(iNode);
    output << sol[iDOF] << std::endl;
  }
}

feExporterVTK::feExporterVTK(std::string vtkFile, feMesh *mesh, feSolution *sol, feMetaNumber *metaNumber,
  const std::vector<feSpace*> &space) : feExporter(mesh, sol, metaNumber, space)
{
  std::filebuf fb;
  if(fb.open(vtkFile,std::ios::out)){
    std::ostream output(&fb);
    writeHeader(output);
    // Write nodes and connectivity : (this assumes all fields are on the same connectivity : fix this in the future)
    writeNodes(output, 0);
    writeElementsConnectivity(output, space[4]->getCncGeo(), 0);
    output << "POINT_DATA " << _metaNumber->getNumbering(0)->getNbNodes() << std::endl;
    // Write each field :
    for(int iField = 0; iField < metaNumber->getNbFields(); ++iField){
      std::string fieldID = metaNumber->getFieldID(iField);
      // writeNodes(output, iField);
      for(auto fS : space){
        if(fS->getFieldID() == fieldID){
          // std::cout<<"Field "<<fieldID<<" is defined on fespace "<<fS->getFieldID()<<" - "<<fS->getCncGeoID()<<std::endl;
          // TODO : C'est de la triche : je teste si la connectivité est TriP1
          // Il faudrait peut-être écrire juste les connectivités non-frontières
          // Peut-être ajouter un flag iBoundary
          if(fS->getCncGeo()->getForme() == "TriP1"){
            // writeElementsConnectivity(output, fS->getCncGeo(), iField);
            writeField(output, fS->getCncGeo(), iField);
          }
        }

      }
    }
    fb.close();
  }
}
