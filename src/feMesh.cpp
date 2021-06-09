#include "feMesh.h"
#include "feNumber.h"

#include "SBoundingBox3d.h"
#include "rtree.h"

#include <iostream>

int feMesh::getCncGeoTag(std::string cncGeoID){
  if(_cncGeoMap.find(cncGeoID) != _cncGeoMap.end())
    return _cncGeoMap[cncGeoID];
  return -1;
}

feCncGeo* feMesh::getCncGeoByName(std::string cncGeoID){
  if(_cncGeoMap.find(cncGeoID) != _cncGeoMap.end())
    return _cncGeo[_cncGeoMap[cncGeoID]];
  return nullptr;
}

feCncGeo* feMesh::getCncGeoByTag(int cncGeoTag){
  if(cncGeoTag >= 0)
    if((unsigned) cncGeoTag < _cncGeo.size())
      return _cncGeo[cncGeoTag];
  return nullptr;
}

// std::vector<double> feMesh::getCoord(std::string cncGeoID, int numElem){
//   feCncGeo* cnc = getCncGeoByName(cncGeoID);
//   int nNodePerElem = cnc->getNbNodePerElem();
//   std::vector<double> coord(nNodePerElem * _dim); // _dim = 3 pour les coordonnees
//   for(int i = 0; i < nNodePerElem; ++i){
//     for(int j = 0; j < _dim; ++j){
//       coord[_dim * i + j] = _coord[_dim * cnc->getNodeConnectivity(numElem, i) + j];
//       // coord[3 * i + 0] = _vertices[3 * cnc->getNodeConnectivity(numElem, i)].x();
//       // coord[3 * i + 1] = _vertices[3 * cnc->getNodeConnectivity(numElem, i)].y();
//       // coord[3 * i + 2] = _vertices[3 * cnc->getNodeConnectivity(numElem, i)].z();
//     }
//   }
//   return coord;
// }

// std::vector<double> feMesh::getCoord(int cncGeoTag, int numElem){
//   feCncGeo* cnc = getCncGeoByTag(cncGeoTag);
//   int nNodePerElem = cnc->getNbNodePerElem();
//   std::vector<double> coord(nNodePerElem * _dim);
//   for(int i = 0; i < nNodePerElem; ++i){
//     for(int j = 0; j < _dim; ++j){
//       coord[_dim * i + j] = _coord[_dim * cnc->getNodeConnectivity(numElem, i) + j];
//     }
//   }
//   return coord;
// }

// numElem est le numero local de l'element dans la connectivité géométrique
std::vector<double> feMesh::getCoord(std::string cncGeoID, int numElem){
  feCncGeo* cnc = getCncGeoByName(cncGeoID);
  int nNodePerElem = cnc->getNbNodePerElem();
  std::vector<double> coord(nNodePerElem * 3); // _dim = 3 pour les coordonnees
  for(int i = 0; i < nNodePerElem; ++i){
    coord[3 * i + 0] = _vertices[cnc->getNodeConnectivity(numElem, i)].x();
    coord[3 * i + 1] = _vertices[cnc->getNodeConnectivity(numElem, i)].y();
    coord[3 * i + 2] = _vertices[cnc->getNodeConnectivity(numElem, i)].z();
  }
  return coord;
}

std::vector<double> feMesh::getCoord(int cncGeoTag, int numElem){
  feCncGeo* cnc = getCncGeoByTag(cncGeoTag);
  int nNodePerElem = cnc->getNbNodePerElem();
  std::vector<double> coord(nNodePerElem * 3);
  for(int i = 0; i < nNodePerElem; ++i){
    coord[3 * i + 0] = _vertices[cnc->getNodeConnectivity(numElem, i)].x();
    coord[3 * i + 1] = _vertices[cnc->getNodeConnectivity(numElem, i)].y();
    coord[3 * i + 2] = _vertices[cnc->getNodeConnectivity(numElem, i)].z();
  }
  return coord;
}

int feMesh::getNbNodePerElem(std::string cncGeoID){
  return getCncGeoByName(cncGeoID)->getNbNodePerElem();
}

int feMesh::getNbNodePerElem(int cncGeoTag){
  return getCncGeoByTag(cncGeoTag)->getNbNodePerElem();
}


int feMesh::getVertex(std::string cncGeoID, int numElem, int numVertex){
  return getCncGeoByName(cncGeoID)->getNodeConnectivity(numElem, numVertex);
}

int feMesh::getVertex(int cncGeoTag, int numElem, int numVertex){
  return getCncGeoByTag(cncGeoTag)->getNodeConnectivity(numElem, numVertex);
}

int feMesh::getElement(std::string cncGeoID, int numElem){
  return getCncGeoByName(cncGeoID)->getGlobalConnectivity(numElem);
}

int feMesh::getElement(int cncGeoTag, int numElem){
  return getCncGeoByTag(cncGeoTag)->getGlobalConnectivity(numElem);
}

int feMesh::getEdge(std::string cncGeoID, int numElem, int numEdge){
  return getCncGeoByName(cncGeoID)->getEdgeConnectivity(numElem, numEdge);
}

int feMesh::getEdge(int cncGeoTag, int numElem, int numEdge){
  return getCncGeoByTag(cncGeoTag)->getEdgeConnectivity(numElem, numEdge);
}

int feMesh::getNbElm(std::string cncGeoID){
  return getCncGeoByName(cncGeoID)->getNbElm();
}

int feMesh::getNbElm(int cncGeoTag){
  return getCncGeoByTag(cncGeoTag)->getNbElm();
}

feSpace* feMesh::getGeometricSpace(std::string cncGeoID){
  return getCncGeoByName(cncGeoID)->getFeSpace();
}

feSpace* feMesh::getGeometricSpace(int cncGeoTag){
  return getCncGeoByTag(cncGeoTag)->getFeSpace();
}

void feMesh::printInfo(){
  std::cout<<"Nombre total d'elements : "<<_nTotalElm<<std::endl;
  std::cout<<"Nombre de connectivites : "<<_nCncGeo<<std::endl;
  for(feCncGeo *cnc : _cncGeo){
    std::cout<< "CncGeo "<<cnc->getID()<<std::endl;
    printf("Connectivité des noeuds :\n");
    for(int i = 0; i < cnc->getNbElm(); ++i){
      for(int j = 0; j < cnc->getNbNodePerElem(); ++j){
        std::cout<< cnc->getLocalConnectivity(i,j)<<" ";
      }
      std::cout<<std::endl;
    }
    printf("Connectivité des éléments :\n");
    for(int i = 0; i < cnc->getNbElm(); ++i){
      std::cout<< cnc->getGlobalConnectivity(i) << std::endl;
    }
    // printf("Connectivité des arêtes :\n");
    // for(int i = 0; i < cnc->getNbElm(); ++i){
    //   for(int j = 0; j < cnc->getNbEdgePerElem(); ++j){
    //     std::cout<< cnc->getEdgeConnectivity(i,j)<<" - ";
    //   }
    //   std::cout<<std::endl;
    // }
  }
}

feMesh1DP1::feMesh1DP1(double xA, double xB, int nElm, std::string bndA_ID, std::string bndB_ID, std::string domID)
  : feMesh(nElm+1, 1, 3, "1D"), _nElm(nElm), _xA(xA), _xB(xB) , _bndA_ID(bndA_ID), _bndB_ID(bndB_ID), _domID(domID),
  _nElmDomain(nElm), _nElmBoundary(1), _nNodDomain(2), _nNodBoundary(1)
{
  // Sommets
  _coord.resize(_nNod);
  for(int i = 0; i < _nNod; ++i)
    _coord[i] = xA + i * (xB-xA)/(_nNod-1);

  _vertices.resize(_nNod);
  for(int i = 0; i < _nNod; ++i)
    _vertices[i] = Vertex(xA + i * (xB-xA)/(_nNod-1), 0., 0., i);

  // Elements 1D
  int dimDomain = 1;
  std::vector<int> connecDomain(_nElmDomain*_nNodDomain, 0);
  for(int i = 0; i < _nElmDomain; ++i){
    connecDomain[_nNodDomain*i+0] = i;
    connecDomain[_nNodDomain*i+1] = i+1;
  }

  int nCncGeo = 0;
  feCncGeo *geoDom = new feCncGeo(nCncGeo, dimDomain, _nNodDomain, _nElmDomain, 0, _domID, "Lg", new feSpace1DP1("xyz"), connecDomain);
  _cncGeo.push_back(geoDom);
  _cncGeoMap[_domID] = nCncGeo;
  geoDom->getFeSpace()->setCncGeoTag(nCncGeo++);

  // Elements 0D
  int dimBoundary = 0;
  std::vector<int> connecBoundaryA(_nElmBoundary*_nNodBoundary, 0);
  std::vector<int> connecBoundaryB(_nElmBoundary*_nNodBoundary, 0);
  connecBoundaryA[0] = 0;
  connecBoundaryB[0] = _nNod-1;

  feCncGeo *geoBndA = new feCncGeo(nCncGeo, dimBoundary, _nNodBoundary, _nElmBoundary, 0, _bndA_ID, "Pt", new feSpace1DP0("xyz"), connecBoundaryA);
  _cncGeo.push_back(geoBndA);
  _cncGeoMap[_bndA_ID] = nCncGeo;
  geoBndA->getFeSpace()->setCncGeoTag(nCncGeo++);

  feCncGeo *geoBndB = new feCncGeo(nCncGeo, dimBoundary, _nNodBoundary, _nElmBoundary, 0, _bndB_ID, "Pt", new feSpace1DP0("xyz"), connecBoundaryB);
  _cncGeo.push_back(geoBndB);
  _cncGeoMap[_bndB_ID] = nCncGeo;
  geoBndB->getFeSpace()->setCncGeoTag(nCncGeo++);

  // Fonction COMPLETER : connectivite globale des elements
  int numElmGlo = 0;
  for(feCncGeo *cnc : _cncGeo){
    for(int i = 0; i < cnc->getNbElm(); ++i)
      cnc->setGlobalConnectivity(i,numElmGlo++);
  }
  _nTotalElm = numElmGlo;

  // Assign pointer to the mesh to cnc and their fespace :
  for(feCncGeo *cnc : _cncGeo){
    cnc->setMeshPtr(this);
    cnc->getFeSpace()->setMeshPtr(this);
  }
};

feMesh1DP1::~feMesh1DP1(){
  for(feCncGeo *cnc : _cncGeo){
    delete cnc->getFeSpace();
    delete cnc;
  }
}

feMesh2DP1::feMesh2DP1(std::string meshName, bool curved, mapType physicalEntitiesDescription) : feMesh(){

  _ID = "myBeautifulMesh";

  if(readGmsh(meshName, curved, physicalEntitiesDescription)){
    printf("In feMesh2DP1::feMesh2DP1 : Error in readGmsh - mesh not finalized.\n");
  }

  // Assign a pointer to this mesh to each of its geometric connectivities and their fespace
  for(feCncGeo *cnc : _cncGeo){
    cnc->setMeshPtr(this);
    cnc->getFeSpace()->setMeshPtr(this);
  }
}

feMesh2DP1::~feMesh2DP1(){
  for(feCncGeo *cnc : _cncGeo){
    delete cnc->getFeSpace();
    delete cnc;
  }
}

static bool rtreeCallback(int id, void *ctx){
  std::vector<int> *vec = reinterpret_cast<std::vector<int> *>(ctx);
  vec->push_back(id);
  return true;
}

void feMesh2DP1::transfer(const feMesh2DP1 *otherMesh, const feMetaNumber *myMN, const feMetaNumber *otherMN, 
  const feSolution *mySol, feSolution* otherSol, const std::vector<feSpace*> &mySpaces, const std::vector<feSpace*> &otherSpaces){

  RTree<int, double, 3>  rtree;

  // Add domain (not boundary) elements to an rtree
  for(feSpace *fS : mySpaces){
    if(fS->getDim() == this->_dim){

      feCncGeo *cnc = fS->getCncGeo();
      int nElm = cnc->getNbElm();
      int nNodePerElm = cnc->getNbNodePerElem(); // All nodes are used to create the bounding box, not just the vertices

      for(int iElm = 0; iElm < nElm; ++iElm){
        SBoundingBox3d bbox;
        printf("element %2d : %2d - %2d - %2d\n", iElm, cnc->getNodeConnectivity(iElm, 0), cnc->getNodeConnectivity(iElm, 1), cnc->getNodeConnectivity(iElm, 2));
        for(int iNode = 0; iNode < nNodePerElm; ++iNode){
          int node = cnc->getNodeConnectivity(iElm, iNode);
          SPoint3 pt(_vertices[node].x(), _vertices[node].y(), _vertices[node].z());
          bbox += pt;
          // std::cout<<"Point added :"<<std::endl;
          // std::cout<<bbox.min().x()<<std::endl;
          // std::cout<<bbox.min().y()<<std::endl;
          // std::cout<<bbox.min().z()<<std::endl;
          // std::cout<<bbox.max().x()<<std::endl;
          // std::cout<<bbox.max().y()<<std::endl;
          // std::cout<<bbox.max().z()<<std::endl;
          rtree.Insert((double*)(bbox.min()), (double*)(bbox.max()), iElm);
          // rtree.Insert(bbox.min().data(), bbox.max().data(), iElm);
        }
      }
    }
  }

  // Find point 0.5,0.5,0
  double tol = 1e-2;
  double min[3] = {0.5-tol, 0.5-tol, 0.0-tol};
  double max[3] = {0.5+tol, 0.5+tol, 0.0+tol};
  std::vector<int> candidates;
  rtree.Search(min, max, rtreeCallback, &candidates);

  for(auto val : candidates){
    std::cout<<"val = "<<val<<std::endl;
  }
  

  SBoundingBox3d foo;

// std::vector<uint64_t> candidates;
//   forestOptions->bndRTree->Search(min, max, rtreeCallback, &candidates);

  // Add bboxes of the surface mesh to rtree
  // HXTBbox bbox_triangle;
  // for(uint64_t i = 0; i < meshBnd->triangles.num; ++i){
  //   hxtBboxInit(&bbox_triangle);
  //   for(uint64_t j = 0; j < 3; ++j) {
  //     double coord[3];
  //     uint32_t node = meshBnd->triangles.node[3*i+j];
  //     for(uint32_t k = 0; k < 3; ++k){ coord[k] = meshBnd->vertices.coord[(size_t) 4*node+k]; }
  //     hxtBboxAddOne(&bbox_triangle, coord);
  //   }
  //   SBoundingBox3d cube_bbox(bbox_triangle.min[0], bbox_triangle.min[1], bbox_triangle.min[2],
  //                            bbox_triangle.max[0], bbox_triangle.max[1], bbox_triangle.max[2]);
  //   bndRTree.Insert((double*)(cube_bbox.min()), (double*)(cube_bbox.max()), i);
  // }

}

// void feNorm::computeL2Norm(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh){
//   double normL2 = 0.0, solInt, solRef, xInt, yInt, zInt, J, t = sol->getCurrentTime();
//   int nElm = _intSpace->getNbElm();

//   // std::cout<<"Computing norm of "<<_intSpace->getFieldID()<<" on "<<nElm<<" elements"<<std::endl;
//   for(int iElm = 0; iElm < nElm; ++iElm){
//     _intSpace->initializeAddressingVector(metaNumber->getNumbering(_intSpace->getFieldID()), iElm);
//     _intSpace->initializeSolution(sol);
//     geoCoord = mesh->getCoord(cncGeoTag, iElm);
//     for(int k = 0; k < _nQuad; ++k){
//       solInt = _intSpace->interpolateSolutionAtQuadNode(k);
//     }
//   }
//   norm = sqrt(normL2);
// }