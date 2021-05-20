#include "feMesh.h"

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
    printf("Connectivité des arêtes :\n");
    for(int i = 0; i < cnc->getNbElm(); ++i){
      for(int j = 0; j < cnc->getNbEdgePerElem(); ++j){
        std::cout<< cnc->getEdgeConnectivity(i,j)<<" - ";
      }
      std::cout<<std::endl;
    }
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
  std::vector<int> connecDomain(_nElmDomain*_nNodDomain, 0);
  for(int i = 0; i < _nElmDomain; ++i){
    connecDomain[_nNodDomain*i+0] = i;
    connecDomain[_nNodDomain*i+1] = i+1;
  }

  int nCncGeo = 0;
  feCncGeo *geoDom = new feCncGeo(_nNodDomain, _nElmDomain, 0, _domID, "Lg", new feSpace1DP1("xyz"), connecDomain);
  _cncGeo.push_back(geoDom);
  _cncGeoMap[_domID] = nCncGeo;
  geoDom->getFeSpace()->setCncGeoTag(nCncGeo++);

  // Elements 0D
  std::vector<int> connecBoundaryA(_nElmBoundary*_nNodBoundary, 0);
  std::vector<int> connecBoundaryB(_nElmBoundary*_nNodBoundary, 0);
  connecBoundaryA[0] = 0;
  connecBoundaryB[0] = _nNod-1;
  feCncGeo *geoBndA = new feCncGeo(_nNodBoundary, _nElmBoundary, 0, _bndA_ID, "Pt", new feSpace1DP0("xyz"), connecBoundaryA);
  feCncGeo *geoBndB = new feCncGeo(_nNodBoundary, _nElmBoundary, 0, _bndB_ID, "Pt", new feSpace1DP0("xyz"), connecBoundaryB);
  _cncGeo.push_back(geoBndA);
  _cncGeoMap[_bndA_ID] = nCncGeo;
  geoBndA->getFeSpace()->setCncGeoTag(nCncGeo++);
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
};

feMesh1DP1::~feMesh1DP1(){
  for(feCncGeo *cnc : _cncGeo){
    delete cnc->getFeSpace();
    delete cnc;
  }
}

feMesh2DP1::feMesh2DP1(std::string meshName, bool curved) : feMesh(){

  _ID = "myBeautifulMesh";
  // _dim = 2; // Should be determined by readGmsh

  if(readGmsh(meshName, curved)){
    printf("In feMesh2DP1::feMesh2DP1 : Error in readGmsh - mesh not finalized.\n");
  }

  // _nCncGeo = _cncGeo.size();

  // for(auto e : _edges)
  //   std::cout<<"Edge "<<e.getTag()<<" : "<<e.getTag(0)<<" - "<<e.getTag(1)<<std::endl;

  // std::cout<<_edges.size()<<std::endl;
  
}

feMesh2DP1::~feMesh2DP1(){
  for(feCncGeo *cnc : _cncGeo){
    delete cnc->getFeSpace();
    delete cnc;
  }
}