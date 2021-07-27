#include "feMesh.h"
#include "feNumber.h"
#include "feTriangle.h"
#include "feSolution.h"

#include "SBoundingBox3d.h"
#include "rtree.h"

#include <iostream>
#include <fstream>

int feMesh::getCncGeoTag(std::string cncGeoID) {
  if(_cncGeoMap.find(cncGeoID) != _cncGeoMap.end()) return _cncGeoMap[cncGeoID];
  return -1;
}

feCncGeo *feMesh::getCncGeoByName(std::string cncGeoID) {
  if(_cncGeoMap.find(cncGeoID) != _cncGeoMap.end()) return _cncGeo[_cncGeoMap[cncGeoID]];
  return nullptr;
}

feCncGeo *feMesh::getCncGeoByTag(int cncGeoTag) {
  if(cncGeoTag >= 0)
    if((unsigned)cncGeoTag < _cncGeo.size()) return _cncGeo[cncGeoTag];
  return nullptr;
}

// numElem est le numero local de l'element dans la connectivité géométrique
std::vector<double> feMesh::getCoord(std::string cncGeoID, int numElem) {
  feCncGeo *cnc = getCncGeoByName(cncGeoID);
  int nNodePerElem = cnc->getNbNodePerElem();
  std::vector<double> coord(nNodePerElem * 3); // _dim = 3 pour les coordonnees
  for(int i = 0; i < nNodePerElem; ++i) {
    coord[3 * i + 0] = _vertices[cnc->getNodeConnectivity(numElem, i)].x();
    coord[3 * i + 1] = _vertices[cnc->getNodeConnectivity(numElem, i)].y();
    coord[3 * i + 2] = _vertices[cnc->getNodeConnectivity(numElem, i)].z();
  }
  return coord;
}

std::vector<double> feMesh::getCoord(int cncGeoTag, int numElem) {
  feCncGeo *cnc = getCncGeoByTag(cncGeoTag);
  int nNodePerElem = cnc->getNbNodePerElem();
  std::vector<double> coord(nNodePerElem * 3);
  for(int i = 0; i < nNodePerElem; ++i) {
    coord[3 * i + 0] = _vertices[cnc->getNodeConnectivity(numElem, i)].x();
    coord[3 * i + 1] = _vertices[cnc->getNodeConnectivity(numElem, i)].y();
    coord[3 * i + 2] = _vertices[cnc->getNodeConnectivity(numElem, i)].z();
  }
  return coord;
}

int feMesh::getNbNodePerElem(std::string cncGeoID) {
  return getCncGeoByName(cncGeoID)->getNbNodePerElem();
}

int feMesh::getNbNodePerElem(int cncGeoTag) {
  return getCncGeoByTag(cncGeoTag)->getNbNodePerElem();
}

int feMesh::getVertex(std::string cncGeoID, int numElem, int numVertex) {
  return getCncGeoByName(cncGeoID)->getNodeConnectivity(numElem, numVertex);
}

int feMesh::getVertex(int cncGeoTag, int numElem, int numVertex) {
  return getCncGeoByTag(cncGeoTag)->getNodeConnectivity(numElem, numVertex);
}

int feMesh::getElement(std::string cncGeoID, int numElem) {
  return getCncGeoByName(cncGeoID)->getElementConnectivity(numElem);
}

int feMesh::getElement(int cncGeoTag, int numElem) {
  return getCncGeoByTag(cncGeoTag)->getElementConnectivity(numElem);
}

int feMesh::getEdge(std::string cncGeoID, int numElem, int numEdge) {
  return getCncGeoByName(cncGeoID)->getEdgeConnectivity(numElem, numEdge);
}

int feMesh::getEdge(int cncGeoTag, int numElem, int numEdge) {
  return getCncGeoByTag(cncGeoTag)->getEdgeConnectivity(numElem, numEdge);
}

int feMesh::getNbElm(std::string cncGeoID) { return getCncGeoByName(cncGeoID)->getNbElm(); }

int feMesh::getNbElm(int cncGeoTag) { return getCncGeoByTag(cncGeoTag)->getNbElm(); }

feSpace *feMesh::getGeometricSpace(std::string cncGeoID) {
  return getCncGeoByName(cncGeoID)->getFeSpace();
}

feSpace *feMesh::getGeometricSpace(int cncGeoTag) {
  return getCncGeoByTag(cncGeoTag)->getFeSpace();
}

void feMesh::printInfo() {
  std::cout << "Nombre total d'elements : " << _nTotalElm << std::endl;
  std::cout << "Nombre de connectivites : " << _nCncGeo << std::endl;
  for(feCncGeo *cnc : _cncGeo) {
    std::cout << "CncGeo " << cnc->getID() << std::endl;
    printf("Connectivité des noeuds :\n");
    for(int i = 0; i < cnc->getNbElm(); ++i) {
      for(int j = 0; j < cnc->getNbNodePerElem(); ++j) {
        std::cout << cnc->getNodeConnectivity(i, j) << " ";
      }
      std::cout << std::endl;
    }
    printf("Connectivité des éléments :\n");
    for(int i = 0; i < cnc->getNbElm(); ++i) {
      std::cout << cnc->getElementConnectivity(i) << std::endl;
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

feMesh1DP1::feMesh1DP1(double xA, double xB, int nElm, std::string bndA_ID, std::string bndB_ID,
                       std::string domID)
  : feMesh(nElm + 1, 1, 3, "1D"), _nElm(nElm), _xA(xA), _xB(xB), _bndA_ID(bndA_ID),
    _bndB_ID(bndB_ID), _domID(domID), _nElmDomain(nElm), _nElmBoundary(1), _nNodDomain(2),
    _nNodBoundary(1) {
  // Sommets
  _coord.resize(_nNod);
  for(int i = 0; i < _nNod; ++i) _coord[i] = xA + i * (xB - xA) / (_nNod - 1);

  _vertices.resize(_nNod);
  for(int i = 0; i < _nNod; ++i) _vertices[i] = Vertex(xA + i * (xB - xA) / (_nNod - 1), 0., 0., i);

  // Elements 1D
  int dimDomain = 1;
  std::vector<int> connecDomain(_nElmDomain * _nNodDomain, 0);
  for(int i = 0; i < _nElmDomain; ++i) {
    connecDomain[_nNodDomain * i + 0] = i;
    connecDomain[_nNodDomain * i + 1] = i + 1;
  }
  int nCncGeo = 0;
  feCncGeo *geoDom = new feCncGeo(nCncGeo, dimDomain, _nNodDomain, _nElmDomain, 0, _domID, "Lg",
                                  new feSpace1DP1("xyz"), connecDomain);
  _cncGeo.push_back(geoDom);
  _cncGeoMap[_domID] = nCncGeo;
  geoDom->getFeSpace()->setCncGeoTag(nCncGeo++);

  // Elements 0D
  int dimBoundary = 0;
  std::vector<int> connecBoundaryA(_nElmBoundary * _nNodBoundary, 0);
  std::vector<int> connecBoundaryB(_nElmBoundary * _nNodBoundary, 0);
  connecBoundaryA[0] = 0;
  connecBoundaryB[0] = _nNod - 1;

  feCncGeo *geoBndA = new feCncGeo(nCncGeo, dimBoundary, _nNodBoundary, _nElmBoundary, 0, _bndA_ID,
                                   "Pt", new feSpace1DP0("xyz"), connecBoundaryA);
  _cncGeo.push_back(geoBndA);
  _cncGeoMap[_bndA_ID] = nCncGeo;
  geoBndA->getFeSpace()->setCncGeoTag(nCncGeo++);

  feCncGeo *geoBndB = new feCncGeo(nCncGeo, dimBoundary, _nNodBoundary, _nElmBoundary, 0, _bndB_ID,
                                   "Pt", new feSpace1DP0("xyz"), connecBoundaryB);
  _cncGeo.push_back(geoBndB);
  _cncGeoMap[_bndB_ID] = nCncGeo;
  geoBndB->getFeSpace()->setCncGeoTag(nCncGeo++);

  // Fonction COMPLETER : connectivite globale des elements
  int numElmGlo = 0;
  for(feCncGeo *cnc : _cncGeo) {
    for(int i = 0; i < cnc->getNbElm(); ++i) cnc->setElementConnectivity(i, numElmGlo++);
  }
  _nTotalElm = numElmGlo;

  // Assign pointer to the mesh to cnc and their fespace :
  for(feCncGeo *cnc : _cncGeo) {
    cnc->setMeshPtr(this);
    cnc->getFeSpace()->setMeshPtr(this);
  }
}

feMesh1DP1::~feMesh1DP1() {
  for(feCncGeo *cnc : _cncGeo) {
    delete cnc->getFeSpace();
    delete cnc;
  }
}

feMesh0DP0::feMesh0DP0(double xA, int nElm, std::string domID)
  : feMesh(1, 0, 1, "0D"), _nElm(nElm), _xA(xA), _domID(domID), _nElmDomain(nElm), _nElmBoundary(1),
    _nNodDomain(1), _nNodBoundary(1) {
  // Sommets
  _coord.resize(_nNod);
  for(int i = 0; i < _nNod; ++i) { _coord[i] = xA; }

  _vertices.resize(_nNod);
  for(int i = 0; i < _nNod; ++i) { _vertices[i] = Vertex(0, 0., 0., 0); }

  // Element 0D
  int dimDomain = 0;
  std::vector<int> connecDomain(1, 0);

  int nCncGeo = 0;
  feCncGeo *geoDom = new feCncGeo(nCncGeo, dimDomain, _nNodDomain, _nElmDomain, 0, _domID, "Pt",
                                  new feSpace1DP0("xyz"), connecDomain);
  _cncGeo.push_back(geoDom);
  _cncGeoMap[_domID] = nCncGeo;
  geoDom->getFeSpace()->setCncGeoTag(nCncGeo++);

  // Fonction COMPLETER : connectivite globale des elements
  int numElmGlo = 0;
  for(feCncGeo *cnc : _cncGeo) {
    for(int i = 0; i < cnc->getNbElm(); ++i) cnc->setElementConnectivity(i, numElmGlo++);
  }
  _nTotalElm = numElmGlo;

  // Assign pointer to the mesh to cnc and their fespace :
  for(feCncGeo *cnc : _cncGeo) {
    cnc->setMeshPtr(this);
    cnc->getFeSpace()->setMeshPtr(this);
  }
}

feMesh0DP0::~feMesh0DP0() {
  for(feCncGeo *cnc : _cncGeo) {
    delete cnc->getFeSpace();
    delete cnc;
  }
}

feMesh2DP1::feMesh2DP1(std::string meshName, bool curved, mapType physicalEntitiesDescription)
  : feMesh() {
  // Check if mesh file exists
  std::ifstream f(meshName.c_str());
  if(!f.good()) {
    printf("In feMesh2DP1::feMesh2DP1 : Error - Mesh file does not exist.\n");
  } else {
    _ID = "myBeautifulMesh";

    if(readGmsh(meshName, curved, physicalEntitiesDescription)) {
      printf("In feMesh2DP1::feMesh2DP1 : Error in readGmsh - mesh not finalized.\n");
    }

    // Assign a pointer to this mesh to each of its geometric connectivities and their fespace
    for(feCncGeo *cnc : _cncGeo) {
      cnc->setMeshPtr(this);
      cnc->getFeSpace()->setMeshPtr(this);
    }
  }
}

feMesh2DP1::~feMesh2DP1() {
  for(Triangle *t : _elements) { delete t; }
  for(feCncGeo *cnc : _cncGeo) {
    delete cnc->getFeSpace();
    delete cnc;
  }
}

static bool rtreeCallback(int id, void *ctx) {
  std::vector<int> *vec = reinterpret_cast<std::vector<int> *>(ctx);
  vec->push_back(id);
  return true;
}

/* Transfers the solution(s) associated to the current mesh to another mesh.
   The solution(s) are stored in solutionContainer and are transferred in place,
   hence the old solutions are lost after transfer. */
void feMesh2DP1::transfer(feMesh2DP1 *otherMesh, feMetaNumber *myMN, feMetaNumber *otherMN,
                          feSolutionContainer *solutionContainer,
                          const std::vector<feSpace *> &mySpaces,
                          const std::vector<feSpace *> &mySpacesEssBC,
                          const std::vector<feSpace *> &otherSpaces) {
  // A temporary solutionContainer in which the interpolated values are stored
  feSolutionContainer *scTmp = new feSolutionContainer(*solutionContainer);
  int nSol = scTmp->getNbSol();
  int nDOF = otherMN->getNbDOFs();
  // Resize to the new number of DOFs
  for(int iSol = 0; iSol < nSol; ++iSol) { scTmp->_sol[iSol].resize(nDOF); }

  // The rtree structure in which the elements of the old mesh are stored
  RTree<int, double, 3> rtree;

  // Add domain (not boundary) elements to the rtree
  for(int i = 0; i < _elements.size(); ++i) {
    Triangle *t = _elements[i];
    SBoundingBox3d bbox;
    // printf("element %2d : %2d - %2d - %2d\n", t->getTag(), t->getVertex(0)->getTag(),
    // t->getVertex(1)->getTag(), t->getVertex(2)->getTag());
    for(int i = 0; i < 3; ++i) {
      Vertex *v = t->getVertex(i);
      SPoint3 pt(v->x(), v->y(), v->z());
      bbox += pt;
    }
    rtree.Insert((double *)(bbox.min()), (double *)(bbox.max()), i);
  }

  double tol = 1e-8;

  // Check which spaces need to be interpolated
  std::vector<std::pair<std::string, std::string> > fieldsToInterpolate;
  for(feSpace *fS1 : mySpaces) {
    bool isBC = false;
    if(fS1->getDim() < this->_dim) {
      printf("Not interpolating boundary field of dimension %d\n", fS1->getDim());
      continue;
    }
    for(feSpace *feEss : mySpacesEssBC) {
      if(fS1->getFieldID() == feEss->getFieldID() && fS1->getCncGeoID() == feEss->getCncGeoID()) {
        isBC = true;
        std::cout << "Not interpolating essential BC field " << fS1->getFieldID()
                  << " on connectivity " << fS1->getCncGeoID() << std::endl;
      }
    }
    if(!isBC) {
      for(feSpace *fS2 : otherSpaces) {
        if(fS1->getFieldID() == fS2->getFieldID() && fS1->getCncGeoID() == fS2->getCncGeoID()) {
          std::cout << "INTERPOLATING FIELD " << fS1->getFieldID() << " on connectivity "
                    << fS1->getCncGeoID() << std::endl;
          fieldsToInterpolate.emplace_back(fS1->getFieldID(), fS1->getCncGeoID());

          int nElm = fS2->getNbElm();
          int nDOFPerElem = fS2->getNbFunctions();
          int cncGeoTag = fS2->getCncGeoTag();
          const std::vector<double> &Lcoor = fS2->getLcoor();
          feNumber *number1 = myMN->getNumbering(fS1->getFieldID());
          feNumber *number2 = otherMN->getNumbering(fS2->getFieldID());
          feSpace *geoSpace2 = fS2->getCncGeo()->getFeSpace();

          for(int iElm = 0; iElm < nElm; ++iElm) {
            fS2->initializeAddressingVector(number2, iElm);
            std::vector<double> geoCoord = otherMesh->getCoord(cncGeoTag, iElm);
            // Loop over the DOFs of the element of the new mesh
            for(int j = 0; j < nDOFPerElem; ++j) {
              // Get coordinates
              std::vector<double> x(3, 0.);
              double r1[3] = {Lcoor[3 * j + 0], Lcoor[3 * j + 1], Lcoor[3 * j + 2]};
              geoSpace2->interpolateVectorField(geoCoord, r1, x);
              // Locate DOF in old mesh
              double min[3] = {x[0] - tol, x[1] - tol, x[2] - tol};
              double max[3] = {x[0] + tol, x[1] + tol, x[2] + tol};
              std::vector<int> candidates;
              rtree.Search(min, max, rtreeCallback, &candidates);

              for(int val : candidates) {
                Triangle *t = _elements[val];
                // printf("element %2d : %2d - %2d - %2d", t->getTag(), t->getVertex(0)->getTag(),
                // t->getVertex(1)->getTag(), t->getVertex(2)->getTag());
                double r2[3];
                t->xyz2uvw(x.data(), r2);
                bool isInside = t->isInside(r2[0], r2[1], r2[2]);
                if(isInside) {
                  fS1->initializeAddressingVector(
                    number1, t->getTag()); // ATTENTION : I give the triangle tag as element number.
                                           // Only works on triangles, not on 1D boundary elements.
                  for(int iSol = 0; iSol < nSol; ++iSol) {
                    fS1->initializeSolution(solutionContainer->getSolution(iSol));
                    double solInt = fS1->interpolateSolution(r2);
                    scTmp->_sol[iSol][fS2->getAddressingVectorAt(j)] = solInt;
                    // TODO : Interpoler le _fresidual du solutionContainer
                  }
                }
              }
            }
          } // for iElm
        } else { // if fields match
          std::cout << "NOT INTERPOLATING : MISMATCH between " << fS1->getFieldID() << " - "
                    << fS1->getCncGeoID() << " and " << fS2->getFieldID() << " - "
                    << fS2->getCncGeoID() << std::endl;
        }
      }
    }
  }

  // Resize the solutionContainer
  solutionContainer->setNbDOFs(nDOF);
  solutionContainer->_d.resize(nDOF);

  // Replace solutionContainer by scTmp
  for(int iSol = 0; iSol < nSol; ++iSol) {
    solutionContainer->_sol[iSol].resize(nDOF);
    solutionContainer->_fResidual[iSol].resize(nDOF);
    for(int iDOF = 0; iDOF < nDOF; ++iDOF)
      solutionContainer->_sol[iSol][iDOF] = scTmp->_sol[iSol][iDOF];
  }

  for(int iDOF = 0; iDOF < 20; ++iDOF) {
    printf("%5d : ", iDOF);
    for(int i = 0; i < solutionContainer->getNbSol(); ++i) {
      printf("%+-10.10e \t", solutionContainer->_sol[i][iDOF]);
    }
    printf("\n");
  }

  delete scTmp;
}