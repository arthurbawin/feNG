#include "feNG.h"
#include "feMesh.h"
#include "feNumber.h"
#include "feTriangle.h"
#include "feSolution.h"

#include "SBoundingBox3d.h"

#include <iostream>
#include <fstream>

int feMesh::getCncGeoTag(std::string const &cncGeoID)
{
#ifdef DEBUG
  if(_cncGeoMap.find(cncGeoID) == _cncGeoMap.end()) return -1;
#endif
return _cncGeoMap[cncGeoID];
}

feCncGeo *feMesh::getCncGeoByName(std::string const &cncGeoID)
{
#ifdef DEBUG
  if(_cncGeoMap.find(cncGeoID) == _cncGeoMap.end()) return nullptr;
#endif
return _cncGeo[_cncGeoMap[cncGeoID]];
}

feCncGeo *feMesh::getCncGeoByTag(int cncGeoTag)
{
#ifdef DEBUG
  if(cncGeoTag < 0)
    if((unsigned)cncGeoTag > _cncGeo.size()) return nullptr;
#endif  
return _cncGeo[cncGeoTag];
}

// Returns a vector containing the physical coordinates of the nodes on the
// element with LOCAL number "numElem" on the connectivity named "cncGeoID" :
// coord = [x1 y1 z1 x2 y2 z2 ... xn yn zn]
void feMesh::getCoord(std::string const &cncGeoID, int numElem, std::vector<double> &geoCoord)
{
  feCncGeo *cnc = getCncGeoByName(cncGeoID);
  int nNodePerElem = cnc->getNbNodePerElem();
  geoCoord.resize(nNodePerElem * 3); // _dim = 3 pour les coordonnees
  for(int i = 0; i < nNodePerElem; ++i) {
    geoCoord[3 * i + 0] = _vertices[cnc->getNodeConnectivity(numElem, i)].x();
    geoCoord[3 * i + 1] = _vertices[cnc->getNodeConnectivity(numElem, i)].y();
    geoCoord[3 * i + 2] = _vertices[cnc->getNodeConnectivity(numElem, i)].z();
  }
}

// Returns a vector containing the physical coordinates of the nodes on the
// element with LOCAL number "numElem" on the connectivity numbered "cncGeoTag" :
// coord = [x1 y1 z1 x2 y2 z2 ... xn yn zn]
void feMesh::getCoord(int cncGeoTag, int numElem, std::vector<double> &geoCoord)
{
  feCncGeo *cnc = getCncGeoByTag(cncGeoTag);
  int nNodePerElem = cnc->getNbNodePerElem();
  geoCoord.resize(nNodePerElem * 3);
  for(int i = 0; i < nNodePerElem; ++i) {
    geoCoord[3 * i + 0] = _vertices[cnc->getNodeConnectivity(numElem, i)].x();
    geoCoord[3 * i + 1] = _vertices[cnc->getNodeConnectivity(numElem, i)].y();
    geoCoord[3 * i + 2] = _vertices[cnc->getNodeConnectivity(numElem, i)].z();
  }
}

int feMesh::getNbNodePerElem(std::string const &cncGeoID)
{
  return getCncGeoByName(cncGeoID)->getNbNodePerElem();
}

int feMesh::getNbNodePerElem(int cncGeoTag)
{
  return getCncGeoByTag(cncGeoTag)->getNbNodePerElem();
}

int feMesh::getVertex(std::string const &cncGeoID, int numElem, int numVertex)
{
  return getCncGeoByName(cncGeoID)->getNodeConnectivity(numElem, numVertex);
}

int feMesh::getVertex(int cncGeoTag, int numElem, int numVertex)
{
  return getCncGeoByTag(cncGeoTag)->getNodeConnectivity(numElem, numVertex);
}

int feMesh::getElement(std::string const &cncGeoID, int numElem)
{
  return getCncGeoByName(cncGeoID)->getElementConnectivity(numElem);
}

int feMesh::getElement(int cncGeoTag, int numElem)
{
  return getCncGeoByTag(cncGeoTag)->getElementConnectivity(numElem);
}

int feMesh::getEdge(std::string const &cncGeoID, int numElem, int numEdge)
{
  return getCncGeoByName(cncGeoID)->getEdgeConnectivity(numElem, numEdge);
}

int feMesh::getEdge(int cncGeoTag, int numElem, int numEdge)
{
  return getCncGeoByTag(cncGeoTag)->getEdgeConnectivity(numElem, numEdge);
}

int feMesh::getNbElm(std::string const &cncGeoID) { return getCncGeoByName(cncGeoID)->getNbElm(); }

int feMesh::getNbElm(int cncGeoTag) { return getCncGeoByTag(cncGeoTag)->getNbElm(); }

feSpace *feMesh::getGeometricSpace(std::string const &cncGeoID)
{
  return getCncGeoByName(cncGeoID)->getFeSpace();
}

feSpace *feMesh::getGeometricSpace(int cncGeoTag)
{
  return getCncGeoByTag(cncGeoTag)->getFeSpace();
}

void feMesh::printInfo(bool printConnectivities)
{
  std::cout << "Nombre total d'elements : " << _nTotalElm << std::endl;
  std::cout << "Nombre de connectivites : " << _nCncGeo << std::endl;
  for(feCncGeo *cnc : _cncGeo) {
    std::cout << "CncGeo " << cnc->getID() << " : " << cnc->getForme() << " - " << cnc->getNbElm()
              << " elements" << std::endl;
    if(printConnectivities) {
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
    _nNodBoundary(1)
{
  // Sommets
  _coord.resize(_nNod);
  for(int i = 0; i < _nNod; ++i) _coord[i] = xA + i * (xB - xA) / (_nNod - 1);

  _vertices.resize(_nNod);
  for(int i = 0; i < _nNod; ++i) _vertices[i] = Vertex(xA + i * (xB - xA) / (_nNod - 1), 0., 0., i);

  // Elements 1D
  int dimDomain = 1;
  std::vector<int> connecNodeDomain(_nElmDomain * _nNodDomain, 0);
  for(int i = 0; i < _nElmDomain; ++i) {
    connecNodeDomain[_nNodDomain * i + 0] = i;
    connecNodeDomain[_nNodDomain * i + 1] = i + 1;
  }
  // Create edges
  int nEdges = 1; // Numbered starting at 1 to match the 2D numbering
  std::vector<int> connecEdgeDomain(_nElmDomain, 0);
  for(int i = 0; i < _nElmDomain; ++i, ++nEdges) {
    Vertex *v0 = &_vertices[i];
    Vertex *v1 = &_vertices[i + 1];
    Edge e(v0, v1, nEdges, 0);
    _edges.insert(e);
    connecEdgeDomain[i] = nEdges; // Trivial edge connectivity but starting at 1
  }
  _nEdg = _nElmDomain;

  int nCncGeo = 0;
  feCncGeo *geoDom =
    new feCncGeo(nCncGeo, dimDomain, _nNodDomain, _nElmDomain, 1, _domID, "Lg",
                 new feSpace1DP1("xyz"), connecNodeDomain, std::vector<int>(), connecEdgeDomain);
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

feMesh1DP1::~feMesh1DP1()
{
  for(feCncGeo *cnc : _cncGeo) {
    delete cnc->getFeSpace();
    delete cnc;
  }
}

feMesh0DP0::feMesh0DP0(double xA, int nElm, std::string domID)
  : feMesh(1, 0, 1, "0D"), _nElm(nElm), _xA(xA), _domID(domID), _nElmDomain(nElm), _nElmBoundary(1),
    _nNodDomain(1), _nNodBoundary(1)
{
  // Sommets
  _coord.resize(_nNod);
  for(int i = 0; i < _nNod; ++i) {
    _coord[i] = xA;
  }

  _vertices.resize(_nNod);
  for(int i = 0; i < _nNod; ++i) {
    _vertices[i] = Vertex(0, 0., 0., 0);
  }

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

feMesh0DP0::~feMesh0DP0()
{
  for(feCncGeo *cnc : _cncGeo) {
    delete cnc->getFeSpace();
    delete cnc;
  }
}

static bool rtreeCallback(int id, void *ctx)
{
  std::vector<int> *vec = reinterpret_cast<std::vector<int> *>(ctx);
  vec->push_back(id);
  return true;
}

feMesh2DP1::feMesh2DP1(std::string meshName, bool curved, bool reversed,
                       mapType physicalEntitiesDescription)
  : feMesh()
{
  feStatus s = readGmsh(meshName, curved, reversed, physicalEntitiesDescription);
  if(s != FE_STATUS_OK) {
    feInfo("Error in readGmsh - mesh not finalized.\n");
    std::exit(1);
  }
}

feMesh2DP1::~feMesh2DP1()
{
  for(Triangle *t : _elements) {
    delete t;
  }
  for(feCncGeo *cnc : _cncGeo) {
    delete cnc->getFeSpace();
    delete cnc;
  }
}

/* Locates the vertex with coordinates x in the mesh using an RTree.
   The search is performed in elements of the highest dimension only.
   The element number is assigned to iElm and the reference coordinates
   are assigned in u. */
bool feMesh2DP1::locateVertex(std::vector<double> &x, int &iElm, std::vector<double> &u, double tol)
{
  double r[3];
  double min[3] = {x[0] - tol, x[1] - tol, x[2] - tol};
  double max[3] = {x[0] + tol, x[1] + tol, x[2] + tol};
  std::vector<int> candidates;
  _rtree.Search(min, max, rtreeCallback, &candidates);
  bool isFound = false;
  for(int val : candidates) {
    Triangle *t = _elements[val];
    t->xyz2uvw(x.data(), r);
    if(t->isInside(r[0], r[1], r[2])) {
      iElm = val;
      u[0] = r[0];
      u[1] = r[1];
      u[2] = r[2];
      isFound = true;
    }
  }
  return isFound;
}

/* Transfers the solution(s) associated to the current mesh to another mesh.
   The solution(s) are stored in solutionContainer and are transferred in place,
   hence the old solutions are lost after transfer. */
void feMesh2DP1::transfer(feMesh2DP1 *otherMesh, feMetaNumber *myMN, feMetaNumber *otherMN,
                          feSolutionContainer *solutionContainer,
                          const std::vector<feSpace *> &mySpaces,
                          const std::vector<feSpace *> &mySpacesEssBC,
                          const std::vector<feSpace *> &otherSpaces)
{
  // A temporary solutionContainer in which the interpolated values are stored
  feSolutionContainer *scTmp = new feSolutionContainer(*solutionContainer);
  int nSol = scTmp->getNbSol();
  int nDOF = otherMN->getNbDOFs();
  // Resize to the new number of DOFs
  for(int iSol = 0; iSol < nSol; ++iSol) {
    scTmp->_sol[iSol].resize(nDOF);
  }

  // The rtree structure in which the elements of the old mesh are stored
  RTree<int, double, 3> rtree;

  // Add domain (not boundary) elements to the rtree
  for(int j = 0; j < _elements.size(); ++j) {
    Triangle *t = _elements[j];
    SBoundingBox3d bbox;
    // printf("element %2d : %2d - %2d - %2d\n", t->getTag(), t->getVertex(0)->getTag(),
    // t->getVertex(1)->getTag(), t->getVertex(2)->getTag());
    for(int i = 0; i < 3; ++i) {
      Vertex *v = t->getVertex(i);
      SPoint3 pt(v->x(), v->y(), v->z());
      bbox += pt;
    }
    rtree.Insert((double *)(bbox.min()), (double *)(bbox.max()), j);
  }

  double tol = 1e-8;

  // Check which spaces need to be interpolated
  std::vector<std::pair<std::string, std::string> > fieldsToInterpolate;
  for(feSpace *fS1 : mySpaces) {
    bool isBC = false;
    if(fS1->getDim() < this->_dim) {
      // printf("Not interpolating boundary field of dimension %d\n", fS1->getDim());
      continue;
    }
    for(feSpace *feEss : mySpacesEssBC) {
      if(fS1->getFieldID() == feEss->getFieldID() && fS1->getCncGeoID() == feEss->getCncGeoID()) {
        isBC = true;
        // std::cout << "Not interpolating essential BC field " << fS1->getFieldID()
        //           << " on connectivity " << fS1->getCncGeoID() << std::endl;
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

          std::vector<feInt> adr1(fS1->getNbFunctions());
          std::vector<feInt> adr2(fS2->getNbFunctions());
          std::vector<double> geoCoord;
          for(int iElm = 0; iElm < nElm; ++iElm) {
            fS2->initializeAddressingVector(number2, iElm, adr2);
            otherMesh->getCoord(cncGeoTag, iElm, geoCoord);
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
              _rtree.Search(min, max, rtreeCallback, &candidates);

              for(int val : candidates) {
                Triangle *t = _elements[val];
                // printf("element %2d : %2d - %2d - %2d", t->getTag(), t->getVertex(0)->getTag(),
                // t->getVertex(1)->getTag(), t->getVertex(2)->getTag());
                double r2[3];
                t->xyz2uvw(x.data(), r2);
                bool isInside = t->isInside(r2[0], r2[1], r2[2]);
                if(isInside) {
                  fS1->initializeAddressingVector(
                    number1, t->getTag(),
                    adr1); // ATTENTION : I give the triangle tag as element number.
                           // Only works on triangles, not on 1D boundary elements.
                  for(int iSol = 0; iSol < nSol; ++iSol) {
                    std::vector<double> sol1(adr1.size());
                    std::vector<double> &solVec1 = solutionContainer->getSolution(iSol);
                    for(size_t i = 0; i < adr1.size(); ++i) {
                      sol1[i] = solVec1[adr1[i]];
                    }

                    double solInt = fS1->interpolateField(sol1, r2);
                    scTmp->_sol[iSol][adr2[j]] = solInt;
                    // TODO : Interpoler le _fresidual du solutionContainer
                  }
                }
              }
            }
          } // for iElm
        } else { // if fields match
          // std::cout << "NOT INTERPOLATING : MISMATCH between " << fS1->getFieldID() << " - "
          //           << fS1->getCncGeoID() << " and " << fS2->getFieldID() << " - "
          //           << fS2->getCncGeoID() << std::endl;
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

  // for(int iDOF = 0; iDOF < 20; ++iDOF) {
  //   printf("%5d : ", iDOF);
  //   for(int i = 0; i < solutionContainer->getNbSol(); ++i) {
  //     printf("%+-10.10e \t", solutionContainer->_sol[i][iDOF]);
  //   }
  //   printf("\n");
  // }

  delete scTmp;
}

// std::vector<vector<int>> feMesh2DP1::getListeElmPerNode()
// {
//   std::vector<vector<int>> ListeElmPerNode(_nNod,vector<int> (max(NbElmPerNode),0));

//   for

//   return(ListeElmPerNode);
// }
