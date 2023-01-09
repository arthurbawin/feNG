#include "feNG.h"
#include "feMesh.h"
#include "feNumber.h"
#include "feTriangle.h"
#include "feSolution.h"

#include "SBoundingBox3d.h"

#include <iostream>
#include <fstream>

extern int FE_VERBOSE;

int feMesh::getCncGeoTag(std::string const &cncGeoID)
{
  auto it = _cncGeoMap.find(cncGeoID);
  return (it == _cncGeoMap.end()) ? -1 : it->second;
}

feCncGeo *feMesh::getCncGeoByName(std::string const &cncGeoID)
{
  auto it = _cncGeoMap.find(cncGeoID); 
  return (it == _cncGeoMap.end()) ? nullptr : _cncGeo[it->second];
}

feCncGeo *feMesh::getCncGeoByTag(int cncGeoTag)
{
  if(cncGeoTag < 0)
    if((unsigned) cncGeoTag > _cncGeo.size()) return nullptr;
  return _cncGeo[cncGeoTag];
}

// Returns a vector containing the physical coordinates of the nodes on the
// element with LOCAL number "numElem" on the connectivity named "cncGeoID" :
// coord = [x1 y1 z1 x2 y2 z2 ... xn yn zn]
void feMesh::getCoord(std::string const &cncGeoID, int numElem, std::vector<double> &geoCoord)
{
  feCncGeo *cnc = getCncGeoByName(cncGeoID);
  int nNodePerElem = cnc->getNbNodePerElem();
#ifdef FENG_DEBUG
  if(geoCoord.size() != 3 * nNodePerElem) {
    printf(" In feMesh::getCoord : Erreur - Wrong Size for vector geoCoord\n");
  }
#endif
  for(int i = 0; i < nNodePerElem; ++i) {
    Vertex V = _vertices[cnc->getNodeConnectivity(numElem, i)];
    geoCoord[3 * i + 0] = V.x();
    geoCoord[3 * i + 1] = V.y();
    geoCoord[3 * i + 2] = V.z();
  }
}

// Returns a vector containing the physical coordinates of the nodes on the
// element with LOCAL number "numElem" on the connectivity numbered "cncGeoTag" :
// coord = [x1 y1 z1 x2 y2 z2 ... xn yn zn]
void feMesh::getCoord(int cncGeoTag, int numElem, std::vector<double> &geoCoord)
{
  feCncGeo *cnc = getCncGeoByTag(cncGeoTag);
  int nNodePerElem = cnc->getNbNodePerElem();
#ifdef FENG_DEBUG
  if(geoCoord.size() != 3 * nNodePerElem) {
    printf(" In feMesh::getCoord : Erreur - Wrong Size for vector geoCoord\n");
  }
#endif
  for(int i = 0; i < nNodePerElem; ++i) {
    Vertex V = _vertices[cnc->getNodeConnectivity(numElem, i)];
    geoCoord[3 * i + 0] = V.x();
    geoCoord[3 * i + 1] = V.y();
    geoCoord[3 * i + 2] = V.z();
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

feMesh1DP1::feMesh1DP1(double xBegin, double xEnd, int numElements,
    std::string boundaryBeginName, std::string boundaryEndName, std::string domainName)
  : feMesh(numElements + 1, 1, 3, "1D")
  , _nElm(numElements)
  , _xA(xBegin)
  , _xB(xEnd)
  , _bndA_ID(boundaryBeginName)
  , _bndB_ID(boundaryEndName)
  , _domID(domainName)
  , _nElmPerBoundary(1)
  , _nNodDomain(2)
  , _nNodBoundary(1)
{
  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "MESH:");
  feInfoCond(FE_VERBOSE > 0, "\t\tCreating a 1D mesh for domain [%f, %f]", xBegin, xEnd);
  feInfoCond(FE_VERBOSE > 0, "\t\t\tNumber of 1D inner elements: %d", numElements);
  feInfoCond(FE_VERBOSE > 0, "\t\t\t1D element size (step): %f", (xEnd - xBegin) / (_nNod - 1));
  feInfoCond(FE_VERBOSE > 0, "\t\t\tNumber of 0D boundary elements: %d", 2);
  feInfoCond(FE_VERBOSE > 0, "\t\t\tNumber of vertices: %d", _nNod);
  feInfoCond(FE_VERBOSE > 0, "\t\t\tNumber of local connectivities: %d", 3);

  _nInteriorElm = numElements;
  _nBoundaryElm = 2;

  // Create the vertices
  _coord.resize(_nNod);
  for(int i = 0; i < _nNod; ++i) _coord[i] = xBegin + i * (xEnd - xBegin) / (_nNod - 1);
  _vertices.resize(_nNod);
  for(int i = 0; i < _nNod; ++i) _vertices[i] = Vertex(xBegin + i * (xEnd - xBegin) / (_nNod - 1), 0., 0., i);

  // Create the (trivial) elements connectivity
  int dimDomain = 1;
  std::vector<int> connecNodeDomain(_nInteriorElm * _nNodDomain, 0);
  for(int i = 0; i < _nInteriorElm; ++i) {
    connecNodeDomain[_nNodDomain * i + 0] = i;
    connecNodeDomain[_nNodDomain * i + 1] = i + 1;
  }

  // Create edges
  int nEdges = 1; // Numbering starts at 1 to match the 2D numbering
  std::vector<int> connecEdgeDomain(_nInteriorElm, 0);
  for(int i = 0; i < _nInteriorElm; ++i, ++nEdges) {
    Vertex *v0 = &_vertices[i];
    Vertex *v1 = &_vertices[i + 1];
    Edge e(v0, v1, nEdges, 0);
    _edges.insert(e);
    connecEdgeDomain[i] = nEdges; // Trivial edge connectivity but starting at 1
  }
  _nEdg = _nInteriorElm;

  // Create the local connectivity associated to the inner 1D domain
  int nCncGeo = 0;
  feCncGeo *geoDom = new feCncGeo(nCncGeo, dimDomain, _nNodDomain, _nInteriorElm, 1, _domID, "Lg",
    new feSpace1DP1("xyz"), connecNodeDomain, std::vector<int>(), connecEdgeDomain);
  _cncGeo.push_back(geoDom);
  _cncGeoMap[_domID] = nCncGeo;
  geoDom->getFeSpace()->setCncGeoTag(nCncGeo++);

  // Create the 0D boundary elements
  int dimBoundary = 0;
  std::vector<int> connecBoundaryA(_nElmPerBoundary * _nNodBoundary, 0);
  std::vector<int> connecBoundaryB(_nElmPerBoundary * _nNodBoundary, 0);
  connecBoundaryA[0] = 0;
  connecBoundaryB[0] = _nNod - 1;

  // Create the local connectivity for the left boundary
  feCncGeo *geoBndA = new feCncGeo(nCncGeo, dimBoundary, _nNodBoundary, _nElmPerBoundary, 0, _bndA_ID,
    "Point0D", new feSpace1DP0("xyz"), connecBoundaryA);
  _cncGeo.push_back(geoBndA);
  _cncGeoMap[_bndA_ID] = nCncGeo;
  geoBndA->getFeSpace()->setCncGeoTag(nCncGeo++);

  // Create the local connectivity for the right boundary
  feCncGeo *geoBndB = new feCncGeo(nCncGeo, dimBoundary, _nNodBoundary, _nElmPerBoundary, 0, _bndB_ID,
    "Point0D", new feSpace1DP0("xyz"), connecBoundaryB);
  _cncGeo.push_back(geoBndB);
  _cncGeoMap[_bndB_ID] = nCncGeo;
  geoBndB->getFeSpace()->setCncGeoTag(nCncGeo++);

  // Set the global element connectivity (a single continuous connectivity for the whole mesh)
  int numElmGlo = 0;
  for(feCncGeo *cnc : _cncGeo) {
    for(int i = 0; i < cnc->getNbElm(); ++i) cnc->setElementConnectivity(i, numElmGlo++);
  }
  _nTotalElm = numElmGlo;

  // Assign pointer to the mesh to cnc and their FE space:
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

// Locate physical point x in 1D mesh aligned with the x-axis.
// Assigns iElm the index of the element in which x lies, and
// u the vector (r, 0, 0) with r the coordinate of x in the reference
// [-1,1] element.
bool feMesh1DP1::locateVertex(const double *x, int &iElm, double *u, double tol)
{
  for(int i = 0; i < _nInteriorElm; ++i)
  {
    double x0 = _vertices[i].x();
    double x1 = _vertices[i + 1].x();
    if(x0 <= x[0] && x[0] <= x1){
      iElm = i;
      // Trivial inversion of the reference to physical transformation
      u[0] = (x[0] - (x0+x1)/2.) * 2. / (x1-x0);
      u[1] = 0.;
      u[2] = 0.;
      return true;
    }
  }
  feErrorMsg(FE_STATUS_ERROR, "Point %f was not found in 1D mesh :/", x[0]);
  return false;
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
  feCncGeo *geoDom = new feCncGeo(nCncGeo, dimDomain, _nNodDomain, _nElmDomain, 0, _domID, "Point0D",
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

// static bool rtreeCallback(int id, void *ctx)
// {
//   std::vector<int> *vec = reinterpret_cast<std::vector<int> *>(ctx);
//   vec->push_back(id);
//   return true;
// }

static bool rtreeCallback(int id, void *ctx)
{
  rtreeSearchCtx *searchCtx = reinterpret_cast<rtreeSearchCtx*>(ctx);
  Triangle *t = (*searchCtx->elements)[id];
  t->xyz2uvw(searchCtx->x, searchCtx->r);
  if(t->isInside(searchCtx->r[0], searchCtx->r[1], searchCtx->r[2])) {
    searchCtx->iElm = id;
    searchCtx->uvw[0] = searchCtx->r[0];
    searchCtx->uvw[1] = searchCtx->r[1];
    searchCtx->uvw[2] = searchCtx->r[2];
    searchCtx->isFound = true;
    return false;
  }
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
bool feMesh2DP1::locateVertex(const double *x, int &iElm, double *u, double tol)
{
  // double min[3] = {x[0] - tol, x[1] - tol, 0. - tol};
  // double max[3] = {x[0] + tol, x[1] + tol, 0. + tol};

  _searchCtx.min[0] = x[0] - tol;
  _searchCtx.min[1] = x[1] - tol;
  _searchCtx.max[0] = x[0] + tol;
  _searchCtx.max[1] = x[1] + tol;

  _searchCtx.x[0] = x[0];
  _searchCtx.x[1] = x[1];
  _searchCtx.x[2] = 0.0;

  _searchCtx.isFound = false;

  _rtree2d.Search(_searchCtx.min, _searchCtx.max, rtreeCallback, &_searchCtx);

  if(_searchCtx.isFound){
    iElm = _searchCtx.iElm;
    u[0] = _searchCtx.uvw[0];
    u[1] = _searchCtx.uvw[1];
    u[2] = _searchCtx.uvw[2];
  }

  return _searchCtx.isFound;
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
  for(size_t j = 0; j < _elements.size(); ++j) {
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
          feInfo("taille %d", 3 * fS2->getCncGeo()->getNbNodePerElem());
          std::vector<double> geoCoord(3 * fS2->getCncGeo()->getNbNodePerElem(),0.);
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
                      feInfo("iSol = %d - Accessing %d in sol1 of size %d from %d in solVec1 of size %d",
                        iSol, i, sol1.size(), adr1[i], solVec1.size());
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
