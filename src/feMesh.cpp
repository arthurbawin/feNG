#include "feNG.h"
#include "feMesh.h"
#include "feNumber.h"
#include "feTriangle.h"
#include "feSolution.h"

#include "SBoundingBox3d.h"

#include <iostream>
#include <fstream>

#if defined(HAVE_MPI)
  #include "mpi.h"
#endif

extern int FE_VERBOSE;

int feMesh::getCncGeoTag(const std::string &cncGeoID) const
{
  auto it = _cncGeoMap.find(cncGeoID);
  return (it == _cncGeoMap.end()) ? -1 : it->second;
}

feCncGeo *feMesh::getCncGeoByName(const std::string &cncGeoID) const
{
  auto it = _cncGeoMap.find(cncGeoID);
  return (it == _cncGeoMap.end()) ? nullptr : _cncGeo[it->second];
}

feCncGeo *feMesh::getCncGeoByTag(int cncGeoTag) const
{
  if(cncGeoTag < 0)
    if((unsigned)cncGeoTag > _cncGeo.size()) return nullptr;
  return _cncGeo[cncGeoTag];
}

void feMesh::recomputeDisplacementVectors()
{
  for(auto p : _edge2midnode) {
    double x0 = p.first->getVertex(0)->x();
    double y0 = p.first->getVertex(0)->y();
    double x1 = p.first->getVertex(1)->x();
    double y1 = p.first->getVertex(1)->y();
    double xMid = (x0+x1)/2.;
    double yMid = (y0+y1)/2.;
    _edge2alpha[p.first][0] = p.second->x() - xMid;
    _edge2alpha[p.first][1] = p.second->y() - yMid;
  }
}

double feMesh::getMaxNormEdgeDisplacement()
{
  if(_edge2alpha.empty()) return 0.;
  double alphaNormMax = 0.;
  for(auto p : _edge2alpha) {
    alphaNormMax = fmax(alphaNormMax, fabs(p.second[0]));
    alphaNormMax = fmax(alphaNormMax, fabs(p.second[1]));
  }
  return alphaNormMax;
}

double feMesh::getLpNormEdgeDisplacement(int p)
{
  if(_edge2alpha.empty()) return 0.;
  double norm = 0.;
  for(auto pair : _edge2alpha) {
    double normAlpha = sqrt(pair.second[0]*pair.second[0] + pair.second[1]*pair.second[1]);
    norm += pow(normAlpha, p);
  }
  return pow(norm, 1./ (double) p);
}

void feMesh::getCoord(const feCncGeo *cnc, const int numElem, std::vector<double> &geoCoord)
{
  int nVerticesPerElem = cnc->getNumVerticesPerElem();
#ifdef FENG_DEBUG
  if(geoCoord.size() != 3 * (size_t) nVerticesPerElem) {
    feErrorMsg(FE_STATUS_ERROR, " In feMesh::getCoord : Wrong size for vector geoCoord\n");
  }
#endif
  for(int i = 0; i < nVerticesPerElem; ++i) {
    Vertex V = _vertices[cnc->getVertexConnectivity(numElem, i)];
    geoCoord[3 * i + 0] = V.x();
    geoCoord[3 * i + 1] = V.y();
    geoCoord[3 * i + 2] = V.z();
  }
}

void feMesh::getCoord(std::string const &cncGeoID, const int numElem, std::vector<double> &geoCoord)
{
  feCncGeo *cnc = getCncGeoByName(cncGeoID);
  int nVerticesPerElem = cnc->getNumVerticesPerElem();
#ifdef FENG_DEBUG
  if(geoCoord.size() != 3 * (size_t) nVerticesPerElem) {
    printf(" In feMesh::getCoord : Wrong size for vector geoCoord\n");
  }
#endif
  for(int i = 0; i < nVerticesPerElem; ++i) {
    Vertex V = _vertices[cnc->getVertexConnectivity(numElem, i)];
    geoCoord[3 * i + 0] = V.x();
    geoCoord[3 * i + 1] = V.y();
    geoCoord[3 * i + 2] = V.z();
  }
}

void feMesh::getCoord(const int cncGeoTag, const int numElem, std::vector<double> &geoCoord)
{
  feCncGeo *cnc = getCncGeoByTag(cncGeoTag);
  int nVerticesPerElem = cnc->getNumVerticesPerElem();
#ifdef FENG_DEBUG
  if(geoCoord.size() != 3 * (size_t) nVerticesPerElem) {
    printf(" In feMesh::getCoord : Wrong Size for vector geoCoord\n");
  }
#endif
  for(int i = 0; i < nVerticesPerElem; ++i) {
    Vertex V = _vertices[cnc->getVertexConnectivity(numElem, i)];
    geoCoord[3 * i + 0] = V.x();
    geoCoord[3 * i + 1] = V.y();
    geoCoord[3 * i + 2] = V.z();
  }
}

void feMesh::getVertexCoord(const int iVertex, double coord[3])
{
#ifdef FENG_DEBUG
  if((size_t) iVertex > _vertices.size()) {
    feErrorMsg(FE_STATUS_ERROR,
               "Out of bounds: accessing entry %d in mesh->_vertices"
               " of size %u",
               iVertex, _vertices.size());
  }
#endif
  Vertex &v = _vertices[iVertex];
  coord[0] = v(0);
  coord[1] = v(1);
  coord[2] = v(2);
}

feStatus feMesh::setVertexCoord(const int iVertex, const double coord[3])
{
#ifdef FENG_DEBUG
  if((size_t) iVertex > _vertices.size()) {
    feErrorMsg(FE_STATUS_ERROR,
               "Out of bounds: accessing entry %d in mesh->_vertices"
               " of size %u",
               iVertex, _vertices.size());
  }
#endif
  Vertex &v = _vertices[iVertex];
  v(0) = coord[0];
  v(1) = coord[1];
  v(2) = coord[2];
  return FE_STATUS_OK;
}

int feMesh::getNumVerticesPerElem(std::string const &cncGeoID)
{
  return getCncGeoByName(cncGeoID)->getNumVerticesPerElem();
}

int feMesh::getNumVerticesPerElem(const int cncGeoTag)
{
  return getCncGeoByTag(cncGeoTag)->getNumVerticesPerElem();
}

int feMesh::getNumVerticesPerElemForDimension(const int dim) const
{
  for(const feCncGeo *cnc : _cncGeo) {
    if(cnc->getDim() == dim) {
      return cnc->getNumVerticesPerElem();
    }
  }
  return -1;
}

int feMesh::getVertex(std::string const &cncGeoID, const int numElem, const int numVertex)
{
  return getCncGeoByName(cncGeoID)->getVertexConnectivity(numElem, numVertex);
}

int feMesh::getVertex(const int cncGeoTag, const int numElem, const int numVertex)
{
  return getCncGeoByTag(cncGeoTag)->getVertexConnectivity(numElem, numVertex);
}

int feMesh::getElement(std::string const &cncGeoID, const int numElem)
{
  return getCncGeoByName(cncGeoID)->getElementConnectivity(numElem);
}

int feMesh::getElement(const int cncGeoTag, const int numElem)
{
  return getCncGeoByTag(cncGeoTag)->getElementConnectivity(numElem);
}

int feMesh::getEdge(std::string const &cncGeoID, const int numElem, const int numEdge)
{
  return getCncGeoByName(cncGeoID)->getEdgeConnectivity(numElem, numEdge);
}

int feMesh::getEdge(const int cncGeoTag, const int numElem, const int numEdge)
{
  return getCncGeoByTag(cncGeoTag)->getEdgeConnectivity(numElem, numEdge);
}

int feMesh::getFace(std::string const &cncGeoID, const int numElem, const int numFace)
{
  return getCncGeoByName(cncGeoID)->getFaceConnectivity(numElem, numFace);
}

int feMesh::getFace(const int cncGeoTag, const int numElem, const int numFace)
{
  return getCncGeoByTag(cncGeoTag)->getFaceConnectivity(numElem, numFace);
}

int feMesh::getNumElements(std::string const &cncGeoID)
{
  return getCncGeoByName(cncGeoID)->getNumElements();
}

int feMesh::getNumElements(int cncGeoTag) { return getCncGeoByTag(cncGeoTag)->getNumElements(); }

feSpace *feMesh::getGeometricSpace(std::string const &cncGeoID)
{
  return getCncGeoByName(cncGeoID)->getFeSpace();
}

feSpace *feMesh::getGeometricSpace(int cncGeoTag)
{
  return getCncGeoByTag(cncGeoTag)->getFeSpace();
}

feMesh1DP1::feMesh1DP1(const double xBegin, const double xEnd, const int numElements,
                       const std::string &boundaryBeginName, const std::string &boundaryEndName,
                       const std::string &domainName)
  : feMesh(numElements + 1, 1, 3, "1D"), _xA(xBegin), _xB(xEnd), _bndA_ID(boundaryBeginName),
    _bndB_ID(boundaryEndName), _domID(domainName), _nElmPerBoundary(1), _nNodDomain(2),
    _nNodBoundary(1)
{
  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "MESH:");
  feInfoCond(FE_VERBOSE > 0, "\t\tCreating a 1D mesh for domain [%f, %f]", xBegin, xEnd);
  feInfoCond(FE_VERBOSE > 0, "\t\t\tNumber of 1D inner elements: %d", numElements);
  feInfoCond(FE_VERBOSE > 0, "\t\t\t1D element size (step): %f",
             (xEnd - xBegin) / (_nVertices - 1));
  feInfoCond(FE_VERBOSE > 0, "\t\t\tNumber of 0D boundary elements: %d", 2);
  feInfoCond(FE_VERBOSE > 0, "\t\t\tNumber of vertices: %d", _nVertices);
  feInfoCond(FE_VERBOSE > 0, "\t\t\tNumber of local connectivities: %d", 3);

  _nInteriorElm = numElements;
  _nBoundaryElm = 2;

  // Create the vertices
  _vertices.resize(_nVertices);
  for(int i = 0; i < _nVertices; ++i)
    _vertices[i] = Vertex(xBegin + i * (xEnd - xBegin) / (_nVertices - 1), 0., 0., i);

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
  _nEdges = _nInteriorElm;

  // Create the local connectivity associated to the inner 1D domain
  int nCncGeo = 0;

  feCncGeo *geoDom =
    new feCncGeo(nCncGeo, dimDomain, 1, _nNodDomain, _nInteriorElm, 1, 0, _domID, geometryType::LINE,
                 geometricInterpolant::LINEP1, new feSpace1DP1("xyz"), connecNodeDomain,
                 std::vector<int>(), connecEdgeDomain);

  _cncGeo.push_back(geoDom);
  _cncGeoMap[_domID] = nCncGeo;
  geoDom->getFeSpace()->setCncGeoTag(nCncGeo++);

  // Create the 0D boundary elements
  int dimBoundary = 0;
  std::vector<int> connecBoundaryA(_nElmPerBoundary * _nNodBoundary, 0);
  std::vector<int> connecBoundaryB(_nElmPerBoundary * _nNodBoundary, 0);
  connecBoundaryA[0] = 0;
  connecBoundaryB[0] = _nVertices - 1;

  // Create the local connectivity for the left boundary
  feCncGeo *geoBndA = new feCncGeo(nCncGeo, dimBoundary, 1, _nNodBoundary, _nElmPerBoundary, 0, 0,
                                   _bndA_ID, geometryType::POINT, geometricInterpolant::POINTP0,
                                   new feSpace0DP0("xyz"), connecBoundaryA);

  _cncGeo.push_back(geoBndA);
  _cncGeoMap[_bndA_ID] = nCncGeo;
  geoBndA->getFeSpace()->setCncGeoTag(nCncGeo++);

  // Create the local connectivity for the right boundary
  feCncGeo *geoBndB = new feCncGeo(nCncGeo, dimBoundary, 1, _nNodBoundary, _nElmPerBoundary, 0, 0,
                                   _bndB_ID, geometryType::POINT, geometricInterpolant::POINTP0,
                                   new feSpace0DP0("xyz"), connecBoundaryB);

  _cncGeo.push_back(geoBndB);
  _cncGeoMap[_bndB_ID] = nCncGeo;
  geoBndB->getFeSpace()->setCncGeoTag(nCncGeo++);

  // Set the global element connectivity (a single continuous connectivity for the whole mesh)
  int numElmGlo = 0;
  for(feCncGeo *cnc : _cncGeo) {
    for(int i = 0; i < cnc->getNumElements(); ++i) cnc->setElementConnectivity(i, numElmGlo++);
  }
  _nTotalElm = numElmGlo;

  // Assign pointer to the mesh to cnc and their FE space:
  for(feCncGeo *cnc : _cncGeo) {
    cnc->setMeshPtr(this);
    cnc->getFeSpace()->setMeshPtr(this);
  }

  // Assign to each geometric FE space a pointer to its connectivity
  for(feCncGeo *cnc : _cncGeo) {
    cnc->getFeSpace()->setCncPtr(cnc);
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
bool feMesh1DP1::locateVertex(const double *x, int &iElm, double *u, double tol,
                              bool returnLocalElmTag, std::string targetConnectivity)
{
  UNUSED(tol, returnLocalElmTag, targetConnectivity);
  for(int i = 0; i < _nInteriorElm; ++i) {
    double x0 = _vertices[i].x();
    double x1 = _vertices[i + 1].x();
    if(x0 <= x[0] && x[0] <= x1) {
      iElm = i;
      // Trivial inversion of the reference to physical transformation
      u[0] = (x[0] - (x0 + x1) / 2.) * 2. / (x1 - x0);
      u[1] = 0.;
      u[2] = 0.;
      return true;
    }
  }
  feErrorMsg(FE_STATUS_ERROR, "Point %f was not found in 1D mesh :/", x[0]);
  return false;
}

bool feMesh1DP1::locateVertexInElements(feCncGeo *cnc, const double *x, const std::vector<int> &elementsToSearch, int &iElm,
                                      double *u, double tol)
{
  UNUSED(cnc, x, elementsToSearch, iElm, u, tol);
  feErrorMsg(FE_STATUS_ERROR, "locateVertexInElements not implemented for 1D meshes :/");
  return false;
}

feMesh0DP0::feMesh0DP0(double xA, int nElm, std::string domID)
  : feMesh(1, 0, 1, "0D"), _nElm(nElm), _xA(xA), _domID(domID), _nElmDomain(nElm), _nElmBoundary(1),
    _nNodDomain(1), _nNodBoundary(1)
{
  // Create the vertex
  _vertices.resize(_nVertices);
  for(int i = 0; i < _nVertices; ++i) {
    _vertices[i] = Vertex(0, 0., 0., 0);
  }

  // Element 0D
  int dimDomain = 0;
  std::vector<int> connecDomain(1, 0);

  int nCncGeo = 0;
  feCncGeo *geoDom =
    new feCncGeo(nCncGeo, dimDomain, 0, _nNodDomain, _nElmDomain, 0, 0, _domID, geometryType::POINT,
                 geometricInterpolant::POINTP0, new feSpace1DP0("xyz"), connecDomain);
  _cncGeo.push_back(geoDom);
  _cncGeoMap[_domID] = nCncGeo;
  geoDom->getFeSpace()->setCncGeoTag(nCncGeo++);

  // Global connectivity
  int numElmGlo = 0;
  for(feCncGeo *cnc : _cncGeo) {
    for(int i = 0; i < cnc->getNumElements(); ++i) cnc->setElementConnectivity(i, numElmGlo++);
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

feMesh2DP1::feMesh2DP1(const std::string &meshName, const bool curved, const bool reversed,
                       const mapType &physicalEntitiesDescription)
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

double feMesh2DP1::computeAverageSliverness(const double LpNorm)
{
  FILE *myfile = fopen("meshSliverness.txt", "w");
  fprintf(myfile, "View \"sliverness\"{\n");

  double res = 0.;
  for(Triangle *t : _elements) {
    double s = t->sliverness();
    res += pow(s, LpNorm);

    fprintf(myfile, "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g){%g,%g,%g};\n",
      t->getVertex(0)->x(), t->getVertex(0)->y(), 0.,
      t->getVertex(1)->x(), t->getVertex(1)->y(), 0.,
      t->getVertex(2)->x(), t->getVertex(2)->y(), 0., s, s, s);
  }

  fprintf(myfile, "};\n"); fclose(myfile);

  return pow(res/(double)_elements.size(), 1./LpNorm);
}

// Callback to give the RTree to locate a physical point in the mesh.
bool rtreeCallback(int id, void *ctx)
{
  rtreeSearchCtx *searchCtx = reinterpret_cast<rtreeSearchCtx *>(ctx);
  Triangle *t = (*(searchCtx->elements))[id];

  // Get the reference coordinates (u,v) such that F_K(u,v) = (x,y),
  // where F_K is the reference-to-physical transformation of element K.
  // For P1 triangle, success is always true.
  // For P2 triangle, success is true if the Newton-Raphson method converged
  // to a pair (u,v).
  // If success is false, the point should be outside of the P2 triangle.
  //
  // When computing (u,v) in P2 triangle, use smaller tolerance for the Newton
  bool success = t->xyz2uvw(searchCtx->x, searchCtx->r, searchCtx->tolerance / 10.);
  if(!success) {
    // Keep looking
    return true;
  }

  if(t->isInsideReference(searchCtx->r[0], searchCtx->r[1], searchCtx->r[2], searchCtx->tolerance)) {
    // Check if vertex was found on prescribed connectivity, if applicable
    if(searchCtx->enforceConnectivity && (t->getPhysicalTag() != searchCtx->targetPhysicalTag)) {
      // Keep looking
      return true;
    }
    searchCtx->iElm = id;
    searchCtx->iElmLocal = t->getLocalTag();
    searchCtx->physicalTag = t->getPhysicalTag();
    searchCtx->uvw[0] = searchCtx->r[0];
    searchCtx->uvw[1] = searchCtx->r[1];
    searchCtx->uvw[2] = searchCtx->r[2];
    searchCtx->isFound = true;
    return false;
  }
  return true;
}

thread_local rtreeSearchCtx SEARCH_CONTEXT;

/* Locates the vertex with coordinates x in the mesh using an RTree.
   The search is performed in elements of the highest dimension only.
   The element number is assigned to iElm and the reference coordinates
   are assigned in u. */
bool feMesh2DP1::locateVertex(const double *x, int &iElm, double *u, double tol,
                              bool returnLocalElmTag, std::string targetConnectivity)
{
  SEARCH_CONTEXT.elements = &(this->_elements);
  SEARCH_CONTEXT.tolerance = tol;
  SEARCH_CONTEXT.min[0] = x[0] - tol;
  SEARCH_CONTEXT.min[1] = x[1] - tol;
  SEARCH_CONTEXT.max[0] = x[0] + tol;
  SEARCH_CONTEXT.max[1] = x[1] + tol;

  SEARCH_CONTEXT.uvw[0] = -1.;
  SEARCH_CONTEXT.uvw[1] = -1.;
  SEARCH_CONTEXT.uvw[2] = -1.;

  SEARCH_CONTEXT.x[0] = x[0];
  SEARCH_CONTEXT.x[1] = x[1];
  SEARCH_CONTEXT.x[2] = 0.0;

  SEARCH_CONTEXT.isFound = false;
  SEARCH_CONTEXT.numFound = 0;

  if(targetConnectivity != "") {
    SEARCH_CONTEXT.enforceConnectivity = true;

    bool OK = false;
    for(auto &p : _physicalEntities) {
      if(p.first.first == 2) {
        if(p.second.name == targetConnectivity) {
          SEARCH_CONTEXT.targetPhysicalTag = p.first.second;
          OK = true;
        }
      }
    }
    if(!OK) {
      feErrorMsg(FE_STATUS_ERROR,
                 "Cannot look for vertex on connectivity %s because it has no matching physical "
                 "tag in the mesh.",
                 targetConnectivity.data());
    }
  }

  _rtree2d.Search(SEARCH_CONTEXT.min, SEARCH_CONTEXT.max, rtreeCallback, &SEARCH_CONTEXT);

  if(SEARCH_CONTEXT.isFound) {
    if(returnLocalElmTag)
      iElm = SEARCH_CONTEXT.iElmLocal;
    else
      iElm = SEARCH_CONTEXT.iElm;

    u[0] = SEARCH_CONTEXT.uvw[0];
    u[1] = SEARCH_CONTEXT.uvw[1];
    u[2] = SEARCH_CONTEXT.uvw[2];

    // Check if vertex was found on the prescribed connectivity
    // if(targetConnectivity != "") {
    //   for(auto p : _physicalEntities) {
    //     if(p.first.first == 2 && p.first.second == SEARCH_CONTEXT.physicalTag) {
    //       feInfo("Found vertex (%f - %f) on physical entity %d named %s", x[0], x[1],
    //              SEARCH_CONTEXT.physicalTag, p.second.name.data());
    //     }
    //   }
    // }
  }

  return SEARCH_CONTEXT.isFound;
}

// From gmsh and feTriangle.cpp
void xyz2uvw(const double v0[2], const double v1[2], const double v2[2], const double xyz[2], double uvw[2])
{
  const double O[3] = {v0[0], v0[1], 0.};

  const double d[3] = {xyz[0] - O[0], xyz[1] - O[1], 0.};
  const double d1[3] = {v1[0] - O[0], v1[1] - O[1], 0.};
  const double d2[3] = {v2[0] - O[0], v2[1] - O[1], 0.};

  const double Jxy = d1[0] * d2[1] - d1[1] * d2[0];
  const double Jxz = d1[0] * d2[2] - d1[2] * d2[0];
  const double Jyz = d1[1] * d2[2] - d1[2] * d2[1];

  if((fabs(Jxy) > fabs(Jxz)) && (fabs(Jxy) > fabs(Jyz))) {
    uvw[0] = (d[0] * d2[1] - d[1] * d2[0]) / Jxy;
    uvw[1] = (d[1] * d1[0] - d[0] * d1[1]) / Jxy;
  } else if(fabs(Jxz) > fabs(Jyz)) {
    uvw[0] = (d[0] * d2[2] - d[2] * d2[0]) / Jxz;
    uvw[1] = (d[2] * d1[0] - d[0] * d1[2]) / Jxz;
  } else {
    uvw[0] = (d[1] * d2[2] - d[2] * d2[1]) / Jyz;
    uvw[1] = (d[2] * d1[1] - d[1] * d1[2]) / Jyz;
  }
}

bool isInsideTriangle(double u, double v, double w)
{
  double tol = 1e-12;
  if(u < (-tol) || v < (-tol) || u > ((1. + tol) - v) || fabs(w) > tol) return false;
  return true;
}

thread_local std::vector<double> ELEM_COORD(18,0.);

bool feMesh2DP1::locateVertexInElements(feCncGeo *cnc, const double *x, const std::vector<int> &elementsToSearch, int &iElm,
                                      double *u, double tol)
{
  UNUSED(tol);
  size_t numElem = elementsToSearch.size();
  for(size_t i = 0; i < numElem; ++i) {
    feInfo("Searching %d/%d elements: %d", i+1, numElem, elementsToSearch[i]);
    this->getCoord(cnc, elementsToSearch[i], ELEM_COORD);
    double v0[2] = {ELEM_COORD[0], ELEM_COORD[1]};
    double v1[2] = {ELEM_COORD[3], ELEM_COORD[4]};
    double v2[2] = {ELEM_COORD[6], ELEM_COORD[7]};
    xyz2uvw(v0, v1, v2, x, u);
    if(isInsideTriangle(u[0], u[1], 0.)) {
      iElm = elementsToSearch[i];
      return true;
    }
  }
  return false;
}

// Transfer from the active feMesh (source) to the targetMesh (destination)
feStatus feMesh2DP1::transfer(feMesh2DP1 *targetMesh, feMetaNumber *prevNumbering, feMetaNumber *targetNumbering,
                          feSolutionContainer *prevSolutionContainer,
                          const std::vector<feSpace *> &mySpaces,
                          const std::vector<feSpace *> &mySpacesEssBC,
                          const std::vector<feSpace *> &targetSpaces)
{
  UNUSED(prevNumbering);
  
  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "SOLUTION PROJECTION:");

  if(targetMesh != targetSpaces[0]->getMeshPtr()) {
    return feErrorMsg(FE_STATUS_ERROR, "The target mesh for solution interpolation does not match the mesh pointer of the first FE space."
                                       "The new FE spaces should be initialized with the targetMesh.");
  }

  // A temporary solutionContainer in which the interpolated values are stored
  feSolutionContainer scTmp(*prevSolutionContainer);

  int nSol = scTmp.getNbSol();
  int nDOF = targetNumbering->getNbDOFs();

  // Resize to the new number of DOFs
  for(int iSol = 0; iSol < nSol; ++iSol) {
    scTmp._sol[iSol].resize(nDOF);
  }

  double tol = 1e-12;

  // Check which spaces need to be interpolated
  std::vector<std::pair<std::string, std::string> > fieldsToInterpolate;
  for(feSpace *fS1 : mySpaces)
  {
    bool isBC = false;

    // Don't transfer the 0D-1D connectivities for now
    if(fS1->getDim() < this->_dim) {
      feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE, "\t\tNot interpolating boundary field of dimension %d", fS1->getDim());
      continue;
    }

    // Don't transfer essential spaces which are re-initialized in the adaptation loop
    for(feSpace *feEss : mySpacesEssBC) {
      if(fS1->getFieldID() == feEss->getFieldID() && fS1->getCncGeoID() == feEss->getCncGeoID()) {
        isBC = true;
        feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE,
          "\t\tNot interpolating essential BC field \"%s\" on connectivity \"%s\"", 
          fS1->getFieldID().data(),
          fS1->getCncGeoID().data());
      }
    }

    std::vector<double> x(3, 0.);
    std::vector<int> candidates;

    if(!isBC)
    {
      for(feSpace *fS2 : targetSpaces)
      {
        if(fS1->getFieldID() == fS2->getFieldID() && fS1->getCncGeoID() == fS2->getCncGeoID())
        {
          // Check that both fields are scalar-valued or vector-valued
          if(fS1->getNumComponents() != fS2->getNumComponents()) {
            return feErrorMsg(FE_STATUS_ERROR, "Cannot interpolate field \"%s\" on connectivity \"%s\""
              " because source has %d vector component(s) and target has %d vector component(s)",
              fS1->getFieldID().data(), fS1->getCncGeoID().data(),
              fS1->getNumComponents(), fS2->getNumComponents());
          }

          feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE, "\t\tInterpolating field \"%s\" on connectivity \"%s\"",
            fS1->getFieldID().data(),
            fS1->getCncGeoID().data());

          fieldsToInterpolate.emplace_back(fS1->getFieldID(), fS1->getCncGeoID());

          // Get the target connectivity from the target mesh
          feCncGeo *targetCnc = targetMesh->getCncGeoByName(fS2->getCncGeoID());
          if(targetCnc == nullptr) {
            return feErrorMsg(FE_STATUS_ERROR, "Target connectivity \"%s\" on target mesh is null :/", fS2->getCncGeoID().data());
          }

          int nElm = targetCnc->getNumElements();
          int nDOFPerElem = fS2->getNumFunctions();
          std::string cncGeoID = fS2->getCncGeoID();
          const std::vector<double> &Lcoor = fS2->getLcoor();
          feSpace *geoSpace2 = fS2->getCncGeo()->getFeSpace();

          std::vector<feInt> adr1(fS1->getNumFunctions());
          std::vector<feInt> adr2(fS2->getNumFunctions());
          std::vector<double> sol1(adr1.size());

          std::vector<double> geoCoord(3 * fS2->getNumVerticesPerElem(), 0.);

          for(int iElm = 0; iElm < nElm; ++iElm)
          {
            fS2->initializeAddressingVector(iElm, adr2);
            targetMesh->getCoord(cncGeoID, iElm, geoCoord);

            // Loop over the DOFs of the element of the new mesh
            for(int j = 0; j < nDOFPerElem; ++j)
            {
              // Locate DOF in old mesh (*this):
              // Interpolate coordinates in new mesh using geoSpace2
              double xsi[3] = {Lcoor[3 * j + 0], Lcoor[3 * j + 1], Lcoor[3 * j + 2]};
              geoSpace2->interpolateVectorField(geoCoord, xsi, x);

              int elmInOldMesh;
              this->locateVertex(x.data(), elmInOldMesh, xsi, tol);

              Triangle *t = _elements[elmInOldMesh];
              fS1->initializeAddressingVector(t->getTag(), adr1);

              // Interpolate each solution in the container (one solution for each BDF level)
              for(int iSol = 0; iSol < nSol; ++iSol) {
                const std::vector<double> &solCopy = prevSolutionContainer->getSolution(iSol);
                for(size_t i = 0; i < adr1.size(); ++i) {
                  sol1[i] = solCopy[adr1[i]];
                }

                // Works for both scalar and vector fields
                int iComponent = j % fS1->getNumComponents();
                double solInt = fS1->interpolateVectorFieldComponent(sol1, iComponent, xsi);
                scTmp._sol[iSol][adr2[j]] = solInt;
                // TODO : Interpoler le _fresidual du solutionContainer
              }
            }
          } // for iElm
        } else { // if fields match
        }
      }
    }
  }

  // Resize the solutionContainer
  prevSolutionContainer->setNbDOFs(nDOF);
  prevSolutionContainer->_d.resize(nDOF);

  // Replace solutionContainer by scTmp
  for(int iSol = 0; iSol < nSol; ++iSol) {
    prevSolutionContainer->_sol[iSol].resize(nDOF);
    prevSolutionContainer->_solDot[iSol].resize(nDOF);
    prevSolutionContainer->_fResidual[iSol].resize(nDOF);
    for(int iDOF = 0; iDOF < nDOF; ++iDOF) {
      prevSolutionContainer->_sol[iSol][iDOF] = scTmp._sol[iSol][iDOF];
      // prevSolutionContainer->_solDot[iSol][iDOF] = scTmp._solDot[iSol][iDOF];
      prevSolutionContainer->_solDot[iSol][iDOF] = nan("");
    }
  }

  return FE_STATUS_OK;
}