#include "feSpace.h"
#include "feMesh.h"
#include "feNumber.h"
#include "fePolynomial.h"

EigenMat I6 = EigenMat::Identity(6, 6);

// Copy array "dim" times by chunks of size "chunksize"
static void duplicateScalarArray(const int scalarSize, const int chunkSize,
                                 const double *scalarArray, const int dim,
                                 std::vector<double> &vectorArray)
{
  int nChunks = scalarSize / chunkSize;
  int cnt = 0;
  for(int j = 0; j < nChunks; ++j) {
    for(int iDim = 0; iDim < dim; ++iDim) {
      for(int i = 0; i < chunkSize; ++i) {
        vectorArray[cnt++] = scalarArray[chunkSize * j + i];
      }
    }
  }
}

static void duplicateScalarArray(const int scalarSize, const int chunkSize,
                                 const double *scalarArray, const int dim, double *vectorArray)
{
  int nChunks = scalarSize / chunkSize;
  int cnt = 0;
  for(int j = 0; j < nChunks; ++j) {
    for(int iDim = 0; iDim < dim; ++iDim) {
      for(int i = 0; i < chunkSize; ++i) {
        feInfo("Writing entry %d of scalararray of size %d at place %d", chunkSize * j + i,
               scalarSize, cnt + 1);
        vectorArray[cnt++] = scalarArray[chunkSize * j + i];
      }
    }
  }
}

static void vectorLagrangeLayout(const double *phi, const int nScalarFunctions,
  const int nComponents, double *res)
{
  int cnt = 0;
  for(int i = 0; i < nScalarFunctions; ++i) {
    for(int j = 0; j < nComponents; ++j) {
      for(int k = 0; k < nComponents; ++k) {
        if(j == k)
          res[cnt++] = phi[i];
        else
          res[cnt++] = 0.;
      }
    }
  }
}

// -----------------------------------------------------------------------------
// Discontinuous Lagrange element of degree 0 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
feSpaceTriP0::feSpaceTriP0(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                           feFunction *fct, const bool useGlobalShapeFunctions)
  : feScalarSpace(2, mesh, fieldID, cncGeoID, fct, useGlobalShapeFunctions)
{
  _isDiscontinuous = true;
  _nFunctions = 1;
  _Lcoor = {1./3., 1./3., 0.};
  _dofLocations = {dofLocation::ELEMENT};
}

std::vector<double> feSpaceTriP0::L(double */*r*/) { return {1.0}; }
void feSpaceTriP0::L(double */*r*/, double *L) { L[0] = 1.0; }
std::vector<double> feSpaceTriP0::dLdr(double */*r*/) { return {0.}; }
std::vector<double> feSpaceTriP0::dLds(double */*r*/) { return {0.}; }
std::vector<double> feSpaceTriP0::d2Ldr2(double */*r*/) { return {0.}; }
std::vector<double> feSpaceTriP0::d2Ldrs(double */*r*/) { return {0.}; }
std::vector<double> feSpaceTriP0::d2Lds2(double */*r*/) { return {0.}; }

void feSpaceTriP0::initializeNumberingUnknowns()
{
  int nDOFPerElem = 1;
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setUnknownElementDOF(_mesh, _cncGeoID, i, nDOFPerElem);
  }
}

void feSpaceTriP0::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialElementDOF(_mesh, _cncGeoID, i);
  }
}

void feSpaceTriP0::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
}

// -----------------------------------------------------------------------------
// Lagrange element of degree 1 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
feSpaceTriP1::feSpaceTriP1(std::string cncGeoID)
  : feScalarSpace(2, nullptr, "GEO", cncGeoID, nullptr)
{
  _isGeometricInterpolant = true;
  _nFunctions = 3;
  _Lcoor = {0., 0., 0., 1., 0., 0., 0., 1., 0.};
  _dofLocations = {dofLocation::VERTEX, dofLocation::VERTEX, dofLocation::VERTEX};
}

feSpaceTriP1::feSpaceTriP1(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                           feFunction *fct, const bool useGlobalShapeFunctions)
  : feScalarSpace(2, mesh, fieldID, cncGeoID, fct, useGlobalShapeFunctions)
{
  _nFunctions = 3;
  _Lcoor = {0., 0., 0., 1., 0., 0., 0., 1., 0.};
  _dofLocations = {dofLocation::VERTEX, dofLocation::VERTEX, dofLocation::VERTEX};
}

std::vector<double> feSpaceTriP1::L(double *r) { return {1.0 - r[0] - r[1], r[0], r[1]}; }
void feSpaceTriP1::L(double *r, double *L)
{
  L[0] = 1.0 - r[0] - r[1];
  L[1] = r[0];
  L[2] = r[1];
}
std::vector<double> feSpaceTriP1::dLdr(double *r) { UNUSED(r); return {-1.0, 1.0, 0.0}; }
std::vector<double> feSpaceTriP1::dLds(double *r) { UNUSED(r); return {-1.0, 0.0, 1.0}; }
std::vector<double> feSpaceTriP1::d2Ldr2(double *r) { UNUSED(r); return {0., 0., 0.}; }
std::vector<double> feSpaceTriP1::d2Ldrs(double *r) { UNUSED(r); return {0., 0., 0.}; }
std::vector<double> feSpaceTriP1::d2Lds2(double *r) { UNUSED(r); return {0., 0., 0.}; }

void feSpaceTriP1::initializeNumberingUnknowns()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 2);
  }
}

void feSpaceTriP1::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 2);
  }
}

void feSpaceTriP1::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2);
}

// -----------------------------------------------------------------------------
// Discontinuous Lagrange element of degree 1 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
feSpaceTriP1_Discontinuous::feSpaceTriP1_Discontinuous(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                           feFunction *fct)
  : feScalarSpace(2, mesh, fieldID, cncGeoID, fct)
{
  _isDiscontinuous = true;
  _nFunctions = 3;
  _Lcoor = {0., 0., 0., 1., 0., 0., 0., 1., 0.};
  _dofLocations = {dofLocation::VERTEX, dofLocation::VERTEX, dofLocation::VERTEX};
}

std::vector<double> feSpaceTriP1_Discontinuous::L(double *r) { return {1.0 - r[0] - r[1], r[0], r[1]}; }
void feSpaceTriP1_Discontinuous::L(double *r, double *L)
{
  L[0] = 1.0 - r[0] - r[1];
  L[1] = r[0];
  L[2] = r[1];
}
std::vector<double> feSpaceTriP1_Discontinuous::dLdr(double *r) { UNUSED(r); return {-1.0, 1.0, 0.0}; }
std::vector<double> feSpaceTriP1_Discontinuous::dLds(double *r) { UNUSED(r); return {-1.0, 0.0, 1.0}; }
std::vector<double> feSpaceTriP1_Discontinuous::d2Ldr2(double *r) { UNUSED(r); return {0., 0., 0.}; }
std::vector<double> feSpaceTriP1_Discontinuous::d2Ldrs(double *r) { UNUSED(r); return {0., 0., 0.}; }
std::vector<double> feSpaceTriP1_Discontinuous::d2Lds2(double *r) { UNUSED(r); return {0., 0., 0.}; }

void feSpaceTriP1_Discontinuous::initializeNumberingUnknowns()
{
  int nDOFPerElem = 3;
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setUnknownElementDOF(_mesh, _cncGeoID, i, nDOFPerElem);
  }
}

void feSpaceTriP1_Discontinuous::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialElementDOF(_mesh, _cncGeoID, i);
  }
}

void feSpaceTriP1_Discontinuous::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  int nDOFPerElem = 3;
  for(int i = 0; i < nDOFPerElem; ++i) {
    adr[i] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, i);
  }
}

void feSpaceTriP1_Discontinuous::synchronizeCodeOfEssentialDOF()
{
  if(_numbering->hasVertexDOFs()) {
    int numDOFAssociatedToVertices = 3;
    // For each element, request the DOF status (unknown or essential)
    // at each vertex, then update the corresponding element DOF
    // if this DOF has been marked as essential by another feSpace.
    for(int iElm = 0; iElm < _mesh->getNumElements(_cncGeoID); ++iElm) {
      for(int j = 0; j < numDOFAssociatedToVertices; ++j) {
        int vert = _mesh->getVertex(_cncGeoID, iElm, j);
        int vertexDOFtype = _numbering->getDOFCodeAtVertex(vert);
        if(vertexDOFtype == dofType::DOF_ESSENTIAL) {
          _numbering->setEssentialElementDOF(_mesh, _cncGeoID, iElm, j);
        }
      }
    }
  }
};

void feSpaceTriP1_Discontinuous::synchronizeNumberingOfEssentialDOF(int &numModifiedDOF)
{
  if(_numbering->hasVertexDOFs()) {
    int numDOFAssociatedToVertices = 3;
    // If vertex DOF is essential, assign its DOF number from the 
    // list of vertex DOF numbers.
    for(int iElm = 0; iElm < _mesh->getNumElements(_cncGeoID); ++iElm) {
      for(int j = 0; j < numDOFAssociatedToVertices; ++j) {
        int vert = _mesh->getVertex(_cncGeoID, iElm, j);
        int vertexDOFtype = _numbering->getDOFCodeAtVertex(vert);
        if(vertexDOFtype == dofType::DOF_ESSENTIAL) {
          int dofNumber = _numbering->getDOFNumberAtVertex(vert, 0); // Is it always 0? Probably for scalar spaces
          _numbering->setElementDOF(_mesh, _cncGeoID, iElm, j, dofNumber);
          numModifiedDOF++;
        }
      }
    }
  }
};

// -----------------------------------------------------------------------------
// Vector Lagrange element of degree 1 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
template <int dim>
feSpaceVecTriP1<dim>::feSpaceVecTriP1(feMesh *mesh, const std::string fieldID,
                                      const std::string cncGeoID, feVectorFunction *fct,
                                      const bool useGlobalShapeFunctions)
  : feVectorSpace(2, mesh, fieldID, cncGeoID, fct, useGlobalShapeFunctions)
{
  _nComponents = dim;
  _nFunctions = 3 * dim;
  double coor[9] = {0., 0., 0., 1., 0., 0., 0., 1., 0.};
  _Lcoor.resize(9 * dim);
  duplicateScalarArray(9, 3, coor, dim, _Lcoor);
  _dofLocations.resize(_nFunctions);
  for(int i = 0; i < _nFunctions; ++i) _dofLocations[i] = dofLocation::VERTEX;
}

template <int dim> std::vector<double> feSpaceVecTriP1<dim>::L(double *r)
{
  double phi[3] = {1.0 - r[0] - r[1], r[0], r[1]};
  std::vector<double> res(_nFunctions*_nComponents);
  // duplicateScalarArray(3, 1, phi, dim, res);
  vectorLagrangeLayout(phi, 3, _nComponents, res.data());
  return res;
}

template <int dim> void feSpaceVecTriP1<dim>::L(double *r, double *res)
{
  double phi[3] = {1.0 - r[0] - r[1], r[0], r[1]};
  // duplicateScalarArray(3, 1, phi, dim, res);
  vectorLagrangeLayout(phi, 3, _nComponents, res);
}

template <int dim> std::vector<double> feSpaceVecTriP1<dim>::dLdr(double */*r*/)
{
  double dldr[3] = {-1.0, 1.0, 0.0};
  std::vector<double> res(_nFunctions*_nComponents);
  // duplicateScalarArray(3, 1, dldr, dim, res);
  vectorLagrangeLayout(dldr, 3, _nComponents, res.data());
  return res;
}

template <int dim> std::vector<double> feSpaceVecTriP1<dim>::dLds(double */*r*/)
{
  double dlds[3] = {-1.0, 0.0, 1.0};
  std::vector<double> res(_nFunctions*_nComponents);
  // duplicateScalarArray(3, 1, dlds, dim, res);
  vectorLagrangeLayout(dlds, 3, _nComponents, res.data());
  return res;
}

template <int dim> std::vector<double> feSpaceVecTriP1<dim>::d2Ldr2(double */*r*/)
{
  double d2ldr2[3] = {0., 0., 0.};
  std::vector<double> res(_nFunctions*_nComponents);
  // duplicateScalarArray(3, 1, d2ldr2, dim, res);
  vectorLagrangeLayout(d2ldr2, 3, _nComponents, res.data());
  return res;
}

template <int dim> std::vector<double> feSpaceVecTriP1<dim>::d2Ldrs(double */*r*/)
{
  double d2ldrs[3] = {0., 0., 0.};
  std::vector<double> res(_nFunctions*_nComponents);
  // duplicateScalarArray(3, 1, d2ldrs, dim, res);
  vectorLagrangeLayout(d2ldrs, 3, _nComponents, res.data());
  return res;
}

template <int dim> std::vector<double> feSpaceVecTriP1<dim>::d2Lds2(double */*r*/)
{
  double d2lds2[3] = {0., 0., 0.};
  std::vector<double> res(_nFunctions*_nComponents);
  // duplicateScalarArray(3, 1, d2lds2, dim, res);
  vectorLagrangeLayout(d2lds2, 3, _nComponents, res.data());
  return res;
}

template <int dim> void feSpaceVecTriP1<dim>::initializeNumberingUnknowns()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 0, dim);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 1, dim);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 2, dim);
  }
}

template <int dim> void feSpaceVecTriP1<dim>::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 2);
  }
}

template <int dim>
void feSpaceVecTriP1<dim>::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 1);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 0);
  adr[3] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 1);
  adr[4] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2, 0);
  adr[5] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2, 1);
}

template class feSpaceVecTriP1<1>;
template class feSpaceVecTriP1<2>;

// -----------------------------------------------------------------------------
// Non-conforming Crouzeix-Raviart element of degree 1 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
feSpaceTri_CR1::feSpaceTri_CR1(feMesh *mesh, std::string fieldID, std::string cncGeoID,
                               feFunction *fct)
  : feScalarSpace(2, mesh, fieldID, cncGeoID, fct)
{
  _nFunctions = 3;
  _Lcoor = {0.5, 0., 0., 0.5, 0.5, 0., 0., 0.5, 0.};
  _dofLocations = {dofLocation::EDGE, dofLocation::EDGE, dofLocation::EDGE};
}

std::vector<double> feSpaceTri_CR1::L(double *r)
{
  return {1 - 2.0 * r[1], -1.0 + 2.0 * r[0] + 2.0 * r[1], 1 - 2.0 * r[0]};
}
void feSpaceTri_CR1::L(double *r, double *L)
{
  L[0] = 1 - 2.0 * r[1];
  L[1] = -1.0 + 2.0 * r[0] + 2.0 * r[1];
  L[2] = 1 - 2.0 * r[0];
}
std::vector<double> feSpaceTri_CR1::dLdr(double *r) { UNUSED(r); return {0.0, 2.0, -2.0}; }
std::vector<double> feSpaceTri_CR1::dLds(double *r) { UNUSED(r); return {-2.0, 2.0, 0.0}; }
std::vector<double> feSpaceTri_CR1::d2Ldr2(double *r) { UNUSED(r); return {0., 0., 0.}; }
std::vector<double> feSpaceTri_CR1::d2Ldrs(double *r) { UNUSED(r); return {0., 0., 0.}; }
std::vector<double> feSpaceTri_CR1::d2Lds2(double *r) { UNUSED(r); return {0., 0., 0.}; }

void feSpaceTri_CR1::initializeNumberingUnknowns()
{
  int nDOFPerEdge = 1;
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, nDOFPerEdge);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 1, nDOFPerEdge);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 2, nDOFPerEdge);
  }
}

void feSpaceTri_CR1::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 2);
  }
}

void feSpaceTri_CR1::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 0);
  adr[1] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 0);
  adr[2] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 0);
}

// -----------------------------------------------------------------------------
// Raviart-Thomas (vector-valued) element of degree 1 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
feSpaceTriRT1::feSpaceTriRT1(feMesh *mesh, const std::string fieldID,
                             const std::string cncGeoID, feVectorFunction *fct,
                             const bool useGlobalShapeFunctions)
  : feVectorSpace(2, mesh, fieldID, cncGeoID, fct, useGlobalShapeFunctions)
{
  _elementType = elementType::RAVIART_THOMAS;
  _isDiscontinuous = true;
  _nComponents = 2;
  _nFunctions = 3;
  // const double coord[9] = {0.5, 0.5, 0.0,
  //                          0.0, 0.5, 0.0,
  //                          0.5, 0.0, 0.0};
  // _Lcoor.resize(9 * _dim);
  // duplicateScalarArray(9, 3, coord, _dim, _Lcoor);
  _Lcoor = {0.5, 0.5, 0.0,
            0.0, 0.5, 0.0,
            0.5, 0.0, 0.0};
  _dofLocations.resize(_nFunctions);
  for(int i = 0; i < _nFunctions; ++i) _dofLocations[i] = dofLocation::EDGE;
}

std::vector<double> feSpaceTriRT1::L(double *r)
{
  // These are the shape functions from Guermond
  // The shape functions from defelement may require
  // sign changes for the normals
  return {r[0], r[1], r[0]-1., r[1], r[0], -1.+r[1]};
}

void feSpaceTriRT1::L(double *r, double *res)
{
  res[0] =   -r[0]; res[1] =   -r[1];
  res[2] = r[0]-1.; res[3] =    r[1];
  res[4] =   -r[0]; res[5] = 1.-r[1];
}

std::vector<double> feSpaceTriRT1::dLdr(double */*r*/)
{
  return {-1., 0., 1., 0., -1., 0.};
}

std::vector<double> feSpaceTriRT1::dLds(double */*r*/)
{
  return {0., -1., 0., 1., 0., -1.};
}

std::vector<double> feSpaceTriRT1::d2Ldr2(double */*r*/)
{
  return {0., 0., 0., 0., 0., 0.};
}

std::vector<double> feSpaceTriRT1::d2Ldrs(double */*r*/)
{
  return {0., 0., 0., 0., 0., 0.};
}

std::vector<double> feSpaceTriRT1::d2Lds2(double */*r*/)
{
  return {0., 0., 0., 0., 0., 0.};
}

void feSpaceTriRT1::initializeNumberingUnknowns()
{
  int nDOFPerEdge = 1;
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, nDOFPerEdge);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 1, nDOFPerEdge);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 2, nDOFPerEdge);
  }
}

void feSpaceTriRT1::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 2);
  }
}

void feSpaceTriRT1::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  // Warning : the standard numbering of the Raviart-Thomas dof
  // is different from the numbering of the mesh edges.
  // For instance : a triangle (v0,v1,v2) has its mesh edges numbered
  // (e0,e1,e2) = (v1-v0, v2-v1, v0-v2), but the first RT dof is 
  // associated to the edge v2-v1, thus e1, and yields the order (e1, e2, e0)
  adr[0] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 0);
  adr[1] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 0);
  adr[2] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 0);
}

// -----------------------------------------------------------------------------
// Lagrange element of degree 2 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
feSpaceTriP2::feSpaceTriP2(std::string cncGeoID)
  : feScalarSpace(2, nullptr, "GEO", cncGeoID, nullptr)
{
  _isGeometricInterpolant = true;
  _nFunctions = 6;
  _Lcoor = {0., 0., 0., 1., 0., 0., 0., 1., 0., 0.5, 0., 0., 0.5, 0.5, 0., 0., 0.5, 0.};
  _dofLocations = {dofLocation::VERTEX, dofLocation::VERTEX, dofLocation::VERTEX,
                   dofLocation::EDGE,   dofLocation::EDGE,   dofLocation::EDGE};
}

feSpaceTriP2::feSpaceTriP2(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                           feFunction *fct, const bool useGlobalShapeFunctions)
  : feScalarSpace(2, mesh, fieldID, cncGeoID, fct, useGlobalShapeFunctions)
{
  _nFunctions = 6;
  _Lcoor = {0., 0., 0., 1., 0., 0., 0., 1., 0., 0.5, 0., 0., 0.5, 0.5, 0., 0., 0.5, 0.};
  _dofLocations = {dofLocation::VERTEX, dofLocation::VERTEX, dofLocation::VERTEX,
                   dofLocation::EDGE,   dofLocation::EDGE,   dofLocation::EDGE};
}

std::vector<double> feSpaceTriP2::L(double *r)
{
  return {(1. - r[0] - r[1]) * (1. - 2. * r[0] - 2. * r[1]),
          r[0] * (2. * r[0] - 1.),
          r[1] * (2. * r[1] - 1.),
          4. * r[0] * (1. - r[0] - r[1]),
          4. * r[0] * r[1],
          4. * r[1] * (1. - r[0] - r[1])};
}
void feSpaceTriP2::L(double *r, double *L)
{
  L[0] = (1. - r[0] - r[1]) * (1. - 2. * r[0] - 2. * r[1]);
  L[1] = r[0] * (2. * r[0] - 1.);
  L[2] = r[1] * (2. * r[1] - 1.);
  L[3] = 4. * r[0] * (1. - r[0] - r[1]);
  L[4] = 4. * r[0] * r[1];
  L[5] = 4. * r[1] * (1. - r[0] - r[1]);
}

/*
  Shape functions in the physical space (x,y). Does not require the reference element.
  The degrees of freedom (functionals of p_i(x,y)) are defined based on the actual vertices
  present in the mesh, i.e. P2 global functions require 6-node triangles. To bypass this,
  one could add 3 fictitious mid-nodes in the geoCoord vector, to define global functions
  on a mesh of 3-node triangles.

  For a better matrix conditioning, the functions are expressed in the local frame centered
  at the element's barycenter (x-xc, y-yc).

  The shape functions and its derivatives are stored in L, dLdx and dLdy.
*/
feStatus feSpaceTriP2::Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                             std::vector<double> &dLdx, std::vector<double> &dLdy)
{
  std::vector<double> geoCoord;
  _mesh->getCoord(_cncGeoTag, iElm, geoCoord);

  EigenMat m = EigenMat::Zero(6, 6);
  int ex[6] = {0, 1, 0, 2, 1, 0}; // Coefficients des monomes
  int ey[6] = {0, 0, 1, 0, 1, 2};
  int dxx[6] = {0, 0, 0, 1, 0, 0}; // Coefficients des derivees en x des monomes
  int dxy[6] = {0, 0, 0, 0, 1, 0};
  int dyx[6] = {0, 0, 0, 0, 1, 0}; // Coefficients des derivees en y des monomes
  int dyy[6] = {0, 0, 0, 0, 0, 1};

  // Coordinates of the barycenter
  std::vector<double> xc(3, 0.0);
  double rc[3] = {1. / 3., 1. / 3., 1. / 3.};
  this->getCncGeo()->getFeSpace()->interpolateVectorField(geoCoord, rc, xc);

  // Polynomial coefficients in the frame centered at the barycenter
  for(int i = 0; i < 6; ++i) {
    for(int j = 0; j < 6; ++j) {
      m(i, j) = pow(geoCoord[i * 3 + 0] - xc[0], ex[j]) * pow(geoCoord[i * 3 + 1] - xc[1], ey[j]);
      // m(i,j) = pow( geoCoord[i*3+0], ex[j]) * pow( geoCoord[i*3+1], ey[j]);
    }
  }
  /* Complete geoCoord if not a 6-node triangle.
  Mid-edge nodes are the average of the P1 nodes. */
  // for(int j = 0; j < 6; ++j) {
  //   m(3,j) = pow( (geoCoord[3*0+0] + geoCoord[3*1+0])/2.0, ex[j]) * pow( (geoCoord[3*0+1] +
  //   geoCoord[3*1+1])/2.0, ey[j]); m(4,j) = pow( (geoCoord[3*1+0] + geoCoord[3*2+0])/2.0, ex[j]) *
  //   pow( (geoCoord[3*1+1] + geoCoord[3*2+1])/2.0, ey[j]); m(5,j) = pow( (geoCoord[3*2+0] +
  //   geoCoord[3*0+0])/2.0, ex[j]) * pow( (geoCoord[3*2+1] + geoCoord[3*0+1])/2.0, ey[j]);
  // }

  m = m.inverse() * I6;

  for(int i = 0; i < 6; ++i) {
    for(int j = 0; j < 6; ++j) {
      if(isnan(m(i, j))) {
        return feErrorMsg(FE_STATUS_ERROR, "NAN");
      }
    }
  }

  for(int i = 0; i < 6; ++i) {
    L[i] = 0.0;
    dLdx[i] = 0.0;
    dLdy[i] = 0.0;
  }

  // x is expressed in the frame centered at the barycenter
  for(int i = 0; i < 6; ++i) {
    for(int j = 0; j < 6; ++j) {
      L[i] += m(j, i) * pow(x[0] - xc[0], ex[j]) * pow(x[1] - xc[1], ey[j]);
      dLdx[i] += ex[j] * m(j, i) * pow(x[0] - xc[0], dxx[j]) * pow(x[1] - xc[1], dxy[j]);
      dLdy[i] += ey[j] * m(j, i) * pow(x[0] - xc[0], dyx[j]) * pow(x[1] - xc[1], dyy[j]);
      //    L[i] +=         m(j, i) * pow(x[0], ex[j])  * pow(x[1], ey[j]);
      // dLdx[i] += ex[j] * m(j, i) * pow(x[0], dxx[j]) * pow(x[1], dxy[j]);
      // dLdy[i] += ey[j] * m(j, i) * pow(x[0], dyx[j]) * pow(x[1], dyy[j]);
    }
  }

  // for(int i = 0; i < 6; ++i) {
  //   std::cout<<L[i]<<std::endl;
  //   if(isnan(L[i]))
  //     feErrorMsg(FE_STATUS_ERROR, "NAN");
  // }

  // std::cout<<"Evaluating at "<<x[0]<<" - "<<x[1]<<" on elem "<<iElm<<std::endl;
  // printf("Barycentre en %f - %f\n", xc[0], xc[1]);
  // double sum = 0.0;
  // double sumx = 0.0;
  // double sumy = 0.0;
  // for(int ii = 0; ii < L.size(); ++ii){
  //   sum += L[ii];
  //   printf("%+-12.12e\n", L[ii]);
  // }
  // std::cout<<std::endl;
  // for(int ii = 0; ii < dLdx.size(); ++ii){
  //   sumx += dLdx[ii];
  //   printf("%+-12.12e\n", dLdx[ii]);
  // }
  // std::cout<<std::endl;
  // for(int ii = 0; ii < dLdy.size(); ++ii){
  //   sumy += dLdy[ii];
  //   printf("%+-12.12e\n", dLdy[ii]);
  // }
  // std::cout<<"Done sum = "<<sum<<std::endl;
  // std::cout<<"Done sum = "<<sumx<<std::endl;
  // std::cout<<"Done sum = "<<sumy<<std::endl;
  return FE_STATUS_OK;

  // Orthogonalisation des monomes d'apres la these/papier de Tesini

  int n = 6;

  std::vector<std::vector<double> > coeffs;
  coeffs.push_back({1, 0, 0, 0, 0, 0});
  coeffs.push_back({0, 1, 0, 0, 0, 0});
  coeffs.push_back({0, 0, 1, 0, 0, 0});
  coeffs.push_back({0, 0, 0, 1, 0, 0});
  coeffs.push_back({0, 0, 0, 0, 1, 0});
  coeffs.push_back({0, 0, 0, 0, 0, 1});

  std::vector<Polynomial> p;
  for(int i = 0; i < n; ++i) {
    p.push_back(Polynomial(2, coeffs[i]));
  }
  printf("Print initial polynomials on elem %d\n", iElm);
  for(auto P : p) P.print();

  std::vector<double> r(n * n, 0.0);

  for(int i = 0; i < n; i++) {
    for(int k = 0; k < 2; k++) {
      for(int j = 0; j < i; j++) {
        r[i * n + j] = p[i].innerProduct(p[j], _mesh, _cncGeoID, iElm);
        Polynomial tmp = p[j];
        tmp *= r[i * n + j];
        p[i] -= tmp;
        p[i].print();
      }
      r[i * n + i] = sqrt(p[i].innerProduct(p[i], _mesh, _cncGeoID, iElm));
      double val = 1.0 / r[i * n + i];
      printf("\n rii = %+-12.12e - inverse = %+-12.12e\n", r[i * n + i], val);
      p[i] *= val;
    }
  }

  // std::vector<double> res(6,0.0);

  // for(int i = 0; i < n; i++){
  //   for(int j = 0; j < i; j++){
  //     Polynomial tmp = p[j];
  //     tmp *= r[i*n+j];
  //     p[i] -= tmp;
  //   }
  //   double val = 1.0/r[i*n+i];
  //   p[i] *= val;
  // }

  for(auto P : p) P.print();

  std::cout << std::endl;

  for(int i = 0; i < 6; ++i) {
    for(int j = 0; j < 6; ++j) {
      std::cout << p[i].innerProduct(p[j], _mesh, _cncGeoID, iElm) << std::endl;
    }
    std::cout << std::endl;
  }

  exit(-1);

  return FE_STATUS_OK;

  // return {p[0].eval(x), p[1].eval(x), p[2].eval(x), p[3].eval(x), p[4].eval(x), p[5].eval(x)};
}

std::vector<double> feSpaceTriP2::dLdr(double *r)
{
  return {4. * (r[0] + r[1]) - 3.,      4. * r[0] - 1., 0.,
          4. * (1. - 2. * r[0] - r[1]), 4. * r[1],      -4. * r[1]};
}

std::vector<double> feSpaceTriP2::dLds(double *r)
{
  return {4. * (r[0] + r[1]) - 3.,     0., 4. * r[1] - 1., -4. * r[0], 4. * r[0],
          4. * (1. - r[0] - 2. * r[1])};
}

std::vector<double> feSpaceTriP2::d2Ldr2(double *r) { UNUSED(r); return {4., 4., 0., -8., 0., 0.}; }
std::vector<double> feSpaceTriP2::d2Ldrs(double *r) { UNUSED(r); return {4., 0., 0., -4., 4., -4.}; }
std::vector<double> feSpaceTriP2::d2Lds2(double *r) { UNUSED(r); return {4., 0., 4., 0., 0., -8.}; }

void feSpaceTriP2::initializeNumberingUnknowns()
{
  int nDOFPerEdge = 1;
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    for(int j = 0; j < this->getCncGeo()->getNumVerticesPerElem(); ++j) {
      _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, j);
    }
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, nDOFPerEdge);
      _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 1, nDOFPerEdge);
      _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 2, nDOFPerEdge);
    }
  }
}

void feSpaceTriP2::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    for(int j = 0; j < this->getCncGeo()->getNumVerticesPerElem(); ++j) {
      _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, j);
    }
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0);
      _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 1);
      _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 2);
    }
  }
}

void feSpaceTriP2::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2);
  if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() == 2) {
    adr[3] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 3);
    adr[4] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 4);
    adr[5] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 5);
  } else {
    adr[3] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 0);
    adr[4] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 0);
    adr[5] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 0);
  }
}

// -----------------------------------------------------------------------------
// Lagrange element of degree 2 enriched with cubic bubble on reference triangle
// -----------------------------------------------------------------------------
feSpaceTriP2Bubble::feSpaceTriP2Bubble(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                           feFunction *fct, const bool useGlobalShapeFunctions)
  : feScalarSpace(2, mesh, fieldID, cncGeoID, fct, useGlobalShapeFunctions)
{
  _nFunctions = 7;
  _Lcoor = {0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.5, 0.0, 0.0,
            0.5, 0.5, 0.0,
            0.0, 0.5, 0.0,
            1./3., 1./3., 0.0};
  _dofLocations = {dofLocation::VERTEX, dofLocation::VERTEX, dofLocation::VERTEX,
                   dofLocation::EDGE,   dofLocation::EDGE,   dofLocation::EDGE,
                   dofLocation::ELEMENT};
}

std::vector<double> feSpaceTriP2Bubble::L(double *r)
{
  return {(1. - r[0] - r[1]) * (1. - 2. * r[0] - 2. * r[1]),
          r[0] * (2. * r[0] - 1.),
          r[1] * (2. * r[1] - 1.),
          4. * r[0] * (1. - r[0] - r[1]),
          4. * r[0] * r[1],
          4. * r[1] * (1. - r[0] - r[1]),
          27. * r[0] * r[1] * (1. - r[0] - r[1])};
}
void feSpaceTriP2Bubble::L(double *r, double *L)
{
  L[0] = (1. - r[0] - r[1]) * (1. - 2. * r[0] - 2. * r[1]);
  L[1] = r[0] * (2. * r[0] - 1.);
  L[2] = r[1] * (2. * r[1] - 1.);
  L[3] = 4. * r[0] * (1. - r[0] - r[1]);
  L[4] = 4. * r[0] * r[1];
  L[5] = 4. * r[1] * (1. - r[0] - r[1]);
  L[6] = 27. * r[0] * r[1] * (1. - r[0] - r[1]);
}

std::vector<double> feSpaceTriP2Bubble::dLdr(double *r)
{
  return {4. * (r[0] + r[1]) - 3.,      4. * r[0] - 1., 0.,
          4. * (1. - 2. * r[0] - r[1]), 4. * r[1], -4. * r[1],
          27. * r[1] * (1. - r[1] - 2.*r[0])};
}

std::vector<double> feSpaceTriP2Bubble::dLds(double *r)
{
  return {4. * (r[0] + r[1]) - 3.,     0., 4. * r[1] - 1., -4. * r[0], 4. * r[0],
          4. * (1. - r[0] - 2. * r[1]),
          27. * r[0] * (1. - r[0] - 2.*r[1])};
}

std::vector<double> feSpaceTriP2Bubble::d2Ldr2(double *r) { return {4., 4., 0., -8., 0.,  0., -54.*r[1]}; }
std::vector<double> feSpaceTriP2Bubble::d2Ldrs(double *r) { return {4., 0., 0., -4., 4., -4., 27. * (1. - 2.*(r[0]+r[1]))}; }
std::vector<double> feSpaceTriP2Bubble::d2Lds2(double *r) { return {4., 0., 4.,  0., 0., -8., -54.*r[0]}; }

void feSpaceTriP2Bubble::initializeNumberingUnknowns()
{
  int nDOFPerEdge = 1;
  int nDOFPerElement = 1;
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    for(int j = 0; j < this->getCncGeo()->getNumVerticesPerElem(); ++j) {
      _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, j);
    }
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, nDOFPerEdge);
      _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 1, nDOFPerEdge);
      _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 2, nDOFPerEdge);
    }
    _numbering->setUnknownElementDOF(_mesh, _cncGeoID, i, nDOFPerElement);
  }
}

void feSpaceTriP2Bubble::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    for(int j = 0; j < this->getCncGeo()->getNumVerticesPerElem(); ++j) {
      _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, j);
    }
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0);
      _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 1);
      _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 2);
    }
    _numbering->setEssentialElementDOF(_mesh, _cncGeoID, i);
  }
}

void feSpaceTriP2Bubble::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2);
  if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() == 2) {
    adr[3] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 3);
    adr[4] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 4);
    adr[5] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 5);
  } else {
    adr[3] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 0);
    adr[4] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 0);
    adr[5] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 0);
  }
  adr[6] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
}

// -----------------------------------------------------------------------------
// Vector Lagrange element of degree 2 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
template <int dim>
feSpaceVecTriP2<dim>::feSpaceVecTriP2(feMesh *mesh, const std::string fieldID,
                                      const std::string cncGeoID, feVectorFunction *fct,
                                      const bool useGlobalShapeFunctions)
  : feVectorSpace(2, mesh, fieldID, cncGeoID, fct, useGlobalShapeFunctions)
{
  _nComponents = dim;
  _nFunctions = 6 * dim;
  double coor[6 * 3] = {0., 0., 0., 1., 0., 0., 0., 1., 0., 0.5, 0., 0., 0.5, 0.5, 0., 0., 0.5, 0.};
  _Lcoor.resize(6 * 3 * dim);
  duplicateScalarArray(6 * 3, 3, coor, dim, _Lcoor);
  _dofLocations.resize(_nFunctions);
  for(int i = 0; i < _nFunctions; ++i) {
    _dofLocations[i] = (i < 3 * dim) ? dofLocation::VERTEX : dofLocation::EDGE;
  }
}

template <int dim> std::vector<double> feSpaceVecTriP2<dim>::L(double *r)
{
  double phi[6] = {(1. - r[0] - r[1]) * (1. - 2. * r[0] - 2. * r[1]),
                   r[0] * (2. * r[0] - 1.),
                   r[1] * (2. * r[1] - 1.),
                   4. * r[0] * (1. - r[0] - r[1]),
                   4. * r[0] * r[1],
                   4. * r[1] * (1. - r[0] - r[1])};
  std::vector<double> res(_nFunctions*_nComponents);
  // duplicateScalarArray(6, 1, phi, dim, res);
  vectorLagrangeLayout(phi, 6, _nComponents, res.data());
  return res;
}

template <int dim> void feSpaceVecTriP2<dim>::L(double *r, double *res)
{
  double phi[6] = {(1. - r[0] - r[1]) * (1. - 2. * r[0] - 2. * r[1]),
                   r[0] * (2. * r[0] - 1.),
                   r[1] * (2. * r[1] - 1.),
                   4. * r[0] * (1. - r[0] - r[1]),
                   4. * r[0] * r[1],
                   4. * r[1] * (1. - r[0] - r[1])};
  // duplicateScalarArray(6, 1, phi, dim, res);
  vectorLagrangeLayout(phi, 6, _nComponents, res);
}

template <int dim> std::vector<double> feSpaceVecTriP2<dim>::dLdr(double *r)
{
  double dldr[6] = {4. * (r[0] + r[1]) - 3.,      4. * r[0] - 1., 0.,
                    4. * (1. - 2. * r[0] - r[1]), 4. * r[1],      -4. * r[1]};
  std::vector<double> res(_nFunctions*_nComponents);
  // duplicateScalarArray(6, 1, dldr, dim, res);
  vectorLagrangeLayout(dldr, 6, _nComponents, res.data());
  return res;
}

template <int dim> std::vector<double> feSpaceVecTriP2<dim>::dLds(double *r)
{
  double dlds[6] = {4. * (r[0] + r[1]) - 3.,     0., 4. * r[1] - 1., -4. * r[0], 4. * r[0],
                    4. * (1. - r[0] - 2. * r[1])};
  std::vector<double> res(_nFunctions*_nComponents);
  // duplicateScalarArray(6, 1, dlds, dim, res);
  vectorLagrangeLayout(dlds, 6, _nComponents, res.data());
  return res;
}

template <int dim> std::vector<double> feSpaceVecTriP2<dim>::d2Ldr2(double */*r*/)
{
  double d2ldr2[6] = {4., 4., 0., -8., 0., 0.};
  std::vector<double> res(_nFunctions*_nComponents);
  // duplicateScalarArray(6, 1, d2ldr2, dim, res);
  vectorLagrangeLayout(d2ldr2, 6, _nComponents, res.data());
  return res;
}

template <int dim> std::vector<double> feSpaceVecTriP2<dim>::d2Ldrs(double */*r*/)
{
  double d2ldrs[6] = {4., 0., 0., -4., 4., -4.};
  std::vector<double> res(_nFunctions*_nComponents);
  // duplicateScalarArray(6, 1, d2ldrs, dim, res);
  vectorLagrangeLayout(d2ldrs, 6, _nComponents, res.data());
  return res;
}

template <int dim> std::vector<double> feSpaceVecTriP2<dim>::d2Lds2(double */*r*/)
{
  double d2lds2[6] = {4., 0., 4., 0., 0., -8.};
  std::vector<double> res(_nFunctions*_nComponents);
  // duplicateScalarArray(6, 1, d2lds2, dim, res);
  vectorLagrangeLayout(d2lds2, 6, _nComponents, res.data());
  return res;
}

template <int dim> void feSpaceVecTriP2<dim>::initializeNumberingUnknowns()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    for(int j = 0; j < this->getCncGeo()->getNumVerticesPerElem(); ++j) {
      _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, j, dim);
    }
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, dim);
      _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 1, dim);
      _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 2, dim);
    }
  }
}

template <int dim> void feSpaceVecTriP2<dim>::initializeNumberingEssential()
{
  // for(int iComp = 0; iComp < _nComponents; ++iComp){
  //   if(_essentialComponents[iComp]){
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    for(int j = 0; j < this->getCncGeo()->getNumVerticesPerElem(); ++j) {
      _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, j); //, iComp);
    }
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0); //, iComp);
      _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 1); //, iComp);
      _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 2); //, iComp);
    }
  }
  //   }
  // }
}

template <int dim>
void feSpaceVecTriP2<dim>::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 1);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 0);
  adr[3] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 1);
  adr[4] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2, 0);
  adr[5] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2, 1);
  if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() == 2) {
    adr[6] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 3, 0);
    adr[7] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 3, 1);
    adr[8] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 4, 0);
    adr[9] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 4, 1);
    adr[10] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 5, 0);
    adr[11] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 5, 1);
  } else {
    adr[6] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 0);
    adr[7] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 1);
    adr[8] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 0);
    adr[9] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 1);
    adr[10] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 0);
    adr[11] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 1);
  }
}

template class feSpaceVecTriP2<1>;
template class feSpaceVecTriP2<2>;

// -----------------------------------------------------------------------------
// Vector Lagrange element of degree 2 with cubic bubble on reference triangle
// -----------------------------------------------------------------------------
template <int dim>
feSpaceVecTriP2Bubble<dim>::feSpaceVecTriP2Bubble(feMesh *mesh, const std::string fieldID,
                                      const std::string cncGeoID, feVectorFunction *fct,
                                      const bool useGlobalShapeFunctions)
  : feVectorSpace(2, mesh, fieldID, cncGeoID, fct, useGlobalShapeFunctions)
{
  _nComponents = dim;
  _nFunctions = 7 * dim;
  double coor[7 * 3] = {0., 0., 0.,
    1., 0., 0.,
    0., 1., 0.,
    0.5, 0., 0.,
    0.5, 0.5, 0.,
    0., 0.5, 0.,
    1./3., 1./3., 0.};
  _Lcoor.resize(7 * 3 * dim);
  duplicateScalarArray(7 * 3, 3, coor, dim, _Lcoor);
  _dofLocations.resize(_nFunctions);
  for(int i = 0; i < _nFunctions; ++i) {
    if     (i < 3 * dim) _dofLocations[i] = dofLocation::VERTEX;
    else if(i < 6 * dim) _dofLocations[i] = dofLocation::EDGE;
    else _dofLocations[i] = dofLocation::ELEMENT;
  }
}

template <int dim> std::vector<double> feSpaceVecTriP2Bubble<dim>::L(double *r)
{
  double phi[7] = {(1. - r[0] - r[1]) * (1. - 2. * r[0] - 2. * r[1]),
                   r[0] * (2. * r[0] - 1.),
                   r[1] * (2. * r[1] - 1.),
                   4. * r[0] * (1. - r[0] - r[1]),
                   4. * r[0] * r[1],
                   4. * r[1] * (1. - r[0] - r[1]),
                   27. * r[0] * r[1] * (1. - r[0] - r[1])};
  std::vector<double> res(_nFunctions);
  duplicateScalarArray(7, 1, phi, dim, res);
  return res;
}

template <int dim> void feSpaceVecTriP2Bubble<dim>::L(double *r, double *res)
{
  double phi[7] = {(1. - r[0] - r[1]) * (1. - 2. * r[0] - 2. * r[1]),
                   r[0] * (2. * r[0] - 1.),
                   r[1] * (2. * r[1] - 1.),
                   4. * r[0] * (1. - r[0] - r[1]),
                   4. * r[0] * r[1],
                   4. * r[1] * (1. - r[0] - r[1]),
                   27. * r[0] * r[1] * (1. - r[0] - r[1])};
  duplicateScalarArray(7, 1, phi, dim, res);
}

template <int dim> std::vector<double> feSpaceVecTriP2Bubble<dim>::dLdr(double *r)
{
  double dldr[7] = {4. * (r[0] + r[1]) - 3.,      4. * r[0] - 1., 0.,
                    4. * (1. - 2. * r[0] - r[1]), 4. * r[1],      -4. * r[1],
                    27. * r[1] * (1. - r[1] - 2.*r[0])};
  std::vector<double> res(_nFunctions);
  duplicateScalarArray(7, 1, dldr, dim, res);
  return res;
}

template <int dim> std::vector<double> feSpaceVecTriP2Bubble<dim>::dLds(double *r)
{
  double dlds[7] = {4. * (r[0] + r[1]) - 3.,     0., 4. * r[1] - 1., -4. * r[0], 4. * r[0],
                    4. * (1. - r[0] - 2. * r[1]),
                    27. * r[0] * (1. - r[0] - 2.*r[1])};
  std::vector<double> res(_nFunctions);
  duplicateScalarArray(7, 1, dlds, dim, res);
  return res;
}

template <int dim> std::vector<double> feSpaceVecTriP2Bubble<dim>::d2Ldr2(double *r)
{
  UNUSED(r);
  double d2ldr2[7] = {4., 4., 0., -8., 0., 0., -54.*r[1]};
  std::vector<double> res(_nFunctions);
  duplicateScalarArray(7, 1, d2ldr2, dim, res);
  return res;
}

template <int dim> std::vector<double> feSpaceVecTriP2Bubble<dim>::d2Ldrs(double *r)
{
  UNUSED(r);
  double d2ldrs[7] = {4., 0., 0., -4., 4., -4., 27. * (1. - 2.*(r[0]+r[1]))};
  std::vector<double> res(_nFunctions);
  duplicateScalarArray(7, 1, d2ldrs, dim, res);
  return res;
}

template <int dim> std::vector<double> feSpaceVecTriP2Bubble<dim>::d2Lds2(double *r)
{
  UNUSED(r);
  double d2lds2[7] = {4., 0., 4., 0., 0., -8., -54.*r[0]};
  std::vector<double> res(_nFunctions);
  duplicateScalarArray(7, 1, d2lds2, dim, res);
  return res;
}

template <int dim> void feSpaceVecTriP2Bubble<dim>::initializeNumberingUnknowns()
{
  int nDOFPerVertex = dim;
  int nDOFPerEdge = dim;
  int nDOFPerElem = dim;
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    for(int j = 0; j < this->getCncGeo()->getNumVerticesPerElem(); ++j) {
      _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, j, nDOFPerVertex);
    }
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, nDOFPerEdge);
      _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 1, nDOFPerEdge);
      _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 2, nDOFPerEdge);
    }
    _numbering->setUnknownElementDOF(_mesh, _cncGeoID, i, nDOFPerElem);
  }
}

template <int dim> void feSpaceVecTriP2Bubble<dim>::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    for(int j = 0; j < this->getCncGeo()->getNumVerticesPerElem(); ++j) {
      _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, j);
    }
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0);
      _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 1);
      _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 2);
    }
    _numbering->setEssentialElementDOF(_mesh, _cncGeoID, i);
  }
}

template <int dim>
void feSpaceVecTriP2Bubble<dim>::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 1);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 0);
  adr[3] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 1);
  adr[4] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2, 0);
  adr[5] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2, 1);
  if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() == 2) {
    adr[6] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 3, 0);
    adr[7] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 3, 1);
    adr[8] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 4, 0);
    adr[9] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 4, 1);
    adr[10] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 5, 0);
    adr[11] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 5, 1);
  } else {
    adr[6] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 0);
    adr[7] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 1);
    adr[8] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 0);
    adr[9] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 1);
    adr[10] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 0);
    adr[11] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 1);
  }
  adr[12] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
  adr[13] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 1);
}

template class feSpaceVecTriP2Bubble<1>;
template class feSpaceVecTriP2Bubble<2>;

// -----------------------------------------------------------------------------
// Non-conforming Crouzeix-Raviart element of degree 2 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
feSpaceTri_CR2::feSpaceTri_CR2(feMesh *mesh, std::string fieldID, std::string cncGeoID,
                               feFunction *fct)
  : feScalarSpace(2, mesh, fieldID, cncGeoID, fct)
{
  _nFunctions = 7;
  _Lcoor = {0., 0.,  0.,  1., 0., 0.,  0., 1.,    0.,    0.5, 0.,
            0., 0.5, 0.5, 0., 0., 0.5, 0., 1 / 3, 1 / 3, 0};
  _dofLocations = {dofLocation::EDGE, dofLocation::EDGE, dofLocation::EDGE,
                   dofLocation::EDGE, dofLocation::EDGE, dofLocation::EDGE};
}

std::vector<double> feSpaceTri_CR2::L(double *r)
{
  return {(1. - r[0] - r[1]) * (1. - 2. * r[0] - 2. * r[1]),
          r[0] * (2. * r[0] - 1.),
          r[1] * (2. * r[1] - 1.),
          4. * r[0] * (1. - r[0] - r[1]),
          4. * r[0] * r[1],
          4. * r[1] * (1. - r[0] - r[1]),
          2. - 3. * (r[0] * r[0] + r[1] * r[1] + (1. - r[0] - r[1]) * (1. - r[0] - r[1]))};
}
void feSpaceTri_CR2::L(double *r, double *L)
{
  L[0] = (1. - r[0] - r[1]) * (1. - 2. * r[0] - 2. * r[1]);
  L[1] = r[0] * (2. * r[0] - 1.);
  L[2] = r[1] * (2. * r[1] - 1.);
  L[3] = 4. * r[0] * (1. - r[0] - r[1]);
  L[4] = 4. * r[0] * r[1];
  L[5] = 4. * r[1] * (1. - r[0] - r[1]);
  L[6] = 2. - 3. * (r[0] * r[0] + r[1] * r[1] + (1. - r[0] - r[1]) * (1. - r[0] - r[1]));
}

std::vector<double> feSpaceTri_CR2::dLdr(double *r)
{
  return {4. * (r[0] + r[1]) - 3.,      4. * r[0] - 1., 0.,
          4. * (1. - 2. * r[0] - r[1]), 4. * r[1],      -4. * r[1],
          6. * (1. - r[1] - 2. * r[0])};
}

std::vector<double> feSpaceTri_CR2::dLds(double *r)
{
  return {4. * (r[0] + r[1]) - 3.,
          0.,
          4. * r[1] - 1.,
          -4. * r[0],
          4. * r[0],
          4. * (1. - r[0] - 2. * r[1]),
          6. * (1. - r[0] - 2. * r[1])};
}

void feSpaceTri_CR2::initializeNumberingUnknowns()
{
  int nDOFPerEdge = 1;
  int nDOFPerElem = 1;
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    for(int j = 0; j < this->getCncGeo()->getNumVerticesPerElem(); ++j) {
      _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, j);
    }
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, nDOFPerEdge);
      _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 1, nDOFPerEdge);
      _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 2, nDOFPerEdge);
      _numbering->setUnknownElementDOF(_mesh, _cncGeoID, i, nDOFPerElem);
    }
  }
}

void feSpaceTri_CR2::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    for(int j = 0; j < this->getCncGeo()->getNumVerticesPerElem(); ++j) {
      _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, j);
    }
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0);
      _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 1);
      _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 2);
      _numbering->setEssentialElementDOF(_mesh, _cncGeoID, i);
    }
  }
}

void feSpaceTri_CR2::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2);
  if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() == 2) {
    adr[3] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 3);
    adr[4] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 4);
    adr[5] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 5);
  } else {
    adr[3] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 0);
    adr[4] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 0);
    adr[5] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 0);
  }
  adr[6] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
}

// -----------------------------------------------------------------------------
// Lagrange element of degree 3 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
feSpaceTriP3::feSpaceTriP3(std::string cncGeoID)
  : feScalarSpace(2, nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 10;
  _Lcoor = {0., 0., 0.,      1., 0., 0.,      0.,      1.,      0.,      1. / 3.,
            0., 0., 2. / 3., 0., 0., 2. / 3., 1. / 3., 0.,      1. / 3., 2. / 3.,
            0., 0., 2. / 3., 0., 0., 1. / 3., 0.,      1. / 3., 1. / 3., 0.};
  _isGeometricInterpolant = true;
  _dofLocations = {dofLocation::VERTEX, dofLocation::VERTEX, dofLocation::VERTEX,
                   dofLocation::EDGE, dofLocation::EDGE, dofLocation::EDGE,
                   dofLocation::EDGE, dofLocation::EDGE, dofLocation::EDGE,
                   dofLocation::ELEMENT};
}

feSpaceTriP3::feSpaceTriP3(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                           feFunction *fct, const bool useGlobalShapeFunctions)
  : feScalarSpace(2, mesh, fieldID, cncGeoID, fct, useGlobalShapeFunctions)
{
  _nFunctions = 10;
  _Lcoor = {0., 0., 0.,      1., 0., 0.,      0.,      1.,      0.,      1. / 3.,
            0., 0., 2. / 3., 0., 0., 2. / 3., 1. / 3., 0.,      1. / 3., 2. / 3.,
            0., 0., 2. / 3., 0., 0., 1. / 3., 0.,      1. / 3., 1. / 3., 0.};
  _dofLocations = {dofLocation::VERTEX, dofLocation::VERTEX, dofLocation::VERTEX,
                   dofLocation::EDGE, dofLocation::EDGE, dofLocation::EDGE,
                   dofLocation::EDGE, dofLocation::EDGE, dofLocation::EDGE,
                   dofLocation::ELEMENT};
}

std::vector<double> feSpaceTriP3::L(double *r)
{
  double R = r[0];
  double S = r[1];
  return {R * (-1.1E1 / 2.0) - S * (1.1E1 / 2.0) + R * S * 1.8E1 - R * (S * S) * (2.7E1 / 2.0) -
            (R * R) * S * (2.7E1 / 2.0) + (R * R) * 9.0 - (R * R * R) * (9.0 / 2.0) +
            (S * S) * 9.0 - (S * S * S) * (9.0 / 2.0) + 1.0,
          R - (R * R) * (9.0 / 2.0) + (R * R * R) * (9.0 / 2.0),
          S - (S * S) * (9.0 / 2.0) + (S * S * S) * (9.0 / 2.0),
          R * 9.0 - R * S * (4.5E1 / 2.0) + R * (S * S) * (2.7E1 / 2.0) + (R * R) * S * 2.7E1 -
            (R * R) * (4.5E1 / 2.0) + (R * R * R) * (2.7E1 / 2.0),
          R * (-9.0 / 2.0) + R * S * (9.0 / 2.0) - (R * R) * S * (2.7E1 / 2.0) + (R * R) * 1.8E1 -
            (R * R * R) * (2.7E1 / 2.0),
          R * S * (-9.0 / 2.0) + (R * R) * S * (2.7E1 / 2.0),
          R * S * (-9.0 / 2.0) + R * (S * S) * (2.7E1 / 2.0),
          S * (-9.0 / 2.0) + R * S * (9.0 / 2.0) - R * (S * S) * (2.7E1 / 2.0) + (S * S) * 1.8E1 -
            (S * S * S) * (2.7E1 / 2.0),
          S * 9.0 - R * S * (4.5E1 / 2.0) + R * (S * S) * 2.7E1 + (R * R) * S * (2.7E1 / 2.0) -
            (S * S) * (4.5E1 / 2.0) + (S * S * S) * (2.7E1 / 2.0),
          R * S * 2.7E1 - R * (S * S) * 2.7E1 - (R * R) * S * 2.7E1};
}
void feSpaceTriP3::L(double *r, double *L)
{
  double R = r[0];
  double S = r[1];
  L[0] = R * (-1.1E1 / 2.0) - S * (1.1E1 / 2.0) + R * S * 1.8E1 - R * (S * S) * (2.7E1 / 2.0) -
         (R * R) * S * (2.7E1 / 2.0) + (R * R) * 9.0 - (R * R * R) * (9.0 / 2.0) + (S * S) * 9.0 -
         (S * S * S) * (9.0 / 2.0) + 1.0;
  L[1] = R - (R * R) * (9.0 / 2.0) + (R * R * R) * (9.0 / 2.0);
  L[2] = S - (S * S) * (9.0 / 2.0) + (S * S * S) * (9.0 / 2.0);
  L[3] = R * 9.0 - R * S * (4.5E1 / 2.0) + R * (S * S) * (2.7E1 / 2.0) + (R * R) * S * 2.7E1 -
         (R * R) * (4.5E1 / 2.0) + (R * R * R) * (2.7E1 / 2.0);
  L[4] = R * (-9.0 / 2.0) + R * S * (9.0 / 2.0) - (R * R) * S * (2.7E1 / 2.0) + (R * R) * 1.8E1 -
         (R * R * R) * (2.7E1 / 2.0);
  L[5] = R * S * (-9.0 / 2.0) + (R * R) * S * (2.7E1 / 2.0);
  L[6] = R * S * (-9.0 / 2.0) + R * (S * S) * (2.7E1 / 2.0);
  L[7] = S * (-9.0 / 2.0) + R * S * (9.0 / 2.0) - R * (S * S) * (2.7E1 / 2.0) + (S * S) * 1.8E1 -
         (S * S * S) * (2.7E1 / 2.0);
  L[8] = S * 9.0 - R * S * (4.5E1 / 2.0) + R * (S * S) * 2.7E1 + (R * R) * S * (2.7E1 / 2.0) -
         (S * S) * (4.5E1 / 2.0) + (S * S * S) * (2.7E1 / 2.0);
  L[9] = R * S * 2.7E1 - R * (S * S) * 2.7E1 - (R * R) * S * 2.7E1;
}

std::vector<double> feSpaceTriP3::dLdr(double *r)
{
  double R = r[0];
  double S = r[1];
  return {R * 1.8E1 + S * 1.8E1 - R * S * 2.7E1 - (R * R) * (2.7E1 / 2.0) -
            (S * S) * (2.7E1 / 2.0) - 1.1E1 / 2.0,
          R * -9.0 + (R * R) * (2.7E1 / 2.0) + 1.0,
          0.,
          R * -4.5E1 - S * (4.5E1 / 2.0) + R * S * 5.4E1 + (R * R) * (8.1E1 / 2.0) +
            (S * S) * (2.7E1 / 2.0) + 9.0,
          R * 3.6E1 + S * (9.0 / 2.0) - R * S * 2.7E1 - (R * R) * (8.1E1 / 2.0) - 9.0 / 2.0,
          S * (-9.0 / 2.0) + R * S * 2.7E1,
          S * (-9.0 / 2.0) + (S * S) * (2.7E1 / 2.0),
          S * (9.0 / 2.0) - (S * S) * (2.7E1 / 2.0),
          S * (-4.5E1 / 2.0) + R * S * 2.7E1 + (S * S) * 2.7E1,
          S * 2.7E1 - R * S * 5.4E1 - (S * S) * 2.7E1};
}

std::vector<double> feSpaceTriP3::dLds(double *r)
{
  double R = r[0];
  double S = r[1];
  return {R * 1.8E1 + S * 1.8E1 - R * S * 2.7E1 - (R * R) * (2.7E1 / 2.0) -
            (S * S) * (2.7E1 / 2.0) - 1.1E1 / 2.0,
          0.,
          S * -9.0 + (S * S) * (2.7E1 / 2.0) + 1.0,
          R * (-4.5E1 / 2.0) + R * S * 2.7E1 + (R * R) * 2.7E1,
          R * (9.0 / 2.0) - (R * R) * (2.7E1 / 2.0),
          R * (-9.0 / 2.0) + (R * R) * (2.7E1 / 2.0),
          R * (-9.0 / 2.0) + R * S * 2.7E1,
          R * (9.0 / 2.0) + S * 3.6E1 - R * S * 2.7E1 - (S * S) * (8.1E1 / 2.0) - 9.0 / 2.0,
          R * (-4.5E1 / 2.0) - S * 4.5E1 + R * S * 5.4E1 + (R * R) * (2.7E1 / 2.0) +
            (S * S) * (8.1E1 / 2.0) + 9.0,
          R * 2.7E1 - R * S * 5.4E1 - (R * R) * 2.7E1};
}

void feSpaceTriP3::initializeNumberingUnknowns()
{
  int nDOFPerEdge = 2;
  int nDOFPerElem = 1;
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 2);
    _numbering->setUnknownElementDOF(_mesh, _cncGeoID, i, nDOFPerElem);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, nDOFPerEdge);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 1, nDOFPerEdge);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 2, nDOFPerEdge);
  }
}

void feSpaceTriP3::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 2);
    _numbering->setEssentialElementDOF(_mesh, _cncGeoID, i);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 2);
  }
}

void feSpaceTriP3::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2);
  int e0 = _mesh->getEdge(_cncGeoID, numElem, 0);
  if(e0 > 0) {
    adr[3] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 0);
    adr[4] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 1);
  } else {
    adr[3] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 1);
    adr[4] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 0);
  }
  int e1 = _mesh->getEdge(_cncGeoID, numElem, 1);
  if(e1 > 0) {
    adr[5] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 0);
    adr[6] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 1);
  } else {
    adr[5] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 1);
    adr[6] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 0);
  }
  int e2 = _mesh->getEdge(_cncGeoID, numElem, 2);
  if(e2 > 0) {
    adr[7] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 0);
    adr[8] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 1);
  } else {
    adr[7] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 1);
    adr[8] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 0);
  }
  adr[9] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
}

// -----------------------------------------------------------------------------
// Vector Lagrange element of degree 3 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
template <int dim>
feSpaceVecTriP3<dim>::feSpaceVecTriP3(feMesh *mesh, const std::string fieldID,
                                      const std::string cncGeoID, feVectorFunction *fct,
                                      const bool useGlobalShapeFunctions)
  : feVectorSpace(2, mesh, fieldID, cncGeoID, fct, useGlobalShapeFunctions)
{
  _nComponents = dim;
  _nFunctions = 10 * dim;
  double coor[10 * 3] = {0., 0., 0.,      1., 0., 0.,      0.,      1.,      0.,      1. / 3.,
                         0., 0., 2. / 3., 0., 0., 2. / 3., 1. / 3., 0.,      1. / 3., 2. / 3.,
                         0., 0., 2. / 3., 0., 0., 1. / 3., 0.,      1. / 3., 1. / 3., 0.};
  _Lcoor.resize(10 * 3 * dim);
  duplicateScalarArray(10 * 3, 3, coor, dim, _Lcoor);
  _dofLocations.resize(_nFunctions);
  for(int i = 0; i < _nFunctions; ++i) {
    if(i < 3 * dim)
      _dofLocations[i] = dofLocation::VERTEX;
    else if(i < 9 * dim)
      _dofLocations[i] = dofLocation::EDGE;
    else
      _dofLocations[i] = dofLocation::ELEMENT;
  }
}

template <int dim> std::vector<double> feSpaceVecTriP3<dim>::L(double *r)
{
  double R = r[0];
  double S = r[1];
  double phi[10] = {R * (-1.1E1 / 2.0) - S * (1.1E1 / 2.0) + R * S * 1.8E1 -
                      R * (S * S) * (2.7E1 / 2.0) - (R * R) * S * (2.7E1 / 2.0) + (R * R) * 9.0 -
                      (R * R * R) * (9.0 / 2.0) + (S * S) * 9.0 - (S * S * S) * (9.0 / 2.0) + 1.0,
                    R - (R * R) * (9.0 / 2.0) + (R * R * R) * (9.0 / 2.0),
                    S - (S * S) * (9.0 / 2.0) + (S * S * S) * (9.0 / 2.0),
                    R * 9.0 - R * S * (4.5E1 / 2.0) + R * (S * S) * (2.7E1 / 2.0) +
                      (R * R) * S * 2.7E1 - (R * R) * (4.5E1 / 2.0) + (R * R * R) * (2.7E1 / 2.0),
                    R * (-9.0 / 2.0) + R * S * (9.0 / 2.0) - (R * R) * S * (2.7E1 / 2.0) +
                      (R * R) * 1.8E1 - (R * R * R) * (2.7E1 / 2.0),
                    R * S * (-9.0 / 2.0) + (R * R) * S * (2.7E1 / 2.0),
                    R * S * (-9.0 / 2.0) + R * (S * S) * (2.7E1 / 2.0),
                    S * (-9.0 / 2.0) + R * S * (9.0 / 2.0) - R * (S * S) * (2.7E1 / 2.0) +
                      (S * S) * 1.8E1 - (S * S * S) * (2.7E1 / 2.0),
                    S * 9.0 - R * S * (4.5E1 / 2.0) + R * (S * S) * 2.7E1 +
                      (R * R) * S * (2.7E1 / 2.0) - (S * S) * (4.5E1 / 2.0) +
                      (S * S * S) * (2.7E1 / 2.0),
                    R * S * 2.7E1 - R * (S * S) * 2.7E1 - (R * R) * S * 2.7E1};
  std::vector<double> res(_nFunctions*_nComponents);
  // duplicateScalarArray(10, 1, phi, dim, res);
  vectorLagrangeLayout(phi, 10, _nComponents, res.data());
  return res;
}

template <int dim> void feSpaceVecTriP3<dim>::L(double *r, double *res)
{
  double R = r[0];
  double S = r[1];
  double phi[10] = {R * (-1.1E1 / 2.0) - S * (1.1E1 / 2.0) + R * S * 1.8E1 -
                      R * (S * S) * (2.7E1 / 2.0) - (R * R) * S * (2.7E1 / 2.0) + (R * R) * 9.0 -
                      (R * R * R) * (9.0 / 2.0) + (S * S) * 9.0 - (S * S * S) * (9.0 / 2.0) + 1.0,
                    R - (R * R) * (9.0 / 2.0) + (R * R * R) * (9.0 / 2.0),
                    S - (S * S) * (9.0 / 2.0) + (S * S * S) * (9.0 / 2.0),
                    R * 9.0 - R * S * (4.5E1 / 2.0) + R * (S * S) * (2.7E1 / 2.0) +
                      (R * R) * S * 2.7E1 - (R * R) * (4.5E1 / 2.0) + (R * R * R) * (2.7E1 / 2.0),
                    R * (-9.0 / 2.0) + R * S * (9.0 / 2.0) - (R * R) * S * (2.7E1 / 2.0) +
                      (R * R) * 1.8E1 - (R * R * R) * (2.7E1 / 2.0),
                    R * S * (-9.0 / 2.0) + (R * R) * S * (2.7E1 / 2.0),
                    R * S * (-9.0 / 2.0) + R * (S * S) * (2.7E1 / 2.0),
                    S * (-9.0 / 2.0) + R * S * (9.0 / 2.0) - R * (S * S) * (2.7E1 / 2.0) +
                      (S * S) * 1.8E1 - (S * S * S) * (2.7E1 / 2.0),
                    S * 9.0 - R * S * (4.5E1 / 2.0) + R * (S * S) * 2.7E1 +
                      (R * R) * S * (2.7E1 / 2.0) - (S * S) * (4.5E1 / 2.0) +
                      (S * S * S) * (2.7E1 / 2.0),
                    R * S * 2.7E1 - R * (S * S) * 2.7E1 - (R * R) * S * 2.7E1};
  // duplicateScalarArray(10, 1, phi, dim, res);
  vectorLagrangeLayout(phi, 10, _nComponents, res);
}

template <int dim> std::vector<double> feSpaceVecTriP3<dim>::dLdr(double *r)
{
  double R = r[0];
  double S = r[1];
  double dldr[10] = {R * 1.8E1 + S * 1.8E1 - R * S * 2.7E1 - (R * R) * (2.7E1 / 2.0) -
                       (S * S) * (2.7E1 / 2.0) - 1.1E1 / 2.0,
                     R * -9.0 + (R * R) * (2.7E1 / 2.0) + 1.0,
                     0.,
                     R * -4.5E1 - S * (4.5E1 / 2.0) + R * S * 5.4E1 + (R * R) * (8.1E1 / 2.0) +
                       (S * S) * (2.7E1 / 2.0) + 9.0,
                     R * 3.6E1 + S * (9.0 / 2.0) - R * S * 2.7E1 - (R * R) * (8.1E1 / 2.0) -
                       9.0 / 2.0,
                     S * (-9.0 / 2.0) + R * S * 2.7E1,
                     S * (-9.0 / 2.0) + (S * S) * (2.7E1 / 2.0),
                     S * (9.0 / 2.0) - (S * S) * (2.7E1 / 2.0),
                     S * (-4.5E1 / 2.0) + R * S * 2.7E1 + (S * S) * 2.7E1,
                     S * 2.7E1 - R * S * 5.4E1 - (S * S) * 2.7E1};
  std::vector<double> res(_nFunctions*_nComponents);
  // duplicateScalarArray(10, 1, dldr, dim, res);
  vectorLagrangeLayout(dldr, 10, _nComponents, res.data());
  return res;
}

template <int dim> std::vector<double> feSpaceVecTriP3<dim>::dLds(double *r)
{
  double R = r[0];
  double S = r[1];
  double dlds[10] = {R * 1.8E1 + S * 1.8E1 - R * S * 2.7E1 - (R * R) * (2.7E1 / 2.0) -
                       (S * S) * (2.7E1 / 2.0) - 1.1E1 / 2.0,
                     0.,
                     S * -9.0 + (S * S) * (2.7E1 / 2.0) + 1.0,
                     R * (-4.5E1 / 2.0) + R * S * 2.7E1 + (R * R) * 2.7E1,
                     R * (9.0 / 2.0) - (R * R) * (2.7E1 / 2.0),
                     R * (-9.0 / 2.0) + (R * R) * (2.7E1 / 2.0),
                     R * (-9.0 / 2.0) + R * S * 2.7E1,
                     R * (9.0 / 2.0) + S * 3.6E1 - R * S * 2.7E1 - (S * S) * (8.1E1 / 2.0) -
                       9.0 / 2.0,
                     R * (-4.5E1 / 2.0) - S * 4.5E1 + R * S * 5.4E1 + (R * R) * (2.7E1 / 2.0) +
                       (S * S) * (8.1E1 / 2.0) + 9.0,
                     R * 2.7E1 - R * S * 5.4E1 - (R * R) * 2.7E1};
  std::vector<double> res(_nFunctions*_nComponents);
  // duplicateScalarArray(10, 1, dlds, dim, res);
  vectorLagrangeLayout(dlds, 10, _nComponents, res.data());
  return res;
}

template <int dim> void feSpaceVecTriP3<dim>::initializeNumberingUnknowns()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 0, dim);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 1, dim);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 2, dim);
    _numbering->setUnknownElementDOF(_mesh, _cncGeoID, i, dim);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, 2 * dim);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 1, 2 * dim);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 2, 2 * dim);
  }
}

template <int dim> void feSpaceVecTriP3<dim>::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 2);
    _numbering->setEssentialElementDOF(_mesh, _cncGeoID, i);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 2);
  }
}

template <int dim>
void feSpaceVecTriP3<dim>::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 1);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 0);
  adr[3] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 1);
  adr[4] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2, 0);
  adr[5] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2, 1);
  int e0 = _mesh->getEdge(_cncGeoID, numElem, 0);
  if(e0 > 0) {
    adr[6] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 0);
    adr[7] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 1);
    adr[8] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 2);
    adr[9] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 3);
  } else {
    adr[6] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 2);
    adr[7] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 3);
    adr[8] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 0);
    adr[9] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 1);
  }
  int e1 = _mesh->getEdge(_cncGeoID, numElem, 1);
  if(e1 > 0) {
    adr[10] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 0);
    adr[11] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 1);
    adr[12] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 2);
    adr[13] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 3);
  } else {
    adr[10] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 2);
    adr[11] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 3);
    adr[12] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 0);
    adr[13] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 1);
  }
  int e2 = _mesh->getEdge(_cncGeoID, numElem, 2);
  if(e2 > 0) {
    adr[14] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 0);
    adr[15] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 1);
    adr[16] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 2);
    adr[17] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 3);
  } else {
    adr[14] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 2);
    adr[15] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 3);
    adr[16] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 0);
    adr[17] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 1);
  }
  adr[18] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
  adr[19] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 1);
}

template class feSpaceVecTriP3<1>;
template class feSpaceVecTriP3<2>;

// -----------------------------------------------------------------------------
// Lagrange element of degree 4 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
/*
   2
   | \
   9   8
   |     \
  10  14   7
   |         \
  11  12  13   6
   |             \
   0---3---4---5---1

*/

feSpaceTriP4::feSpaceTriP4(std::string cncGeoID)
  : feScalarSpace(2, nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 15;
  double x0 = 0.0, x1 = 1. / 4., x2 = 1. / 2., x3 = 3. / 4., x4 = 1.;
  double y0 = 0.0, y1 = 1. / 4., y2 = 1. / 2., y3 = 3. / 4., y4 = 1.;
  _Lcoor = {x0, y0, 0., x4, y0, 0., x0, y4, 0., x1, y0, 0., x2, y0, 0.,
            x3, y0, 0., x3, y1, 0., x2, y2, 0., x1, y3, 0., y3, x0, 0.,
            y2, x0, 0., y1, x0, 0., x1, y1, 0., x2, y1, 0., x1, y2, 0.};
  _isGeometricInterpolant = true;
  _dofLocations = {dofLocation::VERTEX, dofLocation::VERTEX, dofLocation::VERTEX,
                   dofLocation::EDGE, dofLocation::EDGE, dofLocation::EDGE,
                   dofLocation::EDGE, dofLocation::EDGE, dofLocation::EDGE,
                   dofLocation::EDGE, dofLocation::EDGE, dofLocation::EDGE,
                   dofLocation::ELEMENT, dofLocation::ELEMENT, dofLocation::ELEMENT};
}

feSpaceTriP4::feSpaceTriP4(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                           feFunction *fct, const bool useGlobalShapeFunctions)
  : feScalarSpace(2, mesh, fieldID, cncGeoID, fct, useGlobalShapeFunctions)
{
  _nFunctions = 15;
  double x0 = 0.0, x1 = 1. / 4., x2 = 1. / 2., x3 = 3. / 4., x4 = 1.;
  double y0 = 0.0, y1 = 1. / 4., y2 = 1. / 2., y3 = 3. / 4., y4 = 1.;
  _Lcoor = {x0, y0, 0., x4, y0, 0., x0, y4, 0.,
            x1, y0, 0., x2, y0, 0., x3, y0, 0.,
            x3, y1, 0., x2, y2, 0., x1, y3, 0.,
            x0, y3, 0., x0, y2, 0., x0, y1, 0.,
            x1, y1, 0., x2, y1, 0., x1, y2, 0.};
  _dofLocations = {dofLocation::VERTEX, dofLocation::VERTEX, dofLocation::VERTEX,
                   dofLocation::EDGE, dofLocation::EDGE, dofLocation::EDGE,
                   dofLocation::EDGE, dofLocation::EDGE, dofLocation::EDGE,
                   dofLocation::EDGE, dofLocation::EDGE, dofLocation::EDGE,
                   dofLocation::ELEMENT, dofLocation::ELEMENT, dofLocation::ELEMENT};
}

std::vector<double> feSpaceTriP4::L(double *r)
{
  double R = r[0];
  double S = r[1];
  double f1 = 1. / 4. - R - S;
  double f2 = 1. / 2. - R - S;
  double f3 = 3. / 4. - R - S;
  double f4 = 1. - R - S;
  double r1 = R - 1. / 4.;
  double r2 = R - 1. / 2.;
  double r3 = R - 3. / 4.;
  double s1 = S - 1. / 4.;
  double s2 = S - 1. / 2.;
  double s3 = S - 3. / 4.;
  return {f1 * f2 * f3 * f4 * 32. / 3., R * r1 * r2 * r3 * 32. / 3., S * s1 * s2 * s3 * 32. / 3.,
          R * f2 * f3 * f4 * 128. / 3., R * r1 * f3 * f4 * 64.,      R * r1 * r2 * f4 * 128. / 3.,
          S * R * r1 * r2 * 128. / 3.,  R * S * r1 * s1 * 64.,       R * S * s1 * s2 * 128. / 3.,
          S * s1 * s2 * f4 * 128. / 3., S * s1 * f3 * f4 * 64.,      S * f2 * f3 * f4 * 128. / 3.,
          R * S * f3 * f4 * 128.,       R * S * r1 * f4 * 128.,      R * S * s1 * f4 * 128.};
}

void feSpaceTriP4::L(double *r, double *L)
{
  double R = r[0];
  double S = r[1];
  double f1 = 1. / 4. - R - S;
  double f2 = 1. / 2. - R - S;
  double f3 = 3. / 4. - R - S;
  double f4 = 1. - R - S;
  double r1 = R - 1. / 4.;
  double r2 = R - 1. / 2.;
  double r3 = R - 3. / 4.;
  double s1 = S - 1. / 4.;
  double s2 = S - 1. / 2.;
  double s3 = S - 3. / 4.;
  L[0] = f1 * f2 * f3 * f4 * 32. / 3.;
  L[1] =  R * r1 * r2 * r3 * 32. / 3.;
  L[2] =  S * s1 * s2 * s3 * 32. / 3.;
  L[3] = R * f2 * f3 * f4 * 128. / 3.;
  L[4] = R * r1 * f3 * f4 * 64.;
  L[5] = R * r1 * r2 * f4 * 128. / 3.,
  L[6] = S * R * r1 * r2 * 128. / 3.;
  L[7] = R * S * r1 * s1 * 64.;
  L[8] = R * S * s1 * s2 * 128. / 3.,
  L[9] = S * s1 * s2 * f4 * 128. / 3.;
  L[10] = S * s1 * f3 * f4 * 64.;
  L[11] = S * f2 * f3 * f4 * 128. / 3.,
  L[12] = R * S * f3 * f4 * 128.;
  L[13] = R * S * r1 * f4 * 128.;
  L[14] = R * S * s1 * f4 * 128.;
}

std::vector<double> feSpaceTriP4::dLdr(double *r)
{
  double R = r[0];
  double S = r[1];
  return {32. / 3. *
            (4. * R * R * R + 12. * R * R * S - (15. * R * R) / 2. + 12. * R * S * S - 15. * R * S +
             (35. * R) / 8. + 4. * S * S * S - (15. * S * S) / 2. + (35. * S) / 8. - 25. / 32.),
          32. / 3. * (4. * R * R * R - (9. * R * R) / 2. + (11. * R) / 8. - 3. / 32.),
          0.,
          128. / 3. *
            (-4. * R * R * R - 9. * R * R * S + (27. * R * R) / 4. - 6. * R * S * S + 9. * R * S -
             (13. * R) / 4. - S * S * S + (9. * S * S) / 4. - (13. * S) / 8. + 3. / 8.),
          64. * (4. * R * R * R + 6. * R * R * S - 6. * R * R + 2. * R * S * S - (9. * R * S) / 2. +
                 (19. * R) / 8. - S * S / 4. + (7. * S) / 16. - 3. / 16.),
          128. / 3. *
            ((3. * R * S) / 2. - S / 8. - (7. * R) / 4. - 3. * R * R * S + (21. * R * R) / 4. -
             4. * R * R * R + 1. / 8.),
          128. / 3. * ((S * (24. * R * R - 12. * R + 1.)) / 8.),
          64. * ((S * (8. * R - 1.) * (4. * S - 1.)) / 16.),
          128. / 3. * (S * (S - 1. / 2.) * (S - 1. / 4.)),
          128. / 3. * (-S * (S - 1. / 2.) * (S - 1. / 4.)),
          64. * ((S * (4. * S - 1.) * (8. * R + 8. * S - 7.)) / 16.),
          128. / 3. *
            (-(S * (24. * R * R + 48. * R * S - 36. * R + 24. * S * S - 36. * S + 13.)) / 8.),
          128. * ((S * (12. * R * R + 16. * R * S - 14. * R + 4. * S * S - 7. * S + 3.)) / 4.),
          128. * (-(S * (8. * R * S - S - 10. * R + 12. * R * R + 1.)) / 4.),
          128. * (-(S * (4. * S - 1.) * (2. * R + S - 1.)) / 4.)};
}

std::vector<double> feSpaceTriP4::dLds(double *r)
{
  double R = r[0];
  double S = r[1];
  return {32. / 3. *
            (4. * R * R * R + 12. * R * R * S - (15. * R * R) / 2. + 12. * R * S * S - 15. * R * S +
             (35. * R) / 8. + 4. * S * S * S - (15. * S * S) / 2. + (35. * S) / 8. - 25. / 32.),
          0.,
          32. / 3. * (4. * S * S * S - (9. * S * S) / 2. + (11. * S) / 8. - 3. / 32.),
          128. / 3. *
            (-(R * (24. * R * R + 48. * R * S - 36. * R + 24. * S * S - 36. * S + 13.)) / 8.),
          64. * ((R * (4. * R - 1.) * (8. * R + 8. * S - 7.)) / 16.),
          128. / 3. * (-R * (R - 1. / 2.) * (R - 1. / 4.)),
          128. / 3. * (R * (R - 1. / 2.) * (R - 1. / 4.)),
          64. * ((R * (4. * R - 1.) * (8. * S - 1.)) / 16.),
          128. / 3. *
            (R * (S - 1. / 2.) * (S - 1. / 4.) + R * S * (S - 1. / 2.) +
             R * S * (S - 1. / 4.)), //(R*(48.*S*S - 16.*S + 1.))/16.,
          128. / 3. *
            ((3. * R * S) / 2. - (7. * S) / 4. - R / 8. - 3. * R * S * S + (21. * S * S) / 4. -
             4. * S * S * S + 1. / 8.),
          64. * (2. * R * R * S - R * R / 4. + 6. * R * S * S - (9. * R * S) / 2. + (7. * R) / 16. +
                 4. * S * S * S - 6. * S * S + (19. * S) / 8. - 3. / 16.),
          128. / 3. *
            (-R * R * R - 6. * R * R * S + (9. * R * R) / 4. - 9. * R * S * S + 9. * R * S -
             (13. * R) / 8. - 4. * S * S * S + (27. * S * S) / 4. - (13. * S) / 4. + 3. / 8.),
          128. * ((R * (4. * R * R + 16. * R * S - 7. * R + 12. * S * S - 14. * S + 3.)) / 4.),
          128. * (-(R * (4. * R - 1.) * (R + 2. * S - 1.)) / 4.),
          128. * (-(R * (8. * R * S - 10. * S - R + 12. * S * S + 1.)) / 4.)};
}

void feSpaceTriP4::initializeNumberingUnknowns()
{
  int nDOFPerEdge = 3;
  int nDOFPerElem = 3;
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 2);
    _numbering->setUnknownElementDOF(_mesh, _cncGeoID, i, nDOFPerElem);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, nDOFPerEdge);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 1, nDOFPerEdge);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 2, nDOFPerEdge);
  }
}

void feSpaceTriP4::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 2);
    _numbering->setEssentialElementDOF(_mesh, _cncGeoID, i);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 2);
  }
}

void feSpaceTriP4::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2);
  int e0 = _mesh->getEdge(_cncGeoID, numElem, 0);
  if(e0 > 0) {
    adr[3] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 0);
    adr[4] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 1);
    adr[5] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 2);
  } else {
    adr[3] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 2);
    adr[4] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 1);
    adr[5] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 0);
  }
  int e1 = _mesh->getEdge(_cncGeoID, numElem, 1);
  if(e1 > 0) {
    adr[6] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 0);
    adr[7] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 1);
    adr[8] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 2);
  } else {
    adr[6] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 2);
    adr[7] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 1);
    adr[8] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 0);
  }
  int e2 = _mesh->getEdge(_cncGeoID, numElem, 2);
  if(e2 > 0) {
    adr[9] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 0);
    adr[10] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 1);
    adr[11] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 2);
  } else {
    adr[9] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 2);
    adr[10] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 1);
    adr[11] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 0);
  }
  adr[12] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
  adr[13] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 1);
  adr[14] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 2);
}

// -----------------------------------------------------------------------------
// Vector Lagrange element of degree 4 on reference triangle r = [0,1], s = [0,1-r]
// -----------------------------------------------------------------------------
template <int dim>
feSpaceVecTriP4<dim>::feSpaceVecTriP4(feMesh *mesh, const std::string fieldID,
                                      const std::string cncGeoID, feVectorFunction *fct,
                                      const bool useGlobalShapeFunctions)
  : feVectorSpace(2, mesh, fieldID, cncGeoID, fct, useGlobalShapeFunctions)
{
  _nComponents = dim;
  _nFunctions = 15 * dim;
  double x0 = 0.0, x1 = 1. / 4., x2 = 1. / 2., x3 = 3. / 4., x4 = 1.;
  double y0 = 0.0, y1 = 1. / 4., y2 = 1. / 2., y3 = 3. / 4., y4 = 1.;
  double coor[15 * 3] = {x0, y0, 0., x4, y0, 0., x0, y4, 0.,
                         x1, y0, 0., x2, y0, 0., x3, y0, 0.,
                         x3, y1, 0., x2, y2, 0., x1, y3, 0.,
                         x0, y3, 0., x0, y2, 0., x0, y1, 0.,
                         x1, y1, 0., x2, y1, 0., x1, y2, 0.};
  _Lcoor.resize(15 * 3 * dim);
  duplicateScalarArray(15 * 3, 3, coor, dim, _Lcoor);
  _dofLocations.resize(_nFunctions);
  for(int i = 0; i < _nFunctions; ++i) {
    if(i < 3 * dim)
      _dofLocations[i] = dofLocation::VERTEX;
    else if(i < 12 * dim)
      _dofLocations[i] = dofLocation::EDGE;
    else
      _dofLocations[i] = dofLocation::ELEMENT;
  }
}

template <int dim> std::vector<double> feSpaceVecTriP4<dim>::L(double *r)
{
  double R = r[0];
  double S = r[1];
  double f1 = 1. / 4. - R - S;
  double f2 = 1. / 2. - R - S;
  double f3 = 3. / 4. - R - S;
  double f4 = 1. - R - S;
  double r1 = R - 1. / 4.;
  double r2 = R - 1. / 2.;
  double r3 = R - 3. / 4.;
  double s1 = S - 1. / 4.;
  double s2 = S - 1. / 2.;
  double s3 = S - 3. / 4.;
  double phi[15] =  {f1 * f2 * f3 * f4 * 32. / 3., R * r1 * r2 * r3 * 32. / 3., S * s1 * s2 * s3 * 32. / 3.,
                      R * f2 * f3 * f4 * 128. / 3., R * r1 * f3 * f4 * 64.,      R * r1 * r2 * f4 * 128. / 3.,
                      S * R * r1 * r2 * 128. / 3.,  R * S * r1 * s1 * 64.,       R * S * s1 * s2 * 128. / 3.,
                      S * s1 * s2 * f4 * 128. / 3., S * s1 * f3 * f4 * 64.,      S * f2 * f3 * f4 * 128. / 3.,
                      R * S * f3 * f4 * 128.,       R * S * r1 * f4 * 128.,      R * S * s1 * f4 * 128.};
  std::vector<double> res(_nFunctions*_nComponents);
  // duplicateScalarArray(15, 1, phi, dim, res);
  vectorLagrangeLayout(phi, 15, _nComponents, res.data());
  return res;
}

template <int dim> void feSpaceVecTriP4<dim>::L(double *r, double *res)
{
  double R = r[0];
  double S = r[1];
  double phi[15] = {R * (-1.1E1 / 2.0) - S * (1.1E1 / 2.0) + R * S * 1.8E1 -
                      R * (S * S) * (2.7E1 / 2.0) - (R * R) * S * (2.7E1 / 2.0) + (R * R) * 9.0 -
                      (R * R * R) * (9.0 / 2.0) + (S * S) * 9.0 - (S * S * S) * (9.0 / 2.0) + 1.0,
                    R - (R * R) * (9.0 / 2.0) + (R * R * R) * (9.0 / 2.0),
                    S - (S * S) * (9.0 / 2.0) + (S * S * S) * (9.0 / 2.0),
                    R * 9.0 - R * S * (4.5E1 / 2.0) + R * (S * S) * (2.7E1 / 2.0) +
                      (R * R) * S * 2.7E1 - (R * R) * (4.5E1 / 2.0) + (R * R * R) * (2.7E1 / 2.0),
                    R * (-9.0 / 2.0) + R * S * (9.0 / 2.0) - (R * R) * S * (2.7E1 / 2.0) +
                      (R * R) * 1.8E1 - (R * R * R) * (2.7E1 / 2.0),
                    R * S * (-9.0 / 2.0) + (R * R) * S * (2.7E1 / 2.0),
                    R * S * (-9.0 / 2.0) + R * (S * S) * (2.7E1 / 2.0),
                    S * (-9.0 / 2.0) + R * S * (9.0 / 2.0) - R * (S * S) * (2.7E1 / 2.0) +
                      (S * S) * 1.8E1 - (S * S * S) * (2.7E1 / 2.0),
                    S * 9.0 - R * S * (4.5E1 / 2.0) + R * (S * S) * 2.7E1 +
                      (R * R) * S * (2.7E1 / 2.0) - (S * S) * (4.5E1 / 2.0) +
                      (S * S * S) * (2.7E1 / 2.0),
                    R * S * 2.7E1 - R * (S * S) * 2.7E1 - (R * R) * S * 2.7E1};
  // duplicateScalarArray(15, 1, phi, dim, res);
  vectorLagrangeLayout(phi, 15, _nComponents, res);
}

template <int dim> std::vector<double> feSpaceVecTriP4<dim>::dLdr(double *r)
{
  double R = r[0];
  double S = r[1];
  double dldr[15] = {32. / 3. *
                      (4. * R * R * R + 12. * R * R * S - (15. * R * R) / 2. + 12. * R * S * S - 15. * R * S +
                       (35. * R) / 8. + 4. * S * S * S - (15. * S * S) / 2. + (35. * S) / 8. - 25. / 32.),
                    32. / 3. * (4. * R * R * R - (9. * R * R) / 2. + (11. * R) / 8. - 3. / 32.),
                    0.,
                    128. / 3. *
                      (-4. * R * R * R - 9. * R * R * S + (27. * R * R) / 4. - 6. * R * S * S + 9. * R * S -
                       (13. * R) / 4. - S * S * S + (9. * S * S) / 4. - (13. * S) / 8. + 3. / 8.),
                    64. * (4. * R * R * R + 6. * R * R * S - 6. * R * R + 2. * R * S * S - (9. * R * S) / 2. +
                           (19. * R) / 8. - S * S / 4. + (7. * S) / 16. - 3. / 16.),
                    128. / 3. *
                      ((3. * R * S) / 2. - S / 8. - (7. * R) / 4. - 3. * R * R * S + (21. * R * R) / 4. -
                       4. * R * R * R + 1. / 8.),
                    128. / 3. * ((S * (24. * R * R - 12. * R + 1.)) / 8.),
                    64. * ((S * (8. * R - 1.) * (4. * S - 1.)) / 16.),
                    128. / 3. * (S * (S - 1. / 2.) * (S - 1. / 4.)),
                    128. / 3. * (-S * (S - 1. / 2.) * (S - 1. / 4.)),
                    64. * ((S * (4. * S - 1.) * (8. * R + 8. * S - 7.)) / 16.),
                    128. / 3. *
                      (-(S * (24. * R * R + 48. * R * S - 36. * R + 24. * S * S - 36. * S + 13.)) / 8.),
                    128. * ((S * (12. * R * R + 16. * R * S - 14. * R + 4. * S * S - 7. * S + 3.)) / 4.),
                    128. * (-(S * (8. * R * S - S - 10. * R + 12. * R * R + 1.)) / 4.),
                    128. * (-(S * (4. * S - 1.) * (2. * R + S - 1.)) / 4.)};
  std::vector<double> res(_nFunctions*_nComponents);
  // duplicateScalarArray(15, 1, dldr, dim, res);
  vectorLagrangeLayout(dldr, 15, _nComponents, res.data());
  return res;
}

template <int dim> std::vector<double> feSpaceVecTriP4<dim>::dLds(double *r)
{
  double R = r[0];
  double S = r[1];
  double dlds[15] = {32. / 3. *
                      (4. * R * R * R + 12. * R * R * S - (15. * R * R) / 2. + 12. * R * S * S - 15. * R * S +
                       (35. * R) / 8. + 4. * S * S * S - (15. * S * S) / 2. + (35. * S) / 8. - 25. / 32.),
                    32. / 3. * (0.0),
                    32. / 3. * (4. * S * S * S - (9. * S * S) / 2. + (11. * S) / 8. - 3. / 32.),
                    128. / 3. *
                      (-(R * (24. * R * R + 48. * R * S - 36. * R + 24. * S * S - 36. * S + 13.)) / 8.),
                    64. * ((R * (4. * R - 1.) * (8. * R + 8. * S - 7.)) / 16.),
                    128. / 3. * (-R * (R - 1. / 2.) * (R - 1. / 4.)),
                    128. / 3. * (R * (R - 1. / 2.) * (R - 1. / 4.)),
                    64. * ((R * (4. * R - 1.) * (8. * S - 1.)) / 16.),
                    128. / 3. *
                      (R * (S - 1. / 2.) * (S - 1. / 4.) + R * S * (S - 1. / 2.) +
                       R * S * (S - 1. / 4.)), //(R*(48.*S*S - 16.*S + 1.))/16.,
                    128. / 3. *
                      ((3. * R * S) / 2. - (7. * S) / 4. - R / 8. - 3. * R * S * S + (21. * S * S) / 4. -
                       4. * S * S * S + 1. / 8.),
                    64. * (2. * R * R * S - R * R / 4. + 6. * R * S * S - (9. * R * S) / 2. + (7. * R) / 16. +
                           4. * S * S * S - 6. * S * S + (19. * S) / 8. - 3. / 16.),
                    128. / 3. *
                      (-R * R * R - 6. * R * R * S + (9. * R * R) / 4. - 9. * R * S * S + 9. * R * S -
                       (13. * R) / 8. - 4. * S * S * S + (27. * S * S) / 4. - (13. * S) / 4. + 3. / 8.),
                    128. * ((R * (4. * R * R + 16. * R * S - 7. * R + 12. * S * S - 14. * S + 3.)) / 4.),
                    128. * (-(R * (4. * R - 1.) * (R + 2. * S - 1.)) / 4.),
                    128. * (-(R * (8. * R * S - 10. * S - R + 12. * S * S + 1.)) / 4.)};
  std::vector<double> res(_nFunctions*_nComponents);
  // duplicateScalarArray(15, 1, dlds, dim, res);
  vectorLagrangeLayout(dlds, 15, _nComponents, res.data());
  return res;
}

template <int dim> void feSpaceVecTriP4<dim>::initializeNumberingUnknowns()
{
  int nDOFPerEdge = 3;
  int nDOFPerElem = 3;
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 0, dim);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 1, dim);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 2, dim);
    _numbering->setUnknownElementDOF(_mesh, _cncGeoID, i, nDOFPerElem * dim);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, nDOFPerEdge * dim);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 1, nDOFPerEdge * dim);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 2, nDOFPerEdge * dim);
  }
}

template <int dim> void feSpaceVecTriP4<dim>::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 2);
    _numbering->setEssentialElementDOF(_mesh, _cncGeoID, i);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 2);
  }
}

template <int dim>
void feSpaceVecTriP4<dim>::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 1);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 0);
  adr[3] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 1);
  adr[4] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2, 0);
  adr[5] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2, 1);
  int e0 = _mesh->getEdge(_cncGeoID, numElem, 0);
  if(e0 > 0) {
    adr[6] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 0);
    adr[7] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 1);
    adr[8] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 2);
    adr[9] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 3);
    adr[10] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 4);
    adr[11] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 5);
  } else {
    adr[6] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 4);
    adr[7] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 5);
    adr[8] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 2);
    adr[9] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 3);
    adr[10] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 0);
    adr[11] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 0, 1);
  }
  int e1 = _mesh->getEdge(_cncGeoID, numElem, 1);
  if(e1 > 0) {
    adr[12] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 0);
    adr[13] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 1);
    adr[14] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 2);
    adr[15] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 3);
    adr[16] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 4);
    adr[17] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 5);
  } else {
    adr[12] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 4);
    adr[13] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 5);
    adr[14] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 2);
    adr[15] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 3);
    adr[16] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 0);
    adr[17] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 1, 1);
  }
  int e2 = _mesh->getEdge(_cncGeoID, numElem, 2);
  if(e2 > 0) {
    adr[18] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 0);
    adr[19] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 1);
    adr[20] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 2);
    adr[21] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 3);
    adr[22] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 4);
    adr[23] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 5);
  } else {
    adr[18] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 4);
    adr[19] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 5);
    adr[20] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 2);
    adr[21] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 3);
    adr[22] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 0);
    adr[23] = _numbering->getEdgeDOF(_mesh, _cncGeoID, numElem, 2, 1);
  }
  adr[24] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
  adr[25] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 1);
  adr[26] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 2);
  adr[27] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 3);
  adr[28] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 4);
  adr[29] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 5);
}

template class feSpaceVecTriP4<1>;
template class feSpaceVecTriP4<2>;