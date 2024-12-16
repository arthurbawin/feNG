#include "feSpace.h"
#include "feMesh.h"
#include "feNumber.h"
#include "fePolynomial.h"
#include "feTetrahedron.h"

// // Copy array "dim" times by chunks of size "chunksize"
// static void duplicateScalarArray(const int scalarSize, const int chunkSize,
//                                  const double *scalarArray, const int dim,
//                                  std::vector<double> &vectorArray)
// {
//   int nChunks = scalarSize / chunkSize;
//   int cnt = 0;
//   for(int j = 0; j < nChunks; ++j) {
//     for(int iDim = 0; iDim < dim; ++iDim) {
//       for(int i = 0; i < chunkSize; ++i) {
//         vectorArray[cnt++] = scalarArray[chunkSize * j + i];
//       }
//     }
//   }
// }

// static void duplicateScalarArray(const int scalarSize, const int chunkSize,
//                                  const double *scalarArray, const int dim, double *vectorArray)
// {
//   int nChunks = scalarSize / chunkSize;
//   int cnt = 0;
//   for(int j = 0; j < nChunks; ++j) {
//     for(int iDim = 0; iDim < dim; ++iDim) {
//       for(int i = 0; i < chunkSize; ++i) {
//         feInfo("Writing entry %d of scalararray of size %d at place %d", chunkSize * j + i,
//                scalarSize, cnt + 1);
//         vectorArray[cnt++] = scalarArray[chunkSize * j + i];
//       }
//     }
//   }
// }

// -----------------------------------------------------------------------------
// Lagrange element of degree n on reference tetrahedron
// -----------------------------------------------------------------------------
feSpaceTetPn::feSpaceTetPn(int n, std::string cncGeoID)
  // : feScalarSpace(3, nullptr, "GEO", cncGeoID, nullptr), _n(n)
 : feSpaceTetPn(n, nullptr, "GEO", cncGeoID, nullptr, false)
{
  _isGeometricInterpolant = true;
}

feSpaceTetPn::feSpaceTetPn(int n, feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                           feFunction *fct, const bool useGlobalShapeFunctions)
  : feScalarSpace(3, mesh, fieldID, cncGeoID, fct, useGlobalShapeFunctions), _n(n)
{
  _nFunctions = (n+1)*(n+2)*(n+3)/6;
  _Lcoor.clear();
  _Lcoor.reserve(3*(n+1)*(n+2)*(n+3)/6);

  // Get reference barycentric coordinates, see feTetrahedron.cpp
  getTetLagrangeBarycentricCoord(_n, _refBarycentric);
  assert(_refBarycentric.size()/4 == (size_t) _nFunctions);

  // The Lagrange points are the barycentric components [2, ..., n+1] divided by n:
  for(int i = 0; i < _nFunctions; ++i) {
#if defined(FENG_DEBUG)
    double sum = _refBarycentric[4*i+0] + _refBarycentric[4*i+1]
               + _refBarycentric[4*i+2] + _refBarycentric[4*i+3];
    assert(fabs(sum - (double) _n) < 1e-14);
#endif
    
    // Add components [2, ..., n+1]
    _Lcoor.push_back(_refBarycentric[4*i+1] / (double) n);
    _Lcoor.push_back(_refBarycentric[4*i+2] / (double) n);
    _Lcoor.push_back(_refBarycentric[4*i+3] / (double) n);
  }

  // for(int i = 0; i < _nFunctions; ++i) {
  //   feInfo("%1.3f - %1.3f - %1.3f",
  //     _Lcoor[3*i+0],_Lcoor[3*i+1],_Lcoor[3*i+2]);
  // }

  _dofLocations.resize(_nFunctions);
  int start = 0;
  for(int i = 0; i < 4; ++i)
    { _dofLocations[i] = dofLocation::VERTEX; }
  if(n >= 2) {
    start += 4;
    for(int i = start; i < start + 6*(n-1); ++i)
      { _dofLocations[i] = dofLocation::EDGE; }
  }
  if(n >= 3) {
    start += 6*(n-1);
    for(int i = start; i < start + 4*(n-1)*(n-2)/2; ++i)
      { _dofLocations[i] = dofLocation::FACE; }
  }
  if(n >= 4) {
    start += 4*(n-1)*(n-2)/2;
    for(int i = start; i < _nFunctions; ++i)
      { _dofLocations[i] = dofLocation::ELEMENT; }
  }
}

static int factorial(const int n)
{
  int f = 1;
  for(int i = 1; i <= n; ++i)
    f *= i;
  return f;
}

std::vector<double> feSpaceTetPn::L(double *r)
{ 
  std::vector<double> phi(_nFunctions);
  this->L(r,phi.data());
  return phi;
}

void feSpaceTetPn::L(double *r, double *phi)
{
  double u[4] = {1.0 - r[0] - r[1] - r[2], r[0], r[1], r[2]};

  for(int iF = 0; iF < _nFunctions; ++iF) {
    int alpha[4] = {(int) _refBarycentric[4*iF+0], (int) _refBarycentric[4*iF+1], (int) _refBarycentric[4*iF+2], (int) _refBarycentric[4*iF+3]};
    double res = 1., factor = 1.;
    for(int i = 0; i < 4; ++i) {
      factor *= (double) factorial(alpha[i]);
      for(int j = 0; j <= alpha[i]-1; ++j) {
        res *= u[i] - (double) j / (double) _n;
      }
    }
    res *= pow(_n,_n) / factor;
    phi[iF] = res;
  }
}

std::vector<double> feSpaceTetPn::dLdr(double *r) { UNUSED(r); return {-1., 1., 0., 0.}; }
std::vector<double> feSpaceTetPn::dLds(double *r) { UNUSED(r); return {-1., 0., 1., 0.}; }
std::vector<double> feSpaceTetPn::dLdt(double *r) { UNUSED(r); return {-1., 0., 0., 1.}; }
std::vector<double> feSpaceTetPn::d2Ldr2(double *r) { UNUSED(r); return {0., 0., 0., 0.}; }
std::vector<double> feSpaceTetPn::d2Ldrs(double *r) { UNUSED(r); return {0., 0., 0., 0.}; }
std::vector<double> feSpaceTetPn::d2Lds2(double *r) { UNUSED(r); return {0., 0., 0., 0.}; }

void feSpaceTetPn::initializeNumberingUnknowns()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 2);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 3);
  }
}

void feSpaceTetPn::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 2);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 3);
  }
}

void feSpaceTetPn::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 3);
}

// -----------------------------------------------------------------------------
// Lagrange element of degree 1 on reference tetrahedron
// -----------------------------------------------------------------------------
feSpaceTetP1::feSpaceTetP1(std::string cncGeoID)
  : feScalarSpace(3, nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 4;
  _Lcoor = {0., 0., 0., 1., 0., 0., 0., 1., 0., 0., 0., 1.};
  _isGeometricInterpolant = true;
  _dofLocations = {dofLocation::VERTEX, dofLocation::VERTEX,
                   dofLocation::VERTEX, dofLocation::VERTEX};
}

feSpaceTetP1::feSpaceTetP1(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                           feFunction *fct, const bool useGlobalShapeFunctions)
  : feScalarSpace(3, mesh, fieldID, cncGeoID, fct, useGlobalShapeFunctions)
{
  _nFunctions = 4;
  _Lcoor = {0., 0., 0., 1., 0., 0., 0., 1., 0., 0., 0., 1.};
  _dofLocations = {dofLocation::VERTEX, dofLocation::VERTEX,
                   dofLocation::VERTEX, dofLocation::VERTEX};
}

std::vector<double> feSpaceTetP1::L(double *r) { return {1.0 - r[0] - r[1] - r[2], r[0], r[1], r[2]}; }
void feSpaceTetP1::L(double *r, double *L)
{
  L[0] = 1. - r[0] - r[1] - r[2];
  L[1] = r[0];
  L[2] = r[1];
  L[3] = r[2];
}
std::vector<double> feSpaceTetP1::dLdr(double *r) { UNUSED(r); return {-1., 1., 0., 0.}; }
std::vector<double> feSpaceTetP1::dLds(double *r) { UNUSED(r); return {-1., 0., 1., 0.}; }
std::vector<double> feSpaceTetP1::dLdt(double *r) { UNUSED(r); return {-1., 0., 0., 1.}; }
std::vector<double> feSpaceTetP1::d2Ldr2(double *r) { UNUSED(r); return {0., 0., 0., 0.}; }
std::vector<double> feSpaceTetP1::d2Ldrs(double *r) { UNUSED(r); return {0., 0., 0., 0.}; }
std::vector<double> feSpaceTetP1::d2Lds2(double *r) { UNUSED(r); return {0., 0., 0., 0.}; }

void feSpaceTetP1::initializeNumberingUnknowns()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 2);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 3);
  }
}

void feSpaceTetP1::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 2);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 3);
  }
}

void feSpaceTetP1::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 3);
}