#include "feSpace.h"
#include "feMesh.h"
#include "feNumber.h"
#include "fePolynomial.h"

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

