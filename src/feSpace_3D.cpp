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
                           const feFunction *fct, const bool useGlobalShapeFunctions)
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

std::vector<double> feSpaceTetPn::L(const double *r) const
{ 
  std::vector<double> phi(_nFunctions);
  this->L(r, phi.data());
  return phi;
}

void feSpaceTetPn::L(const double *r, double *phi) const
{
  double u[4] = {1.0 - r[0] - r[1] - r[2], r[0], r[1], r[2]};

#if defined(FENG_DEBUG)
  double sum = 0.;
#endif

  for(int iF = 0; iF < _nFunctions; ++iF)
  {
    int alpha[4] = {(int) _refBarycentric[4*iF+0],
                    (int) _refBarycentric[4*iF+1],
                    (int) _refBarycentric[4*iF+2],
                    (int) _refBarycentric[4*iF+3]};

    double res = 1., factor = 1.;
    for(int i = 0; i < 4; ++i) {
      factor *= (double) factorial(alpha[i]);
      for(int j = 0; j <= alpha[i]-1; ++j) {
        res *= u[i] - (double) j / (double) _n;
      }
    }
    res *= pow((double) _n,_n) / factor;
    phi[iF] = res;
#if defined(FENG_DEBUG)
    sum += res;
#endif
  }

#if defined(FENG_DEBUG)
  assert(fabs(sum - 1.) < 1e-14);
#endif
}

// Jacobian matrix dg_ij/dxsi^m where g_ij(xsi) = (u^i(xsi) - j/k)
// and u(xsi) are the barycentric coordinates written in terms of the
// reference coordinates xsi.
static double jacg[4][3] = {{-1., -1., -1.}, {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

//
// Compute the derivative of all Lagrange shape functions w.r.t. the m-th component of xsi.
// Obtained by applying the product rule to the general expression of the shape functions:
//
//                        n^n         dim+1 alpha_i-1
// phi_alpha(xsi) = --------------- * prod  prod      (u_i(xsi) - j/k)
//                  dim +1            i = 1 j = 0
//                   prod  alpha_i!
//                  i = 1
//
// For the product rule of an arbitrary number of functions,
// see e.g. https://en.wikipedia.org/wiki/Product_rule#Generalizations
//
std::vector<double> feSpaceTetPn::shapeTetDerivatives(const double xsi[3], const int m) const
{
  std::vector<double> res(_nFunctions);

  // Barycentric coordinates
  double u[4] = {1. - xsi[0] - xsi[1] - xsi[2], xsi[0], xsi[1], xsi[2]};

  double nn = pow((double) _n, _n);

#if defined(FENG_DEBUG)
  double sum = 0.;
#endif

  for(int iPhi = 0; iPhi < _nFunctions; ++iPhi)
  {
    int alpha[4] = {(int) _refBarycentric[4*iPhi+0],
                    (int) _refBarycentric[4*iPhi+1],
                    (int) _refBarycentric[4*iPhi+2],
                    (int) _refBarycentric[4*iPhi+3]};

    double factor = 1.;
    for(int i = 0; i < 4; ++i) {
      factor *= (double) factorial(alpha[i]);
    }
    double C_alpha = nn / factor;

    double h[4] = {1., 1., 1., 1.};

    for(int i = 0; i < 4; ++i) {
      for(int j = 0; j <= alpha[i]-1; ++j) {
        h[i] *= u[i] - (double) j / (double) _n;
      }
    }

    double dhdxsi_m[4] = {0., 0., 0., 0.};

    for(int i = 0; i < 4; ++i) {
      // j is constant in gij, so derivatives of all gij with variable j are the same
      double dgijdxsi_m = jacg[i][m];
      for(int j = 0; j <= alpha[i]-1; ++j) {
        double prod = 1.;
        for(int s = 0; s <= alpha[i]-1; ++s) {
          if(s != j) {
            prod *= u[i] - (double) s / (double) _n;
          }
        }
        dhdxsi_m[i] += dgijdxsi_m * prod;
      }
    }

    double dphidxsi_m = 0.;
    for(int i = 0; i < 4; ++i) {
      double prod = 1.;
      for(int s = 0; s < 4; ++s) {
        if(s != i) {
          prod *= h[s];
        }
      }
      dphidxsi_m += dhdxsi_m[i] * prod;
    }

    res[iPhi] = C_alpha * dphidxsi_m;
#if defined(FENG_DEBUG)
    sum += res[iPhi];
#endif
  }

#if defined(FENG_DEBUG)
  assert(fabs(sum) < 1e-14);
#endif

  return res;
}

std::vector<double> feSpaceTetPn::dLdr(const double *r) const { return this->shapeTetDerivatives(r, 0); }
std::vector<double> feSpaceTetPn::dLds(const double *r) const { return this->shapeTetDerivatives(r, 1); }
std::vector<double> feSpaceTetPn::dLdt(const double *r) const { return this->shapeTetDerivatives(r, 2); }
std::vector<double> feSpaceTetPn::d2Ldr2(const double */*r*/) const { std::vector<double> res(_nFunctions, 0.); return res; }
std::vector<double> feSpaceTetPn::d2Ldrs(const double */*r*/) const { std::vector<double> res(_nFunctions, 0.); return res; }
std::vector<double> feSpaceTetPn::d2Lds2(const double */*r*/) const { std::vector<double> res(_nFunctions, 0.); return res; }

void feSpaceTetPn::initializeNumberingUnknowns()
{
  int nDOFPerVert = 1;
  int nDOFPerEdge = _n-1;
  int nDOFPerFace = (_n-1)*(_n-2)/2;
  int nDOFPerElem = _nFunctions - 4 - 6 * nDOFPerEdge - 4 * nDOFPerFace;

  for(int i = 0; i < _mesh->getNumElements(_cncGeoTag); ++i) {
    for(int j = 0; j < 4; ++j) {
      _numbering->setUnknownVertexDOF(_mesh, _cncGeoTag, i, j, nDOFPerVert);
    }
    if(_n >= 2) {
      for(int j = 0; j < 6; ++j) {
        _numbering->setUnknownEdgeDOF(_mesh, _cncGeoTag, i, j, nDOFPerEdge);
      }
    }
    if(_n >= 3) {
      for(int j = 0; j < 4; ++j) {
        _numbering->setUnknownFaceDOF(_mesh, _cncGeoTag, i, j, nDOFPerFace);
      }
    }
    if(_n >= 4) {
      _numbering->setUnknownElementDOF(_mesh, _cncGeoTag, i, nDOFPerElem);
    }
  }
}

void feSpaceTetPn::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoTag); ++i) {
    for(int j = 0; j < 4; ++j) {
      _numbering->setEssentialVertexDOF(_mesh, _cncGeoTag, i, j);
    }
    if(_n >= 2) {
      for(int j = 0; j < 6; ++j) {
        _numbering->setEssentialEdgeDOF(_mesh, _cncGeoTag, i, j);
      }
    }
    if(_n >= 3) {
      for(int j = 0; j < 4; ++j) {
        _numbering->setEssentialFaceDOF(_mesh, _cncGeoTag, i, j);
      }
    }
    if(_n >= 4) {
      _numbering->setEssentialElementDOF(_mesh, _cncGeoTag, i);
    }
  }
}

void feSpaceTetPn::initializeAddressingVector(int numElem, std::vector<feInt> &adr) const
{
  for(int j = 0; j < 4; ++j) {
    adr[j] = _numbering->getVertexDOF(_mesh, _cncGeoTag, numElem, j);
  }
  int start = 4;
  if(_n >= 2) {
    // Loop over the edges
    for(int iE = 0; iE < 6; ++iE) {
      int e = _mesh->getEdge(_cncGeoTag, numElem, iE);

      if(e > 0) {
        // Edge orientation is positive: Number edge DOF in default order
        for(int j = 0; j < _n-1; ++j) {
          adr[start + iE * (_n-1) + j] = _numbering->getEdgeDOF(_mesh, _cncGeoTag, numElem, iE, j);
        }
      } else {
        // Edge orientation is negative: Number edge DOF in reverse order
        for(int j = 0; j < _n-1; ++j) {
          adr[start + (iE+1) * (_n-1) - j - 1] = _numbering->getEdgeDOF(_mesh, _cncGeoTag, numElem, iE, j);
        }
      }
    }
  }
  if(_n >= 3) {
    start += 6 * (_n-1);
    int nDOFPerFace = (_n-1)*(_n-2)/2;
    // Loop over the edges
    for(int iF = 0; iF < 4; ++iF) {

      int f = _mesh->getFace(_cncGeoTag, numElem, iF);

      // MODIFY: Rotate when necessary
      if(f > 0) {
        // Face orientation is positive: Number face DOF in default order
        for(int j = 0; j < nDOFPerFace; ++j) {
          adr[start + iF * nDOFPerFace + j] = _numbering->getFaceDOF(_mesh, _cncGeoTag, numElem, iF, j);
        }
      } else {
        // Face orientation is negative: Number face DOF in reverse order
        for(int j = 0; j < nDOFPerFace; ++j) {
          adr[start + (iF+1) * (nDOFPerFace) - j - 1] = _numbering->getFaceDOF(_mesh, _cncGeoTag, numElem, iF, j);
        }
      }
    }
  }
  if(_n >= 4) {
    start += 4 * (_n-1)*(_n-2)/2;
    int nDOFPerElem = _nFunctions - 4 - 6 * (_n-1) - 4 * (_n-1)*(_n-2)/2;
    for(int i = 0; i < nDOFPerElem; ++i) {
      adr[start + i] = _numbering->getElementDOF(_mesh, _cncGeoTag, numElem, i);
    }
  }
}