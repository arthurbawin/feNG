#include "feSpace.h"
#include "feNumber.h"

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
        vectorArray[cnt++] = scalarArray[chunkSize * j + i];
      }
    }
  }
}

// -----------------------------------------------------------------------------
// Lagrange interpolation functions of degree 0 on 0D point
// -----------------------------------------------------------------------------
feSpace0DP0::feSpace0DP0(const std::string &cncGeoID)
  : feScalarSpace(0, nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 1;
  _Lcoor = {1., 0., 0.};
  _isGeometricInterpolant = true;
};

feSpace0DP0::feSpace0DP0(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                         feFunction *fct)
  : feScalarSpace(0, mesh, fieldID, cncGeoID, fct)
{
  _nFunctions = 1;
  _Lcoor = {1., 0., 0.};
};

std::vector<double> feSpace0DP0::L(double *r) { return {1.}; };

void feSpace0DP0::L(double *r, double *L) { L[0] = 1.; };

void feSpace0DP0::initializeNumberingUnknowns()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i)
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 0);
}

void feSpace0DP0::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i)
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 0);
}

void feSpace0DP0::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0);
}

// -----------------------------------------------------------------------------
// 0D Hermite point to fix boundary conditions
// -----------------------------------------------------------------------------
feSpace0D_Hermite::feSpace0D_Hermite(feMesh *mesh, const std::string fieldID,
                                     const std::string cncGeoID, feFunction *fct)
  : feScalarSpace(0, mesh, fieldID, cncGeoID, fct)
{
  _nFunctions = 2;
  _Lcoor = {1., 0., 0., 1., 0., 0.};
};

std::vector<double> feSpace0D_Hermite::L(double *r) { return {1., 1.}; };

void feSpace0D_Hermite::L(double *r, double *L)
{
  L[0] = 1.;
  L[1] = 1.;
};

void feSpace0D_Hermite::initializeNumberingUnknowns()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 0, 2);
  }
}

void feSpace0D_Hermite::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 0);
  }
}

void feSpace0D_Hermite::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 1);
}

// -----------------------------------------------------------------------------
// Lagrange interpolation functions of degree 0 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
feSpace1DP0::feSpace1DP0(const std::string &cncGeoID)
  : feScalarSpace(1, nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 1;
  _Lcoor = {1., 0., 0.};
  _isGeometricInterpolant = true;
};

feSpace1DP0::feSpace1DP0(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                         feFunction *fct)
  : feScalarSpace(1, mesh, fieldID, cncGeoID, fct)
{
  _nFunctions = 1;
  _Lcoor = {1., 0., 0.};
};

std::vector<double> feSpace1DP0::L(double *r) { return {1.}; };

void feSpace1DP0::L(double *r, double *L) { L[0] = 1.; };

void feSpace1DP0::initializeNumberingUnknowns()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i)
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 0);
}

void feSpace1DP0::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i)
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 0);
}

void feSpace1DP0::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0);
}

// -----------------------------------------------------------------------------
// Lagrange interpolation functions of degree 1 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
feSpace1DP1::feSpace1DP1(std::string cncGeoID) : feScalarSpace(1, nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 2;
  _Lcoor = {-1., 0., 0., 1., 0., 0.};
  _isGeometricInterpolant = true;
};

feSpace1DP1::feSpace1DP1(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                         feFunction *fct)
  : feScalarSpace(1, mesh, fieldID, cncGeoID, fct)
{
  _nFunctions = 2;
  _Lcoor = {-1., 0., 0., 1., 0., 0.};
};

std::vector<double> feSpace1DP1::L(double *r) { return {(1. - r[0]) / 2., (1. + r[0]) / 2.}; };

void feSpace1DP1::L(double *r, double *L)
{
  L[0] = (1. - r[0]) / 2.;
  L[1] = (1. + r[0]) / 2.;
};

std::vector<double> feSpace1DP1::dLdr(double *r) { return {-1. / 2., 1. / 2.}; };

void feSpace1DP1::initializeNumberingUnknowns()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 1);
  }
}

void feSpace1DP1::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 1);
  }
}

void feSpace1DP1::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1);
}

// -----------------------------------------------------------------------------
// Vector Lagrange interpolation functions of degree 1 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
template <int dim>
feSpaceVecP1<dim>::feSpaceVecP1(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                                feVectorFunction *fct)
  : feVectorSpace(1, mesh, fieldID, cncGeoID, fct)
{
  _nComponents = dim;
  _nFunctions = 2 * dim;
  double coor[6] = {-1., 0., 0., 1., 0., 0.};
  _Lcoor.resize(6 * dim);
  duplicateScalarArray(6, 3, coor, dim, _Lcoor);
};

template <int dim> std::vector<double> feSpaceVecP1<dim>::L(double *r)
{
  double phi[2] = {(1. - r[0]) / 2., (1. + r[0]) / 2.};
  std::vector<double> res(_nFunctions);
  duplicateScalarArray(2, 1, phi, dim, res);
  return res;
};

template <int dim> void feSpaceVecP1<dim>::L(double *r, double *res)
{
  double phi[2] = {(1. - r[0]) / 2., (1. + r[0]) / 2.};
  duplicateScalarArray(2, 1, phi, dim, res);
};

template <int dim> std::vector<double> feSpaceVecP1<dim>::dLdr(double *r)
{
  double dldr[2] = {-1. / 2., 1. / 2.};
  std::vector<double> res(_nFunctions);
  duplicateScalarArray(2, 1, dldr, dim, res);
  return res;
};

template <int dim> void feSpaceVecP1<dim>::initializeNumberingUnknowns()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 0, dim);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 1, dim);
  }
}

template <int dim> void feSpaceVecP1<dim>::initializeNumberingEssential()
{
  // for(int iComp = 0; iComp < _nComponents; ++iComp){
  //   if(_essentialComponents[iComp]){
  // feInfo("Imposing essential condition on compo %d and space %s - %s",
  //   iComp, _fieldID.data(), _cncGeoID.data());
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 0); //, iComp);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 1); //, iComp);
  }
  //   }
  // }
}

template <int dim>
void feSpaceVecP1<dim>::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 1);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 0);
  adr[3] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 1);
}

template class feSpaceVecP1<1>;
template class feSpaceVecP1<2>;

// -----------------------------------------------------------------------------
// Discontinuous Lagrange interpolation functions of degree 1 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
feSpace1D_DG_P1::feSpace1D_DG_P1(feMesh *mesh, const std::string fieldID,
                                 const std::string cncGeoID, feFunction *fct)
  : feScalarSpace(1, mesh, fieldID, cncGeoID, fct)
{
  _nFunctions = 2;
  _Lcoor = {-1., 0., 0., 1., 0., 0.};
  _DOFinitialization = dofInitialization::NODEWISE;
};

std::vector<double> feSpace1D_DG_P1::L(double *r) { return {(1. - r[0]) / 2., (1. + r[0]) / 2.}; };

void feSpace1D_DG_P1::L(double *r, double *L)
{
  L[0] = (1. - r[0]) / 2.;
  L[1] = (1. + r[0]) / 2.;
};

std::vector<double> feSpace1D_DG_P1::dLdr(double *r) { return {-1. / 2., 1. / 2.}; };

void feSpace1D_DG_P1::initializeNumberingUnknowns()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setUnknownElementDOF(_mesh, _cncGeoID, i, 2);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, 2);
  }
}

void feSpace1D_DG_P1::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialElementDOF(_mesh, _cncGeoID, i);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0);
  }
}

void feSpace1D_DG_P1::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
  adr[1] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 1);
}

// -----------------------------------------------------------------------------
// Crouzeiv_raviart element of degree 0 on line
// -----------------------------------------------------------------------------
feSpace1D_CR0::feSpace1D_CR0(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                             feFunction *fct)
  : feScalarSpace(1, mesh, fieldID, cncGeoID, fct)
{
  _nFunctions = 1;
  _Lcoor = {0., 0., 0.};
};

std::vector<double> feSpace1D_CR0::L(double *r) { return {1.}; };

void feSpace1D_CR0::L(double *r, double *L) { L[0] = 1.; };

void feSpace1D_CR0::initializeNumberingUnknowns()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setUnknownElementDOF(_mesh, _cncGeoID, i, 1);
    // If the line is a boundary element, the edge should be set by the interior element
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, 1);
  }
}

void feSpace1D_CR0::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialElementDOF(_mesh, _cncGeoID, i);
    // Set essential BC on the edge if the line is a boundary element
    // If the line is an interior element, there is nothing in cncGeo->connecEdges
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0);
  }
}

void feSpace1D_CR0::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
}

// -----------------------------------------------------------------------------
// Lagrange interpolation functions of degree 2 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
feSpace1DP2::feSpace1DP2(std::string cncGeoID) : feScalarSpace(1, nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 3;
  _Lcoor = {-1., 0., 0., 1., 0., 0., 0., 0., 0.};
  _isGeometricInterpolant = true;
};

feSpace1DP2::feSpace1DP2(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                         feFunction *fct)
  : feScalarSpace(1, mesh, fieldID, cncGeoID, fct)
{
  _nFunctions = 3;
  _Lcoor = {-1., 0., 0., 1., 0., 0., 0., 0., 0.};
};

std::vector<double> feSpace1DP2::L(double *r)
{
  return {-r[0] * (1. - r[0]) / 2., r[0] * (1. + r[0]) / 2., -(r[0] + 1.) * (r[0] - 1.)};
};

void feSpace1DP2::L(double *r, double *L)
{
  L[0] = -r[0] * (1. - r[0]) / 2.;
  L[1] = r[0] * (1. + r[0]) / 2.;
  L[2] = -(r[0] + 1.) * (r[0] - 1.);
};

std::vector<double> feSpace1DP2::dLdr(double *r)
{
  return {(2. * r[0] - 1.) / 2., (2. * r[0] + 1.) / 2., -2. * r[0]};
};

void feSpace1DP2::initializeNumberingUnknowns()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    // Loop over the elements nodes on geometric interpolant (> 2 if curved)
    for(int j = 0; j < this->getCncGeo()->getNumVerticesPerElem(); ++j) {
      _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, j);
    }
    // If the edge is curved (P2), the middle vertex is already numbered.
    // This is true for Pn geometries with n even, where the middle vertex matches
    // the middle node of the interpolant, but so far we will limit to P2 geometries.
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      _numbering->setUnknownElementDOF(_mesh, _cncGeoID, i, 1);
      // If the line is a boundary element, the edge should be set by the interior element

      // FIXME : When should this be commented ? :/
      _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, 1);
    }
  }
}

void feSpace1DP2::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    for(int j = 0; j < this->getCncGeo()->getNumVerticesPerElem(); ++j) {
      _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, j);
    }
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      _numbering->setEssentialElementDOF(_mesh, _cncGeoID, i);
      // Set essential BC on the edge if the line is a boundary element
      // If the line is an interior element, there is nothing in cncGeo->connecEdges
      // This is only true is there is a surface feSpace associated with the current feSpace !
      _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0);
    }
  }
}

void feSpace1DP2::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1);
  if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() == 2) {
    adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2);
  } else {
    adr[2] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
  }
}

// -----------------------------------------------------------------------------
// Vector Lagrange interpolation functions of degree 2 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
template <int dim>
feSpaceVecP2<dim>::feSpaceVecP2(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                                feVectorFunction *fct)
  : feVectorSpace(1, mesh, fieldID, cncGeoID, fct)
{
  _nComponents = dim;
  _nFunctions = 3 * dim;
  double coor[3 * 3] = {-1., 0., 0., 1., 0., 0., 0., 0., 0.};
  _Lcoor.resize(3 * 3 * dim);
  duplicateScalarArray(3 * 3, 3, coor, dim, _Lcoor);
};

template <int dim> std::vector<double> feSpaceVecP2<dim>::L(double *r)
{
  double phi[3] = {-r[0] * (1. - r[0]) / 2., r[0] * (1. + r[0]) / 2., -(r[0] + 1.) * (r[0] - 1.)};
  std::vector<double> res(_nFunctions);
  duplicateScalarArray(3, 1, phi, dim, res);
  return res;
};

template <int dim> void feSpaceVecP2<dim>::L(double *r, double *res)
{
  double phi[3] = {-r[0] * (1. - r[0]) / 2., r[0] * (1. + r[0]) / 2., -(r[0] + 1.) * (r[0] - 1.)};
  duplicateScalarArray(3, 1, phi, dim, res);
};

template <int dim> std::vector<double> feSpaceVecP2<dim>::dLdr(double *r)
{
  double dldr[3] = {(2. * r[0] - 1.) / 2., (2. * r[0] + 1.) / 2., -2. * r[0]};
  std::vector<double> res(_nFunctions);
  duplicateScalarArray(3, 1, dldr, dim, res);
  return res;
};

template <int dim> void feSpaceVecP2<dim>::initializeNumberingUnknowns()
{
  // See feSpace1DP2 for comments
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    for(int j = 0; j < this->getCncGeo()->getNumVerticesPerElem(); ++j) {
      _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, j, dim);
    }
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      _numbering->setUnknownElementDOF(_mesh, _cncGeoID, i, dim);
      _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, dim);
    }
  }
}

template <int dim> void feSpaceVecP2<dim>::initializeNumberingEssential()
{
  // for(int iComp = 0; iComp < _nComponents; ++iComp){
  //   if(_essentialComponents[iComp]){
  //     feInfo("Imposing essential condition on compo %d and space %s - %s",
  //       iComp, _fieldID.data(), _cncGeoID.data());
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    for(int j = 0; j < this->getCncGeo()->getNumVerticesPerElem(); ++j) {
      _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, j);
    }
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      _numbering->setEssentialElementDOF(_mesh, _cncGeoID, i);
      _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0);
    }
  }
  //   }
  // }
}

template <int dim>
void feSpaceVecP2<dim>::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 1);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 0);
  adr[3] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 1);
  if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() == 2) {
    adr[4] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2, 0);
    adr[5] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 2, 1);
  } else {
    adr[4] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
    adr[5] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 1);
  }
}

template class feSpaceVecP2<1>;
template class feSpaceVecP2<2>;

// -----------------------------------------------------------------------------
// Lagrange interpolation functions of degree 3 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
feSpace1DP3::feSpace1DP3(std::string cncGeoID) : feScalarSpace(1, nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 4;
  _Lcoor = {-1., 0., 0., 1., 0., 0., -1. / 3., 0., 0., 1. / 3., 0., 0.};
  _isGeometricInterpolant = true;
};

feSpace1DP3::feSpace1DP3(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                         feFunction *fct)
  : feScalarSpace(1, mesh, fieldID, cncGeoID, fct)
{
  _nFunctions = 4;
  _Lcoor = {-1., 0., 0., 1., 0., 0., -1. / 3., 0., 0., 1. / 3., 0., 0.};
};

std::vector<double> feSpace1DP3::L(double *r)
{
  return {-9. / 16. * (r[0] + 1. / 3.) * (r[0] - 1. / 3.) * (r[0] - 1.),
          9. / 16. * (r[0] + 1. / 3.) * (r[0] - 1. / 3.) * (r[0] + 1.),
          27. / 16. * (r[0] + 1.) * (r[0] - 1. / 3.) * (r[0] - 1.),
          -27. / 16. * (r[0] + 1.) * (r[0] + 1. / 3.) * (r[0] - 1.)};
};

void feSpace1DP3::L(double *r, double *L)
{
  L[0] = -9. / 16. * (r[0] + 1. / 3.) * (r[0] - 1. / 3.) * (r[0] - 1.);
  L[1] = 9. / 16. * (r[0] + 1. / 3.) * (r[0] - 1. / 3.) * (r[0] + 1.);
  L[2] = 27. / 16. * (r[0] + 1.) * (r[0] - 1. / 3.) * (r[0] - 1.);
  L[3] = -27. / 16. * (r[0] + 1.) * (r[0] + 1. / 3.) * (r[0] - 1.);
};

std::vector<double> feSpace1DP3::dLdr(double *r)
{
  return {r[0] * (9. / 8.) - r[0] * r[0] * (27. / 16.) + 1. / 16.,
          r[0] * (9. / 8.) + r[0] * r[0] * (27. / 16.) - 1. / 16.,
          r[0] * (-9. / 8.) + r[0] * r[0] * (81. / 16.) - 27. / 16.,
          r[0] * (-9. / 8.) - r[0] * r[0] * (81. / 16.) + 27. / 16.};
};

void feSpace1DP3::initializeNumberingUnknowns()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    // Do not loop over all geometric vertices : only 0 and 1
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 1);
    /* The mid-edge vertex does not match any of the element dofs,
    so they must both be added. */
    _numbering->setUnknownElementDOF(_mesh, _cncGeoID, i, 2);
    // TODO : add a test to check if it is in the situation of Verwer
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, 2);
  }
}

void feSpace1DP3::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setEssentialElementDOF(_mesh, _cncGeoID, i);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0);
  }
}

void feSpace1DP3::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1);
  int e = _mesh->getEdge(_cncGeoID, numElem, 0);
  if(e > 0) {
    adr[2] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
    adr[3] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 1);
  } else {
    adr[2] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 1);
    adr[3] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
  }
}

// -----------------------------------------------------------------------------
// Vector Lagrange interpolation functions of degree 3 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
template <int dim>
feSpaceVecP3<dim>::feSpaceVecP3(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                                feVectorFunction *fct)
  : feVectorSpace(1, mesh, fieldID, cncGeoID, fct)
{
  _nComponents = dim;
  _nFunctions = 4 * dim;
  double coor[4 * 3] = {-1., 0., 0., 1., 0., 0., -1. / 3., 0., 0., 1. / 3., 0., 0.};
  _Lcoor.resize(4 * 3 * dim);
  duplicateScalarArray(4 * 3, 3, coor, dim, _Lcoor);
};

template <int dim> std::vector<double> feSpaceVecP3<dim>::L(double *r)
{
  double phi[4] = {-9. / 16. * (r[0] + 1. / 3.) * (r[0] - 1. / 3.) * (r[0] - 1.),
                   9. / 16. * (r[0] + 1. / 3.) * (r[0] - 1. / 3.) * (r[0] + 1.),
                   27. / 16. * (r[0] + 1.) * (r[0] - 1. / 3.) * (r[0] - 1.),
                   -27. / 16. * (r[0] + 1.) * (r[0] + 1. / 3.) * (r[0] - 1.)};
  std::vector<double> res(_nFunctions);
  duplicateScalarArray(4, 1, phi, dim, res);
  return res;
};

template <int dim> void feSpaceVecP3<dim>::L(double *r, double *res)
{
  double phi[4] = {-9. / 16. * (r[0] + 1. / 3.) * (r[0] - 1. / 3.) * (r[0] - 1.),
                   9. / 16. * (r[0] + 1. / 3.) * (r[0] - 1. / 3.) * (r[0] + 1.),
                   27. / 16. * (r[0] + 1.) * (r[0] - 1. / 3.) * (r[0] - 1.),
                   -27. / 16. * (r[0] + 1.) * (r[0] + 1. / 3.) * (r[0] - 1.)};
  duplicateScalarArray(4, 1, phi, dim, res);
};

template <int dim> std::vector<double> feSpaceVecP3<dim>::dLdr(double *r)
{
  double dldr[4] = {r[0] * (9. / 8.) - r[0] * r[0] * (27. / 16.) + 1. / 16.,
                    r[0] * (9. / 8.) + r[0] * r[0] * (27. / 16.) - 1. / 16.,
                    r[0] * (-9. / 8.) + r[0] * r[0] * (81. / 16.) - 27. / 16.,
                    r[0] * (-9. / 8.) - r[0] * r[0] * (81. / 16.) + 27. / 16.};
  std::vector<double> res(_nFunctions);
  duplicateScalarArray(4, 1, dldr, dim, res);
  return res;
};

template <int dim> void feSpaceVecP3<dim>::initializeNumberingUnknowns()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    // Do not loop over all geometric vertices : only 0 and 1
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 0, dim);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 1, dim);
    /* The mid-edge vertex does not match any of the element dofs,
    so they must both be added. */
    _numbering->setUnknownElementDOF(_mesh, _cncGeoID, i, 2 * dim);
    // TODO : add a test to check if it is in the situation of Verwer
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, 2 * dim);
  }
}

template <int dim> void feSpaceVecP3<dim>::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setEssentialElementDOF(_mesh, _cncGeoID, i);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0);
  }
}

template <int dim>
void feSpaceVecP3<dim>::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 1);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 0);
  adr[3] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 1);
  int e = _mesh->getEdge(_cncGeoID, numElem, 0);
  if(e > 0) {
    adr[4] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
    adr[5] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 1);
    adr[6] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 2);
    adr[7] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 3);
  } else {
    adr[4] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 2);
    adr[5] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 3);
    adr[6] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
    adr[7] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 1);
  }
}

template class feSpaceVecP3<1>;
template class feSpaceVecP3<2>;

// -----------------------------------------------------------------------------
// Hermite interpolation functions of degree 3 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
feSpace1D_H3::feSpace1D_H3(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                           feFunction *fct)
  : feScalarSpace(1, mesh, fieldID, cncGeoID, fct)
{
  _nFunctions = 4;
  _Lcoor = {-1., 0., 0., -1., 0., 0., 1., 0., 0., 1., 0., 0.};
};

std::vector<double> feSpace1D_H3::L(double *r)
{
  double x = r[0];
  return {x * (-3.0 / 4.0) + (x * x * x) / 4.0 + 1.0 / 2.0,
          x * (-1.0 / 4.0) - (x * x) / 4.0 + (x * x * x) / 4.0 + 1.0 / 4.0,
          x * (3.0 / 4.0) - (x * x * x) / 4.0 + 1.0 / 2.0,
          x * (-1.0 / 4.0) + (x * x) / 4.0 + (x * x * x) / 4.0 - 1.0 / 4.0};
};

void feSpace1D_H3::L(double *r, double *L)
{
  double x = r[0];
  L[0] = x * (-3.0 / 4.0) + (x * x * x) / 4.0 + 1.0 / 2.0;
  L[1] = x * (-1.0 / 4.0) - (x * x) / 4.0 + (x * x * x) / 4.0 + 1.0 / 4.0;
  L[2] = x * (3.0 / 4.0) - (x * x * x) / 4.0 + 1.0 / 2.0;
  L[3] = x * (-1.0 / 4.0) + (x * x) / 4.0 + (x * x * x) / 4.0 - 1.0 / 4.0;
};

std::vector<double> feSpace1D_H3::dLdr(double *r)
{
  double x = r[0];
  return {(x * x) * (3.0 / 4.0) - 3.0 / 4.0, x * (-1.0 / 2.0) + (x * x) * (3.0 / 4.0) - 1.0 / 4.0,
          (x * x) * (-3.0 / 4.0) + 3.0 / 4.0, x / 2.0 + (x * x) * (3.0 / 4.0) - 1.0 / 4.0};
};

std::vector<double> feSpace1D_H3::d2Ldr2(double *r)
{
  double x = r[0];
  return {x * (3.0 / 2.0), x * (3.0 / 2.0) - 1.0 / 2.0, x * (-3.0 / 2.0),
          x * (3.0 / 2.0) + 1.0 / 2.0};
};

void feSpace1D_H3::initializeNumberingUnknowns()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 0, 2);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 1, 2);
  }
}

void feSpace1D_H3::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 1);
  }
}

void feSpace1D_H3::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 1);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 0);
  adr[3] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 1);
}

// -----------------------------------------------------------------------------
// Lagrange interpolation functions of degree 4 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
feSpace1DP4::feSpace1DP4(std::string cncGeoID) : feScalarSpace(1, nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 5;
  _Lcoor = {-1., 0., 0., 1., 0., 0., -1. / 2., 0., 0., 0., 0., 0., 1. / 2., 0., 0.};
  _isGeometricInterpolant = true;
};

feSpace1DP4::feSpace1DP4(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                         feFunction *fct)
  : feScalarSpace(1, mesh, fieldID, cncGeoID, fct)
{
  _nFunctions = 5;
  _Lcoor = {-1., 0., 0., 1., 0., 0., -1. / 2., 0., 0., 0., 0., 0., 1. / 2., 0., 0.};
};

std::vector<double> feSpace1DP4::L(double *r)
{
  return {r[0] * (r[0] - 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0),
          r[0] * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0),
          r[0] * (r[0] - 1.0) * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (-8.0 / 3.0),
          (r[0] - 1.0) * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * 4.0,
          r[0] * (r[0] - 1.0) * (r[0] + 1.0) * (r[0] + 1.0 / 2.0) * (-8.0 / 3.0)};
};

void feSpace1DP4::L(double *r, double *L)
{
  L[0] = r[0] * (r[0] - 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0);
  L[1] = r[0] * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0);
  L[2] = r[0] * (r[0] - 1.0) * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (-8.0 / 3.0);
  L[3] = (r[0] - 1.0) * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * 4.0;
  L[4] = r[0] * (r[0] - 1.0) * (r[0] + 1.0) * (r[0] + 1.0 / 2.0) * (-8.0 / 3.0);
}

std::vector<double> feSpace1DP4::dLdr(double *r)
{
  return {r[0] * (r[0] - 1.0) * (r[0] - 1.0 / 2.0) * (2.0 / 3.0) +
            r[0] * (r[0] - 1.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0) +
            r[0] * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0) +
            (r[0] - 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0),
          r[0] * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (2.0 / 3.0) +
            r[0] * (r[0] + 1.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0) +
            r[0] * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0) +
            (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0),
          r[0] * (r[0] - 1.0) * (r[0] + 1.0) * (-8.0 / 3.0) -
            r[0] * (r[0] - 1.0) * (r[0] - 1.0 / 2.0) * (8.0 / 3.0) -
            r[0] * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (8.0 / 3.0) -
            (r[0] - 1.0) * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (8.0 / 3.0),
          (r[0] - 1.0) * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * 4.0 +
            (r[0] - 1.0) * (r[0] + 1.0) * (r[0] + 1.0 / 2.0) * 4.0 +
            (r[0] - 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * 4.0 +
            (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * 4.0,
          r[0] * (r[0] - 1.0) * (r[0] + 1.0) * (-8.0 / 3.0) -
            r[0] * (r[0] - 1.0) * (r[0] + 1.0 / 2.0) * (8.0 / 3.0) -
            r[0] * (r[0] + 1.0) * (r[0] + 1.0 / 2.0) * (8.0 / 3.0) -
            (r[0] - 1.0) * (r[0] + 1.0) * (r[0] + 1.0 / 2.0) * (8.0 / 3.0)};
};

void feSpace1DP4::initializeNumberingUnknowns()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    // TODO : Modifier pour elements courbes
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setUnknownElementDOF(_mesh, _cncGeoID, i, 3);
    // TODO : add a test to check if it is in the situation of Verwer
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, 1);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, 2);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, 3);
  }
}

void feSpace1DP4::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    // TODO : Modifier pour elements courbes
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setEssentialElementDOF(_mesh, _cncGeoID, i);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0);
  }
}

void feSpace1DP4::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1);
  int e = _mesh->getEdge(_cncGeoID, numElem, 0);
  if(e > 0) {
    adr[2] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
    adr[3] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 1);
    adr[4] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 2);
  } else {
    adr[4] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
    adr[3] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 1);
    adr[2] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 2);
  }
}

// -----------------------------------------------------------------------------
// Vector Lagrange interpolation functions of degree 4 on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
template <int dim>
feSpaceVecP4<dim>::feSpaceVecP4(feMesh *mesh, const std::string fieldID, const std::string cncGeoID,
                                feVectorFunction *fct)
  : feVectorSpace(1, mesh, fieldID, cncGeoID, fct)
{
  _nComponents = dim;
  _nFunctions = 5 * dim;
  double coor[5 * 3] = {-1., 0., 0., 1., 0., 0., -1. / 2., 0., 0., 0., 0., 0., 1. / 2., 0., 0.};
  _Lcoor.resize(5 * 3 * dim);
  duplicateScalarArray(5 * 3, 3, coor, dim, _Lcoor);
};

template <int dim> std::vector<double> feSpaceVecP4<dim>::L(double *r)
{
  double phi[5] = {r[0] * (r[0] - 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0),
                  r[0] * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0),
                  r[0] * (r[0] - 1.0) * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (-8.0 / 3.0),
                  (r[0] - 1.0) * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * 4.0,
                  r[0] * (r[0] - 1.0) * (r[0] + 1.0) * (r[0] + 1.0 / 2.0) * (-8.0 / 3.0)};
  std::vector<double> res(_nFunctions);
  duplicateScalarArray(5, 1, phi, dim, res);
  return res;
};

template <int dim> void feSpaceVecP4<dim>::L(double *r, double *res)
{
  double phi[4] = {-9. / 16. * (r[0] + 1. / 3.) * (r[0] - 1. / 3.) * (r[0] - 1.),
                   9. / 16. * (r[0] + 1. / 3.) * (r[0] - 1. / 3.) * (r[0] + 1.),
                   27. / 16. * (r[0] + 1.) * (r[0] - 1. / 3.) * (r[0] - 1.),
                   -27. / 16. * (r[0] + 1.) * (r[0] + 1. / 3.) * (r[0] - 1.)};
  duplicateScalarArray(4, 1, phi, dim, res);
};

template <int dim> std::vector<double> feSpaceVecP4<dim>::dLdr(double *r)
{
  double dldr[5] = {r[0] * (r[0] - 1.0) * (r[0] - 1.0 / 2.0) * (2.0 / 3.0) +
            r[0] * (r[0] - 1.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0) +
            r[0] * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0) +
            (r[0] - 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0),
          r[0] * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (2.0 / 3.0) +
            r[0] * (r[0] + 1.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0) +
            r[0] * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0) +
            (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * (2.0 / 3.0),
          r[0] * (r[0] - 1.0) * (r[0] + 1.0) * (-8.0 / 3.0) -
            r[0] * (r[0] - 1.0) * (r[0] - 1.0 / 2.0) * (8.0 / 3.0) -
            r[0] * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (8.0 / 3.0) -
            (r[0] - 1.0) * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (8.0 / 3.0),
          (r[0] - 1.0) * (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * 4.0 +
            (r[0] - 1.0) * (r[0] + 1.0) * (r[0] + 1.0 / 2.0) * 4.0 +
            (r[0] - 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * 4.0 +
            (r[0] + 1.0) * (r[0] - 1.0 / 2.0) * (r[0] + 1.0 / 2.0) * 4.0,
          r[0] * (r[0] - 1.0) * (r[0] + 1.0) * (-8.0 / 3.0) -
            r[0] * (r[0] - 1.0) * (r[0] + 1.0 / 2.0) * (8.0 / 3.0) -
            r[0] * (r[0] + 1.0) * (r[0] + 1.0 / 2.0) * (8.0 / 3.0) -
            (r[0] - 1.0) * (r[0] + 1.0) * (r[0] + 1.0 / 2.0) * (8.0 / 3.0)};
  std::vector<double> res(_nFunctions);
  duplicateScalarArray(5, 1, dldr, dim, res);
  return res;
};

template <int dim> void feSpaceVecP4<dim>::initializeNumberingUnknowns()
{
  int nDOFPerEdge = 3;
  int nDOFPerElem = nDOFPerEdge;
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    // Do not loop over all geometric vertices : only 0 and 1
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 0, dim);
    _numbering->setUnknownVertexDOF(_mesh, _cncGeoID, i, 1, dim);
    /* The mid-edge vertex does not match any of the element dofs,
    so they must both be added. */
    _numbering->setUnknownElementDOF(_mesh, _cncGeoID, i, nDOFPerElem * dim);
    // TODO : add a test to check if it is in the situation of Verwer
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, nDOFPerEdge * dim);
  }
}

template <int dim> void feSpaceVecP4<dim>::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 0);
    _numbering->setEssentialVertexDOF(_mesh, _cncGeoID, i, 1);
    _numbering->setEssentialElementDOF(_mesh, _cncGeoID, i);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0);
  }
}

template <int dim>
void feSpaceVecP4<dim>::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  adr[0] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 0);
  adr[1] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 0, 1);
  adr[2] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 0);
  adr[3] = _numbering->getVertexDOF(_mesh, _cncGeoID, numElem, 1, 1);
  int e = _mesh->getEdge(_cncGeoID, numElem, 0);
  if(e > 0) {
    adr[4] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
    adr[5] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 1);
    adr[6] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 2);
    adr[7] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 3);
    adr[8] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 4);
    adr[9] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 5);
  } else {
    adr[4] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 4);
    adr[5] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 5);
    adr[6] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 2);
    adr[7] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 3);
    adr[8] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 0);
    adr[9] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, 1);
  }
}

template class feSpaceVecP4<1>;
template class feSpaceVecP4<2>;

// -----------------------------------------------------------------------------
// Legendre interpolation functions of arbitrary degree on 1D reference element [-1,1]
// -----------------------------------------------------------------------------
feSpace1D_Legendre::feSpace1D_Legendre(const int degree, feMesh *mesh, const std::string fieldID,
                                       const std::string cncGeoID, feFunction *fct)
  : feScalarSpace(1, mesh, fieldID, cncGeoID, fct)
{
  _degree = degree;
  _nFunctions = degree + 1;
  _Lcoor = {};
  _DOFinitialization = dofInitialization::LEAST_SQUARES;
};

void legendrePolynomials(double x, int p, double *phi)
{
  phi[0] = 1.0;
  if(p >= 1) phi[1] = x;
  for(int i = 1; i < p; i++) phi[i + 1] = ((2 * i + 1) * x * phi[i] - i * phi[i - 1]) / (i + 1);
}

void legendrePolynomialsDerivative(double x, int p, double *phi, double *dphi)
{
  dphi[0] = 0.0;
  if(p >= 1) dphi[1] = 1.0;
  for(int i = 1; i < p; i++)
    dphi[i + 1] = ((2 * i + 1) * (phi[i] + x * dphi[i]) - i * dphi[i - 1]) / (i + 1);
}

std::vector<double> feSpace1D_Legendre::L(double *r)
{
  std::vector<double> phi(_nFunctions);
  legendrePolynomials(r[0], _degree, phi.data());
  return phi;
}

void feSpace1D_Legendre::L(double *r, double *L) { legendrePolynomials(r[0], _degree, L); }

std::vector<double> feSpace1D_Legendre::dLdr(double *r)
{
  // Make sure phi is computed because it's needed for the derivatives
  std::vector<double> phi(_nFunctions);
  std::vector<double> dphi(_nFunctions);
  legendrePolynomials(r[0], _degree, phi.data());
  legendrePolynomialsDerivative(r[0], _degree, phi.data(), dphi.data());
  return dphi;
};

void feSpace1D_Legendre::initializeNumberingUnknowns()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setUnknownElementDOF(_mesh, _cncGeoID, i, _degree + 1);
    _numbering->setUnknownEdgeDOF(_mesh, _cncGeoID, i, 0, _degree + 1);
  }
}

void feSpace1D_Legendre::initializeNumberingEssential()
{
  for(int i = 0; i < _mesh->getNumElements(_cncGeoID); ++i) {
    _numbering->setEssentialElementDOF(_mesh, _cncGeoID, i);
    _numbering->setEssentialEdgeDOF(_mesh, _cncGeoID, i, 0);
  }
}

void feSpace1D_Legendre::initializeAddressingVector(int numElem, std::vector<feInt> &adr)
{
  for(int i = 0; i < _degree + 1; ++i) {
    adr[i] = _numbering->getElementDOF(_mesh, _cncGeoID, numElem, i);
  }
}