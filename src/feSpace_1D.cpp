#include "feSpace.h"
// #include "feMesh.h"
#include "feNumber.h"

//
// Lagrange interpolation functions of degree 0 on 1D reference element [-1,1]
//
feSpace1DP0::feSpace1DP0(const std::string &cncGeoID)
  : feSpace(nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 1;
  _nQuad = 1;
  _Lcoor = {1., 0., 0.};
};

feSpace1DP0::feSpace1DP0(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct)
  : feSpace(mesh, fieldID, cncGeoID, fct)
{
  _nFunctions = 1;
  _nQuad = 1;
  _Lcoor = {1., 0., 0.};
};

std::vector<double> feSpace1DP0::L(double *r)
{
  return {1.};
};

void feSpace1DP0::L(double *r, double *L)
{ 
  L[0] = 1.;
};

std::vector<double> feSpace1DP0::dLdr(double *r)
{ 
  return {0.};
};

std::vector<double> feSpace1DP0::dLds(double *r)
{ 
  return {0.};
};

std::vector<double> feSpace1DP0::dLdt(double *r)
{ 
  return {0.};
};

feStatus feSpace1DP0::Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
  std::vector<double> &dLdx, std::vector<double> &dLdy)
{
  printf("Not implemented\n");
  exit(-1);
};

void feSpace1DP0::initializeNumberingUnknowns(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) number->defDDLSommet(_mesh, _cncGeoID, i, 0);
}

void feSpace1DP0::initializeNumberingEssential(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i)
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 0);
}

void feSpace1DP0::initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr)
{
  adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
}

//
// Lagrange interpolation functions of degree 1 on 1D reference element [-1,1]
//
feSpace1DP1::feSpace1DP1(std::string cncGeoID) : feSpace(nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 2;
  _Lcoor = {-1., 0., 0., 1., 0., 0.};
};

feSpace1DP1::feSpace1DP1(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct)
  : feSpace(mesh, fieldID, cncGeoID, fct)
{
  _nFunctions = 2;
  _Lcoor = {-1., 0., 0., 1., 0., 0.};
};

std::vector<double> feSpace1DP1::L(double *r)
{ 
  return {(1. - r[0]) / 2., (1. + r[0]) / 2.};
};

void feSpace1DP1::L(double *r, double *L)
{
  L[0] = (1. - r[0]) / 2.;
  L[1] = (1. + r[0]) / 2.;
};

std::vector<double> feSpace1DP1::dLdr(double *r)
{
  return {-1. / 2., 1. / 2.};
};

std::vector<double> feSpace1DP1::dLds(double *r)
{
  return {0., 0.};
};

std::vector<double> feSpace1DP1::dLdt(double *r)
{
  return {0., 0.};
};

feStatus feSpace1DP1::Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                       std::vector<double> &dLdx, std::vector<double> &dLdy)
{
  printf("Not implemented\n");
  exit(-1);
};

void feSpace1DP1::initializeNumberingUnknowns(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    number->defDDLSommet(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet(_mesh, _cncGeoID, i, 1);
  }
}

void feSpace1DP1::initializeNumberingEssential(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 1);
  }
}

void feSpace1DP1::initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr)
{
  adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
}

//
// Non-conformal interpolation functions of degree 1
//
feSpace1DP1_nonConsistant::feSpace1DP1_nonConsistant(std::string cncGeoID)
  : feSpace(nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 1;
  _Lcoor = {0., 0., 0.};
};

feSpace1DP1_nonConsistant::feSpace1DP1_nonConsistant(feMesh *mesh, 
  const std::string fieldID, const std::string cncGeoID, feFunction *fct)
  : feSpace(mesh, fieldID, cncGeoID, fct)
{
  _nFunctions = 1;
  _Lcoor = {0., 0., 0.};
};

std::vector<double> feSpace1DP1_nonConsistant::L(double *r)
{
  return {1.};
};

void feSpace1DP1_nonConsistant::L(double *r, double *L)
{
  L[0] = 1.;
};

std::vector<double> feSpace1DP1_nonConsistant::dLdr(double *r)
{
  return {0.};
};

std::vector<double> feSpace1DP1_nonConsistant::dLds(double *r)
{
  return {0.};
};

std::vector<double> feSpace1DP1_nonConsistant::dLdt(double *r)
{
  return {0.};
};

feStatus feSpace1DP1_nonConsistant::Lphys(int iElm, std::vector<double> &x,
  std::vector<double> &L, std::vector<double> &dLdx, std::vector<double> &dLdy)
{
  printf("Not implemented\n");
  exit(-1);
};

void feSpace1DP1_nonConsistant::initializeNumberingUnknowns(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    number->defDDLElement(_mesh, _cncGeoID, i, 1);
    // If the line is a boundary element, the edge should be set by the interior element
    number->defDDLEdge(_mesh, _cncGeoID, i, 0, 1);
  }
}

void feSpace1DP1_nonConsistant::initializeNumberingEssential(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    number->defDDLElement_essentialBC(_mesh, _cncGeoID, i);
    // Set essential BC on the edge if the line is a boundary element
    // If the line is an interior element, there is nothing in cncGeo->connecEdges
    number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 0);
  }
}

void feSpace1DP1_nonConsistant::initializeAddressingVector(feNumber *number, int numElem,
                                                           std::vector<feInt> &adr)
{
  adr[0] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
}

//
// Lagrange interpolation functions of degree 2 on 1D reference element [-1,1]
//
feSpace1DP2::feSpace1DP2(std::string cncGeoID)
  : feSpace(nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 3;
  _Lcoor = {-1., 0., 0., 1., 0., 0., 0., 0., 0.};
};

feSpace1DP2::feSpace1DP2(feMesh *mesh, const std::string fieldID,
  const std::string cncGeoID, feFunction *fct)
  : feSpace(mesh, fieldID, cncGeoID, fct)
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
  L[1] =  r[0] * (1. + r[0]) / 2.;
  L[2] = -(r[0] + 1.) * (r[0] - 1.);
};

std::vector<double> feSpace1DP2::dLdr(double *r)
{
  return {(2. * r[0] - 1.) / 2., (2. * r[0] + 1.) / 2., -2. * r[0]};
};

std::vector<double> feSpace1DP2::dLds(double *r)
{ 
  return {0., 0., 0.};
};

std::vector<double> feSpace1DP2::dLdt(double *r)
{ 
  return {0., 0., 0.};
};

feStatus feSpace1DP2::Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                       std::vector<double> &dLdx, std::vector<double> &dLdy)
{
  printf("Not implemented\n");
  exit(-1);
};

void feSpace1DP2::initializeNumberingUnknowns(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    // Loop over the elements nodes on geometric interpolant (> 2 if curved)
    for(int j = 0; j < this->getCncGeo()->getNbNodePerElem(); ++j) {
      number->defDDLSommet(_mesh, _cncGeoID, i, j);
    }
    // If the edge is curved (P2), the middle vertex is already numbered.
    // This is true for Pn geometries with n even, where the middle vertex matches
    // the middle node of the interpolant, but so far we will limit to P2 geometries.
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      number->defDDLElement(_mesh, _cncGeoID, i, 1);
      // If the line is a boundary element, the edge should be set by the interior element

      // FIXME : When should this be commented ? :/
      number->defDDLEdge(_mesh, _cncGeoID, i, 0, 1);
    }
  }
}

void feSpace1DP2::initializeNumberingEssential(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    for(int j = 0; j < this->getCncGeo()->getNbNodePerElem(); ++j) {
      number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, j);
    }
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      number->defDDLElement_essentialBC(_mesh, _cncGeoID, i);
      // Set essential BC on the edge if the line is a boundary element
      // If the line is an interior element, there is nothing in cncGeo->connecEdges
      // This is only true is there is a surface feSpace associated with the current feSpace !
      number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 0);
    }
  }
}

void feSpace1DP2::initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr)
{
  adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
  if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() == 2) {
    adr[2] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 2);
  } else {
    adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
  }
}

//
// Lagrange interpolation functions of degree 3 on 1D reference element [-1,1]
//
feSpace1DP3::feSpace1DP3(std::string cncGeoID)
  : feSpace(nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 4;
  _Lcoor = {-1., 0., 0., 1., 0., 0., -1. / 3., 0., 0., 1. / 3., 0., 0.};
};

feSpace1DP3::feSpace1DP3(feMesh *mesh, const std::string fieldID,
  const std::string cncGeoID, feFunction *fct)
  : feSpace(mesh, fieldID, cncGeoID, fct)
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
  L[0] =  -9. / 16. * (r[0] + 1. / 3.) * (r[0] - 1. / 3.) * (r[0] - 1.);
  L[1] =   9. / 16. * (r[0] + 1. / 3.) * (r[0] - 1. / 3.) * (r[0] + 1.);
  L[2] =  27. / 16. * (r[0] + 1.) * (r[0] - 1. / 3.) * (r[0] - 1.);
  L[3] = -27. / 16. * (r[0] + 1.) * (r[0] + 1. / 3.) * (r[0] - 1.);
};

std::vector<double> feSpace1DP3::dLdr(double *r)
{
  return {r[0] * (9. / 8.) - r[0] * r[0] * (27. / 16.) + 1. / 16.,
          r[0] * (9. / 8.) + r[0] * r[0] * (27. / 16.) - 1. / 16.,
          r[0] * (-9. / 8.) + r[0] * r[0] * (81. / 16.) - 27. / 16.,
          r[0] * (-9. / 8.) - r[0] * r[0] * (81. / 16.) + 27. / 16.};
};

std::vector<double> feSpace1DP3::dLds(double *r)
{
  return {0., 0., 0., 0.};
};

std::vector<double> feSpace1DP3::dLdt(double *r)
{
  return {0., 0., 0., 0.};
};

feStatus feSpace1DP3::Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                       std::vector<double> &dLdx, std::vector<double> &dLdy)
{
  printf("Not implemented\n");
  exit(-1);
};

void feSpace1DP3::initializeNumberingUnknowns(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    // Do not loop over all geometric vertices : only 0 and 1
    number->defDDLSommet(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet(_mesh, _cncGeoID, i, 1);
    /* The mid-edge vertex does not match any of the element dofs,
    so they must both be added. */
    number->defDDLElement(_mesh, _cncGeoID, i, 2);
    // TODO : add a test to check if it is in the situation of Verwer
    number->defDDLEdge(_mesh, _cncGeoID, i, 0, 2);
  }
}

void feSpace1DP3::initializeNumberingEssential(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 1);
    number->defDDLElement_essentialBC(_mesh, _cncGeoID, i);
    number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 0);
  }
}

void feSpace1DP3::initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr)
{
  adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
  int e = _mesh->getEdge(_cncGeoID, numElem, 0);
  if(e > 0) {
    adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
    adr[3] = number->getDDLElement(_mesh, _cncGeoID, numElem, 1);
  } else {
    adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 1);
    adr[3] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
  }
}

//
// Lagrange interpolation functions of degree 4 on 1D reference element [-1,1]
//
feSpace1DP4::feSpace1DP4(std::string cncGeoID)
  : feSpace(nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 5;
  _Lcoor = {-1., 0., 0., 1., 0., 0., -1. / 2., 0., 0., 0., 0., 0., 1. / 2., 0., 0.};
};

feSpace1DP4::feSpace1DP4(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct)
  : feSpace(mesh, fieldID, cncGeoID, fct)
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

std::vector<double> feSpace1DP4::dLds(double *r)
{
  return {0., 0., 0., 0., 0.};
};

std::vector<double> feSpace1DP4::dLdt(double *r)
{
  return {0., 0., 0., 0., 0.};
};

feStatus feSpace1DP4::Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                       std::vector<double> &dLdx, std::vector<double> &dLdy)
{
  printf("Not implemented\n");
  exit(-1);
};

void feSpace1DP4::initializeNumberingUnknowns(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    // TODO : Modifier pour elements courbes
    number->defDDLSommet(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet(_mesh, _cncGeoID, i, 1);
    number->defDDLElement(_mesh, _cncGeoID, i, 3);
    // TODO : add a test to check if it is in the situation of Verwer
    number->defDDLEdge(_mesh, _cncGeoID, i, 0, 1);
    number->defDDLEdge(_mesh, _cncGeoID, i, 0, 2);
    number->defDDLEdge(_mesh, _cncGeoID, i, 0, 3);
  }
}

void feSpace1DP4::initializeNumberingEssential(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    // TODO : Modifier pour elements courbes
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 1);
    number->defDDLElement_essentialBC(_mesh, _cncGeoID, i);
    number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 0);
  }
}

void feSpace1DP4::initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr)
{
  adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
  int e = _mesh->getEdge(_cncGeoID, numElem, 0);
  if(e > 0) {
    adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
    adr[3] = number->getDDLElement(_mesh, _cncGeoID, numElem, 1);
    adr[4] = number->getDDLElement(_mesh, _cncGeoID, numElem, 2);
  } else {
    adr[4] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
    adr[3] = number->getDDLElement(_mesh, _cncGeoID, numElem, 1);
    adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 2);
  }
}

//
// Legendre interpolation functions of arbitrary degree on 1D reference element [-1,1]
//
feSpace1D_Legendre::feSpace1D_Legendre(const int degree, feMesh *mesh,
  const std::string fieldID, const std::string cncGeoID, feFunction *fct)
  : feSpace(mesh, fieldID, cncGeoID, fct)
{
  _degree = degree;
  _nFunctions = degree + 1;
  _Lcoor = {};
};

feStatus feSpace1D_Legendre::Lphys(int iElm, std::vector<double> &x, std::vector<double> &L,
                         std::vector<double> &dLdx, std::vector<double> &dLdy)
{
  printf("Not implemented\n");
  exit(-1);
};

void legendrePolynomials(double x, int p, double *phi)
{
  phi[0] = 1.0;
  if (p >= 1) phi[1] = x;
  for (int i=1; i < p; i++)
    phi[i+1] = ((2*i+1)*x*phi[i] - i*phi[i-1]) / (i+1); 
}

std::vector<double> feSpace1D_Legendre::L(double *r)
{
  std::vector<double> phi(_nFunctions);
  legendrePolynomials(r[0], _degree, phi.data());
  return phi;
}

void feSpace1D_Legendre::L(double *r, double *L)
{
  legendrePolynomials(r[0], _degree, L);
}

std::vector<double> feSpace1D_Legendre::dLdr(double *r)
{
  // FIXME
  printf("Implement derivatives for arbitray degree\n");
  exit(-1);
  return {0., 1., 3.*r[0]};
};

std::vector<double> feSpace1D_Legendre::dLds(double *r)
{
  return {0., 0.};
};

std::vector<double> feSpace1D_Legendre::dLdt(double *r)
{
  return {0., 0.};
};

void feSpace1D_Legendre::initializeNumberingUnknowns(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    // number->defDDLSommet(_mesh, _cncGeoID, i, 0, 2);
    // number->defDDLSommet(_mesh, _cncGeoID, i, 1, 2);
    // number->defDDLElement(_mesh, _cncGeoID, i, _degree - 1);
    number->defDDLElement(_mesh, _cncGeoID, i, _degree + 1);
    number->defDDLEdge(_mesh, _cncGeoID, i, 0, _degree + 1);

    // Autre option: 2 aux sommets et le reste sur l'élément?
  }
}

void feSpace1D_Legendre::initializeNumberingEssential(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    // number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 0);
    // number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 1);
    number->defDDLElement_essentialBC(_mesh, _cncGeoID, i);
  }
}

void feSpace1D_Legendre::initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr)
{
  // adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0, 1); // DOF^+ (right) on first vertex
  // adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1, 0); // DOF^- (left)  on end   vertex
  // for(int i = 0; i < _degree - 1; ++i){
  //   adr[i+2] = number->getDDLElement(_mesh, _cncGeoID, numElem, i);
  // }
  for(int i = 0; i < _degree + 1; ++i){
    adr[i] = number->getDDLElement(_mesh, _cncGeoID, numElem, i);
  }
  // adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0, 1); // DOF^+ (right) on first vertex
  // adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1, 0); // DOF^- (left)  on end   vertex
}