#include "feSpace.h"
#include "feMesh.h"
#include "feNumber.h"
#include "fePolynomial.h"

EigenMat I6 = EigenMat::Identity(6, 6);

// feSpace used to interpolate on a geometric connectivity
feSpaceTriP1::feSpaceTriP1(std::string cncGeoID) : feSpace(nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 3;
  _Lcoor = {0., 0., 0., 1., 0., 0., 0., 1., 0.};
}

// feSpace used to compute (bi-)linear forms
feSpaceTriP1::feSpaceTriP1(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct,
                           const bool useGlobalShapeFunctions)
  : feSpace(mesh, fieldID, cncGeoID, fct, useGlobalShapeFunctions)
{
  _nFunctions = 3;
  _Lcoor = {0., 0., 0., 1., 0., 0., 0., 1., 0.};
  // _adr.resize(_nFunctions);
}

std::vector<double> feSpaceTriP1::L(double *r) { return {1.0 - r[0] - r[1], r[0], r[1]}; }
void feSpaceTriP1::L(double *r, double *L)
{
  L[0] = 1.0 - r[0] - r[1];
  L[1] = r[0];
  L[2] = r[1];
}
std::vector<double> feSpaceTriP1::dLdr(double *r) { return {-1.0, 1.0, 0.0}; }
std::vector<double> feSpaceTriP1::dLds(double *r) { return {-1.0, 0.0, 1.0}; }
std::vector<double> feSpaceTriP1::dLdt(double *r) { return {0., 0., 0.}; }

void feSpaceTriP1::initializeNumberingUnknowns(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    number->defDDLSommet(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet(_mesh, _cncGeoID, i, 1);
    number->defDDLSommet(_mesh, _cncGeoID, i, 2);
  }
}

void feSpaceTriP1::initializeNumberingEssential(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 1);
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 2);
  }
}

void feSpaceTriP1::initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr)
{
  adr.resize(3);
  adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
  adr[2] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 2);
}

// feSpace used to interpolate on a geometric connectivity
feSpaceTriP1_nonConsistant::feSpaceTriP1_nonConsistant(std::string cncGeoID)
  : feSpace(nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 3;
  _Lcoor = {0.5, 0., 0., 0.5, 0.5, 0., 0., 0.5, 0.};
}

// feSpace used to compute (bi-)linear forms
feSpaceTriP1_nonConsistant::feSpaceTriP1_nonConsistant(feMesh *mesh, std::string fieldID,
                                                       std::string cncGeoID, feFunction *fct)
  : feSpace(mesh, fieldID, cncGeoID, fct)
{
  _nFunctions = 3;
  _Lcoor = {0.5, 0., 0., 0.5, 0.5, 0., 0., 0.5, 0.};
  // _adr.resize(_nFunctions);
}

std::vector<double> feSpaceTriP1_nonConsistant::L(double *r)
{
  return {1 - 2.0 * r[1], -1.0 + 2.0 * r[0] + 2.0 * r[1], 1 - 2.0 * r[0]};
}
void feSpaceTriP1_nonConsistant::L(double *r, double *L)
{
  L[0] = 1 - 2.0 * r[1];
  L[1] = -1.0 + 2.0 * r[0] + 2.0 * r[1];
  L[2] = 1 - 2.0 * r[0];
}
std::vector<double> feSpaceTriP1_nonConsistant::dLdr(double *r) { return {0.0, 2.0, -2.0}; }
std::vector<double> feSpaceTriP1_nonConsistant::dLds(double *r) { return {-2.0, 2.0, 0.0}; }
std::vector<double> feSpaceTriP1_nonConsistant::dLdt(double *r) { return {0., 0., 0.}; }

void feSpaceTriP1_nonConsistant::initializeNumberingUnknowns(feNumber *number)
{
  int nDOFPerEdge = 1;
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    number->defDDLEdge(_mesh, _cncGeoID, i, 0, nDOFPerEdge);
    number->defDDLEdge(_mesh, _cncGeoID, i, 1, nDOFPerEdge);
    number->defDDLEdge(_mesh, _cncGeoID, i, 2, nDOFPerEdge);
  }
}

void feSpaceTriP1_nonConsistant::initializeNumberingEssential(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 0);
    number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 1);
    number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 2);
  }
}

void feSpaceTriP1_nonConsistant::initializeAddressingVector(feNumber *number, int numElem,
                                                            std::vector<feInt> &adr)
{
  adr[0] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 0, 0);
  adr[1] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 1, 0);
  adr[2] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 2, 0);
}

// feSpace used to interpolate on a geometric connectivity
feSpaceTriP2::feSpaceTriP2(std::string cncGeoID) : feSpace(nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 6;
  _Lcoor = {0., 0., 0., 1., 0., 0., 0., 1., 0., 0.5, 0., 0., 0.5, 0.5, 0., 0., 0.5, 0.};
}

// feSpace used to compute (bi-)linear forms
feSpaceTriP2::feSpaceTriP2(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct,
                           const bool useGlobalShapeFunctions)
  : feSpace(mesh, fieldID, cncGeoID, fct, useGlobalShapeFunctions)
{
  _nFunctions = 6;
  _Lcoor = {0., 0., 0., 1., 0., 0., 0., 1., 0., 0.5, 0., 0., 0.5, 0.5, 0., 0., 0.5, 0.};
  // _adr.resize(_nFunctions);
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

std::vector<double> feSpaceTriP2::dLdt(double *r) { return {0., 0., 0., 0., 0., 0.}; }

void feSpaceTriP2::initializeNumberingUnknowns(feNumber *number)
{
  int nDOFPerEdge = 1;
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    for(int j = 0; j < this->getCncGeo()->getNbNodePerElem(); ++j) {
      number->defDDLSommet(_mesh, _cncGeoID, i, j);
    }
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      number->defDDLEdge(_mesh, _cncGeoID, i, 0, nDOFPerEdge);
      number->defDDLEdge(_mesh, _cncGeoID, i, 1, nDOFPerEdge);
      number->defDDLEdge(_mesh, _cncGeoID, i, 2, nDOFPerEdge);
    }
  }
}

void feSpaceTriP2::initializeNumberingEssential(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    for(int j = 0; j < this->getCncGeo()->getNbNodePerElem(); ++j) {
      number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, j);
    }
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 0);
      number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 1);
      number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 2);
    }
  }
}

void feSpaceTriP2::initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr)
{
  // adr.resize(6);
  adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
  adr[2] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 2);
  if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() == 2) {
    adr[3] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 3);
    adr[4] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 4);
    adr[5] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 5);
  } else {
    adr[3] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 0, 0);
    adr[4] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 1, 0);
    adr[5] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 2, 0);
  }
}

// feSpace used to interpolate on a geometric connectivity
feSpaceTriP2_nonConsistant::feSpaceTriP2_nonConsistant(std::string cncGeoID)
  : feSpace(nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 7;
  _Lcoor = {0., 0.,  0.,  1., 0., 0.,  0., 1.,      0.,      0.5, 0.,
            0., 0.5, 0.5, 0., 0., 0.5, 0., 1. / 3., 1. / 3., 0};
}

// feSpace used to compute (bi-)linear forms
feSpaceTriP2_nonConsistant::feSpaceTriP2_nonConsistant(feMesh *mesh, std::string fieldID,
                                                       std::string cncGeoID, feFunction *fct)
  : feSpace(mesh, fieldID, cncGeoID, fct)
{
  _nFunctions = 7;
  _Lcoor = {0., 0.,  0.,  1., 0., 0.,  0., 1.,    0.,    0.5, 0.,
            0., 0.5, 0.5, 0., 0., 0.5, 0., 1 / 3, 1 / 3, 0};
  // _adr.resize(_nFunctions);
}

std::vector<double> feSpaceTriP2_nonConsistant::L(double *r)
{
  return {(1. - r[0] - r[1]) * (1. - 2. * r[0] - 2. * r[1]),
          r[0] * (2. * r[0] - 1.),
          r[1] * (2. * r[1] - 1.),
          4. * r[0] * (1. - r[0] - r[1]),
          4. * r[0] * r[1],
          4. * r[1] * (1. - r[0] - r[1]),
          2. - 3. * (r[0] * r[0] + r[1] * r[1] + (1. - r[0] - r[1]) * (1. - r[0] - r[1]))};
}
void feSpaceTriP2_nonConsistant::L(double *r, double *L)
{
  L[0] = (1. - r[0] - r[1]) * (1. - 2. * r[0] - 2. * r[1]);
  L[1] = r[0] * (2. * r[0] - 1.);
  L[2] = r[1] * (2. * r[1] - 1.);
  L[3] = 4. * r[0] * (1. - r[0] - r[1]);
  L[4] = 4. * r[0] * r[1];
  L[5] = 4. * r[1] * (1. - r[0] - r[1]);
  L[6] = 2. - 3. * (r[0] * r[0] + r[1] * r[1] + (1. - r[0] - r[1]) * (1. - r[0] - r[1]));
}

std::vector<double> feSpaceTriP2_nonConsistant::dLdr(double *r)
{
  return {4. * (r[0] + r[1]) - 3.,      4. * r[0] - 1., 0.,
          4. * (1. - 2. * r[0] - r[1]), 4. * r[1],      -4. * r[1],
          6. * (1. - r[1] - 2. * r[0])};
}

std::vector<double> feSpaceTriP2_nonConsistant::dLds(double *r)
{
  return {4. * (r[0] + r[1]) - 3.,
          0.,
          4. * r[1] - 1.,
          -4. * r[0],
          4. * r[0],
          4. * (1. - r[0] - 2. * r[1]),
          6. * (1. - r[0] - 2. * r[1])};
}

std::vector<double> feSpaceTriP2_nonConsistant::dLdt(double *r)
{
  return {0., 0., 0., 0., 0., 0., 0.};
}

void feSpaceTriP2_nonConsistant::initializeNumberingUnknowns(feNumber *number)
{
  int nDOFPerEdge = 1;
  int nDOFPerElem = 1;
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    for(int j = 0; j < this->getCncGeo()->getNbNodePerElem(); ++j) {
      number->defDDLSommet(_mesh, _cncGeoID, i, j);
    }
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      number->defDDLEdge(_mesh, _cncGeoID, i, 0, nDOFPerEdge);
      number->defDDLEdge(_mesh, _cncGeoID, i, 1, nDOFPerEdge);
      number->defDDLEdge(_mesh, _cncGeoID, i, 2, nDOFPerEdge);
      number->defDDLElement(_mesh, _cncGeoID, i, nDOFPerElem);
    }
  }
}

void feSpaceTriP2_nonConsistant::initializeNumberingEssential(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    for(int j = 0; j < this->getCncGeo()->getNbNodePerElem(); ++j) {
      number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, j);
    }
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 0);
      number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 1);
      number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 2);
      number->defDDLElement_essentialBC(_mesh, _cncGeoID, i);
    }
  }
}

void feSpaceTriP2_nonConsistant::initializeAddressingVector(feNumber *number, int numElem,
                                                            std::vector<feInt> &adr)
{
  adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
  adr[2] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 2);
  if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() == 2) {
    adr[3] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 3);
    adr[4] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 4);
    adr[5] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 5);
  } else {
    adr[3] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 0, 0);
    adr[4] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 1, 0);
    adr[5] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 2, 0);
  }
  adr[6] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
}

// feSpace used to interpolate on a geometric connectivity
feSpaceTriP3::feSpaceTriP3(std::string cncGeoID) : feSpace(nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 10;
  _Lcoor = {0., 0., 0.,      1., 0., 0.,      0.,      1.,      0.,      1. / 3.,
            0., 0., 2. / 3., 0., 0., 2. / 3., 1. / 3., 0.,      1. / 3., 2. / 3.,
            0., 0., 2. / 3., 0., 0., 1. / 3., 0.,      1. / 3., 1. / 3., 0.};
}

// feSpace used to compute (bi-)linear forms
feSpaceTriP3::feSpaceTriP3(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct,
                           const bool useGlobalShapeFunctions)
  : feSpace(mesh, fieldID, cncGeoID, fct, useGlobalShapeFunctions)
{
  _nFunctions = 10;
  _Lcoor = {0., 0., 0.,      1., 0., 0.,      0.,      1.,      0.,      1. / 3.,
            0., 0., 2. / 3., 0., 0., 2. / 3., 1. / 3., 0.,      1. / 3., 2. / 3.,
            0., 0., 2. / 3., 0., 0., 1. / 3., 0.,      1. / 3., 1. / 3., 0.};
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
        (R * R) * S * (2.7E1 / 2.0) + (R * R) * 9.0 - (R * R * R) * (9.0 / 2.0) +
        (S * S) * 9.0 - (S * S * S) * (9.0 / 2.0) + 1.0;
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

std::vector<double> feSpaceTriP3::dLdt(double *r)
{
  return {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
}

void feSpaceTriP3::initializeNumberingUnknowns(feNumber *number)
{
  int nDOFPerEdge = 2;
  int nDOFPerElem = 1;
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    number->defDDLSommet(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet(_mesh, _cncGeoID, i, 1);
    number->defDDLSommet(_mesh, _cncGeoID, i, 2);
    number->defDDLElement(_mesh, _cncGeoID, i, nDOFPerElem);
    number->defDDLEdge(_mesh, _cncGeoID, i, 0, nDOFPerEdge);
    number->defDDLEdge(_mesh, _cncGeoID, i, 1, nDOFPerEdge);
    number->defDDLEdge(_mesh, _cncGeoID, i, 2, nDOFPerEdge);
  }
}

void feSpaceTriP3::initializeNumberingEssential(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 1);
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 2);
    number->defDDLElement_essentialBC(_mesh, _cncGeoID, i);
    number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 0);
    number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 1);
    number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 2);
  }
}

void feSpaceTriP3::initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr)
{
  adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
  adr[2] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 2);
  int e0 = _mesh->getEdge(_cncGeoID, numElem, 0);
  if(e0 > 0) {
    adr[3] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 0, 0);
    adr[4] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 0, 1);
  } else {
    adr[3] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 0, 1);
    adr[4] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 0, 0);
  }
  int e1 = _mesh->getEdge(_cncGeoID, numElem, 1);
  if(e1 > 0) {
    adr[5] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 1, 0);
    adr[6] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 1, 1);
  } else {
    adr[5] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 1, 1);
    adr[6] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 1, 0);
  }
  int e2 = _mesh->getEdge(_cncGeoID, numElem, 2);
  if(e2 > 0) {
    adr[7] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 2, 0);
    adr[8] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 2, 1);
  } else {
    adr[7] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 2, 1);
    adr[8] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 2, 0);
  }
  adr[9] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
}

/* P4 Lagrange interpolant on triangle :

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

// feSpace used to interpolate on a geometric connectivity
feSpaceTriP4::feSpaceTriP4(std::string cncGeoID) : feSpace(nullptr, "GEO", cncGeoID, nullptr)
{
  _nFunctions = 15;
  double x0 = 0.0, x1 = 1. / 4., x2 = 1. / 2., x3 = 3. / 4., x4 = 1.;
  double y0 = 0.0, y1 = 1. / 4., y2 = 1. / 2., y3 = 3. / 4., y4 = 1.;
  _Lcoor = {x0, y0, 0., x4, y0, 0., x0, y4, 0., x1, y0, 0., x2, y0, 0.,
            x3, y0, 0., x3, y1, 0., x2, y2, 0., x1, y3, 0., y3, x0, 0.,
            y2, x0, 0., y1, x0, 0., x1, y1, 0., x2, y1, 0., x1, y2, 0.};
}

// feSpace used to compute (bi-)linear forms
feSpaceTriP4::feSpaceTriP4(feMesh *mesh, const std::string fieldID, const std::string cncGeoID, feFunction *fct,
                           const bool useGlobalShapeFunctions)
  : feSpace(mesh, fieldID, cncGeoID, fct, useGlobalShapeFunctions)
{
  _nFunctions = 15;
  double x0 = 0.0, x1 = 1. / 4., x2 = 1. / 2., x3 = 3. / 4., x4 = 1.;
  double y0 = 0.0, y1 = 1. / 4., y2 = 1. / 2., y3 = 3. / 4., y4 = 1.;
  _Lcoor = {x0, y0, 0., x4, y0, 0., x0, y4, 0., x1, y0, 0., x2, y0, 0.,
            x3, y0, 0., x3, y1, 0., x2, y2, 0., x1, y3, 0., y3, x0, 0.,
            y2, x0, 0., y1, x0, 0., x1, y1, 0., x2, y1, 0., x1, y2, 0.};
  // _adr.resize(_nFunctions);
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

std::vector<double> feSpaceTriP4::dLdr(double *r)
{
  feWarning("BASIS FUNCTIONS MAY BE WRONG : There is a 0 shape function...");
  double R = r[0];
  double S = r[1];
  return {32. / 3. *
            (4. * R * R * R + 12. * R * R * S - (15. * R * R) / 2. + 12. * R * S * S - 15. * R * S +
             (35. * R) / 8. + 4. * S * S * S - (15. * S * S) / 2. + (35. * S) / 8. - 25. / 32.),
          32. / 3. * (4. * R * R * R - (9. * R * R) / 2. + (11. * R) / 8. - 3. / 32.),
          32. / 3. * (0.0),
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
void feSpaceTriP4::L(double *r, double *L)
{
  double R = r[0];
  double S = r[1];
  L[0] = 32. / 3. * (4. * R * R * R + 12. * R * R * S - (15. * R * R) / 2. + 12. * R * S * S - 15. * R * S +
             (35. * R) / 8. + 4. * S * S * S - (15. * S * S) / 2. + (35. * S) / 8. - 25. / 32.);
  L[1] = 32. / 3. * (4. * R * R * R - (9. * R * R) / 2. + (11. * R) / 8. - 3. / 32.);
  L[2] = 32. / 3. * (0.0),
  L[3] = 128. / 3. * (-4. * R * R * R - 9. * R * R * S + (27. * R * R) / 4. - 6. * R * S * S + 9. * R * S -
             (13. * R) / 4. - S * S * S + (9. * S * S) / 4. - (13. * S) / 8. + 3. / 8.);
  L[4] = 64. * (4. * R * R * R + 6. * R * R * S - 6. * R * R + 2. * R * S * S - (9. * R * S) / 2. +
                 (19. * R) / 8. - S * S / 4. + (7. * S) / 16. - 3. / 16.);
  L[5] = 128. / 3. * ((3. * R * S) / 2. - S / 8. - (7. * R) / 4. - 3. * R * R * S + (21. * R * R) / 4. -
             4. * R * R * R + 1. / 8.);
  L[6] = 128. / 3. * ((S * (24. * R * R - 12. * R + 1.)) / 8.);
  L[7] = 64. * ((S * (8. * R - 1.) * (4. * S - 1.)) / 16.);
  L[8] = 128. / 3. * (S * (S - 1. / 2.) * (S - 1. / 4.));
  L[9] = 128. / 3. * (-S * (S - 1. / 2.) * (S - 1. / 4.));
  L[10] = 64. * ((S * (4. * S - 1.) * (8. * R + 8. * S - 7.)) / 16.);
  L[11] = 128. / 3. * (-(S * (24. * R * R + 48. * R * S - 36. * R + 24. * S * S - 36. * S + 13.)) / 8.);
  L[12] = 128. * ((S * (12. * R * R + 16. * R * S - 14. * R + 4. * S * S - 7. * S + 3.)) / 4.);
  L[13] = 128. * (-(S * (8. * R * S - S - 10. * R + 12. * R * R + 1.)) / 4.);
  L[14] = 128. * (-(S * (4. * S - 1.) * (2. * R + S - 1.)) / 4.);
}

std::vector<double> feSpaceTriP4::dLds(double *r)
{
  double R = r[0];
  double S = r[1];
  return {32. / 3. *
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
}

std::vector<double> feSpaceTriP4::dLdt(double *r)
{
  return {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
}

void feSpaceTriP4::initializeNumberingUnknowns(feNumber *number)
{
  int nDOFPerEdge = 3;
  int nDOFPerElem = 3;
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    number->defDDLSommet(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet(_mesh, _cncGeoID, i, 1);
    number->defDDLSommet(_mesh, _cncGeoID, i, 2);
    number->defDDLElement(_mesh, _cncGeoID, i, nDOFPerElem);
    number->defDDLEdge(_mesh, _cncGeoID, i, 0, nDOFPerEdge);
    number->defDDLEdge(_mesh, _cncGeoID, i, 1, nDOFPerEdge);
    number->defDDLEdge(_mesh, _cncGeoID, i, 2, nDOFPerEdge);
  }
}

void feSpaceTriP4::initializeNumberingEssential(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 1);
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 2);
    number->defDDLElement_essentialBC(_mesh, _cncGeoID, i);
    number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 0);
    number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 1);
    number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 2);
  }
}

void feSpaceTriP4::initializeAddressingVector(feNumber *number, int numElem, std::vector<feInt> &adr)
{
  adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
  adr[2] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 2);
  int e0 = _mesh->getEdge(_cncGeoID, numElem, 0);
  if(e0 > 0) {
    adr[3] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 0, 0);
    adr[4] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 0, 1);
    adr[5] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 0, 2);
  } else {
    adr[3] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 0, 2);
    adr[4] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 0, 1);
    adr[5] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 0, 0);
  }
  int e1 = _mesh->getEdge(_cncGeoID, numElem, 1);
  if(e1 > 0) {
    adr[6] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 1, 0);
    adr[7] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 1, 1);
    adr[8] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 1, 2);
  } else {
    adr[6] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 1, 2);
    adr[7] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 1, 1);
    adr[8] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 1, 0);
  }
  int e2 = _mesh->getEdge(_cncGeoID, numElem, 2);
  if(e2 > 0) {
    adr[9] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 2, 0);
    adr[10] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 2, 1);
    adr[11] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 2, 2);
  } else {
    adr[9] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 2, 2);
    adr[10] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 2, 1);
    adr[11] = number->getDDLEdge(_mesh, _cncGeoID, numElem, 2, 0);
  }
  adr[12] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
  adr[13] = number->getDDLElement(_mesh, _cncGeoID, numElem, 1);
  adr[14] = number->getDDLElement(_mesh, _cncGeoID, numElem, 2);
}