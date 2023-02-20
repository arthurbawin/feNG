
#include "fePolynomial.h"

Polynomial::Polynomial(int deg, std::vector<double> &coeff) : _deg(deg), _nMonomials(coeff.size())
{
  _c.resize(coeff.size());
  for(unsigned i = 0; i < coeff.size(); ++i) _c[i] = coeff[i];

  _exponentsX.resize(_nMonomials);
  _exponentsY.resize(_nMonomials);

  int cnt = 0;
  for(int i = 0; i < deg + 1; ++i) {
    for(int j = i; j >= 0; --j) {
      _exponentsX[cnt] = j;
      ++cnt;
    }
  }
  cnt = 0;
  for(int i = 0; i < deg + 1; ++i) {
    for(int j = 0; j <= i; ++j) {
      _exponentsY[cnt] = j;
      ++cnt;
    }
  }
}

Polynomial::Polynomial(const Polynomial &p) : _deg(p.deg()), _nMonomials(p.nMonomials())
{
  _c.resize(_nMonomials);
  for(int i = 0; i < _nMonomials; ++i) _c[i] = p.coeff(i);
  _exponentsX.resize(_nMonomials);
  _exponentsY.resize(_nMonomials);
  int cnt = 0;
  for(int i = 0; i < _deg + 1; ++i) {
    for(int j = i; j >= 0; --j) {
      _exponentsX[cnt] = j;
      ++cnt;
    }
  }
  cnt = 0;
  for(int i = 0; i < _deg + 1; ++i) {
    for(int j = 0; j <= i; ++j) {
      _exponentsY[cnt] = j;
      ++cnt;
    }
  }
}

double Polynomial::eval(std::vector<double> &x)
{
  double val, res = 0.0;
  for(int i = 0; i < _nMonomials; ++i) {
    val = _c[i];
    for(int ix = 0; ix < _exponentsX[i]; ++ix) val *= x[0];
    for(int iy = 0; iy < _exponentsY[i]; ++iy) val *= x[1];
    res += val;
  }
  return res;
}

double Polynomial::integrate(feMesh *mesh, std::string cncGeoID, int iElm)
{
  double res = 0.0;

  std::vector<double> geoCoord;
  mesh->getCoord(cncGeoID, iElm, geoCoord);
  feSpace *geoSpace = mesh->getGeometricSpace(cncGeoID);

  int n = geoSpace->getNumQuadPoints();
  const std::vector<double> &J = geoSpace->getCncGeo()->getJacobians();
  std::vector<double> &w = geoSpace->getQuadratureWeights();
  std::vector<double> x(3, 0.0);
  for(int i = 0; i < n; ++i) {
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, i, x);
    res += J[n * iElm + i] * w[i] * this->eval(x);
  }
  return res;
}

double Polynomial::innerProduct(Polynomial &other, feMesh *mesh, std::string cncGeoID, int iElm)
{
  double res = 0.0;

  std::vector<double> geoCoord;
  mesh->getCoord(cncGeoID, iElm, geoCoord);
  feSpace *geoSpace = mesh->getGeometricSpace(cncGeoID);

  // std::cout<<cncGeoID<<" - "<<iElm<<std::endl;
  // std::cout<<geoCoord.size()<<" - "<<iElm<<std::endl;
  // std::cout<<geoSpace->getCncGeoID()<<" - "<<iElm<<std::endl;

  // printf("Integrating polynoms :\n");
  // this->print();
  // other.print();
  // printf("over elem %d\n", iElm);

  int n = geoSpace->getNumQuadPoints();
  // printf("looping on n = %d\n", n);
  const std::vector<double> &J = geoSpace->getCncGeo()->getJacobians();
  std::vector<double> &w = geoSpace->getQuadratureWeights();
  // printf("Jsize = %d\n", J.size());
  // printf("wsize = %d\n", w.size());
  std::vector<double> x(3, 0.0);
  for(int i = 0; i < n; ++i) {
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, i, x);
    res += J[n * iElm + i] * w[i] * this->eval(x) * other.eval(x);
  }
  // printf("returning res = %f\n", res);
  return res;
}
