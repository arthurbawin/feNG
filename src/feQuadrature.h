#ifndef _FEQUADRATURE_
#define _FEQUADRATURE_

#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <utility>

class feQuadrature{

protected:
  int _degQuad;
  int _dimQuad;
  int _nQuad1D;
  int _nQuadLocal;
  std::vector<double> _w1D;
  std::vector<double> _x1D;
  std::vector<double> _w;
  std::vector<double> _xr;
  std::vector<double> _yr;
  std::vector<double> _zr;

  const double EPSILON = 1e-16;
  double result[2]; 

public:
  feQuadrature(int degQuad, int dimQuad, std::string formGeo);
  ~feQuadrature() {}

private : 
  void calculateWeightAndRoot(int nQuadLocal); 
  void calculatePolynomialValueAndDerivative(double x, double *result); 
  void calculateWeightAndRootTri(int method);
  void calculateWeightAndRootSquare();
  void calculateWeightAndRootCube();
  void calculateWeightAndRootTetra();

public :
  std::vector<double> getXPoints(){ return _xr; }
  std::vector<double> getYPoints(){ return _yr; }
  std::vector<double> getZPoints(){ return _zr; }
  std::vector<double> getWeights(){ return _w; }
  int getNQuad() { return _w.size(); }
  int getDim(){ return _dimQuad; }
};

#endif