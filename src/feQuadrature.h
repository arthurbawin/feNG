#ifndef _FEQUADRATURE_
#define _FEQUADRATURE_

#include <vector>
#include <string>

#include "feCncGeo.h"

// From Gmsh's GaussIntegration.h
struct IntPt {
  double pt[3];
  double weight;
};

// Hardcoded quadrature points and weights for triangles
int getNGQTPts(int order);
IntPt *getGQTPts(int order);

// Hardcoded quadrature points and weights for tetrahedra
int getNGQTetPts(int order);
IntPt *getGQTetPts(int order);

class feQuadrature
{
protected:
  int _dim;
  int _degQuad;
  int _nQuad1D;
  std::vector<double> _w1D;
  std::vector<double> _x1D;

  int _nQuad;
  std::vector<double> _w;
  std::vector<double> _xr;
  std::vector<double> _yr;
  std::vector<double> _zr;

  double result[2];

public:
  feQuadrature(int degree, geometryType geometry);
  ~feQuadrature(){};

private:
  void getGaussPoints(int nQuadLocal);
  void computePolynomialValueAndDerivative(double x, double *result);
  void getGaussPointsTri(int degree, int method);
  void getGaussPointsSquare();
  void getGaussPointsCube();
  void getGaussPointsTetra(int degree);

public:
  std::vector<double> getXPoints() { return _xr; }
  std::vector<double> getYPoints() { return _yr; }
  std::vector<double> getZPoints() { return _zr; }
  std::vector<double> getWeights() { return _w; }
  int getNumQuadPoints() { return _nQuad; }
};

#endif