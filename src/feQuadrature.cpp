#include "feQuadrature.h"

#include <cmath>

const double EPSILON = 1e-16;

feQuadrature::feQuadrature(int degree, geometryType geometry) : _degQuad(degree)
{
  _nQuad1D = (_degQuad + 2) / 2;
  _w1D.resize(_nQuad1D);
  _x1D.resize(_nQuad1D);

  if(geometry == geometryType::POINT) _dim = 0;
  if(geometry == geometryType::LINE) _dim = 1;
  if(geometry == geometryType::TRI || geometry == geometryType::QUAD) _dim = 2;
  if(geometry == geometryType::TET || geometry == geometryType::HEX) _dim = 3;

  if(geometry == geometryType::POINT) {
    _w.resize(1);
    _w[0] = 1.;
    _xr.resize(1);
    _xr[0] = 0.;
    _yr.resize(1);
    _zr.resize(1);
  } else if(geometry == geometryType::LINE) {
    computeWeightAndRoot(_nQuad1D);
    _w = _w1D;
    _xr = _x1D;
    _yr.resize(_nQuad1D);
    _zr.resize(_nQuad1D);
  } else if(geometry == geometryType::TRI) {
    computeWeightAndRootTri(degree, 1);
  } else if(geometry == geometryType::QUAD) {
    computeWeightAndRootSquare();
  } else if(geometry == geometryType::HEX) {
    computeWeightAndRootCube();
  } else if(geometry == geometryType::TET) {
    computeWeightAndRootTetra();
  } else {
    feErrorMsg(FE_STATUS_ERROR, "Could not create a quadrature rule for geometry \"%s\".",
               toString(geometry).data());
    exit(-1);
  }

  _nQuad = _w.size();
}

void feQuadrature::computeWeightAndRoot(int nQuadLocal)
{
  _w1D.resize(nQuadLocal);
  _x1D.resize(nQuadLocal);

  if(nQuadLocal == 1) {
    _x1D[0] = 0.0;
    _w1D[0] = 2.0;
  } else {
    _nQuad1D = nQuadLocal; // On met a jour le nombre de point que l'on veut pour la règle
    for(int i = 0; i < _nQuad1D; i++) {
      double root = cos(M_PI * ((i + 1) - 0.25) / (_nQuad1D + 0.5));
      computePolynomialValueAndDerivative(root, result);
      double newtonRaphsonRatio;
      do {
        newtonRaphsonRatio = result[0] / result[1];
        root -= newtonRaphsonRatio;
        computePolynomialValueAndDerivative(root, result);
      } while(fabs(newtonRaphsonRatio) > EPSILON);
      _x1D[_nQuad1D - 1 - i] = root;
      _w1D[_nQuad1D - 1 - i] = 2.0 / ((1 - root * root) * result[1] * result[1]);
    }
  }
}

void feQuadrature::computePolynomialValueAndDerivative(double x, double *result)
{
  result[0] = x;
  result[1] = 0.0;
  double value_minus_1 = 1;
  const double f = 1 / (x * x - 1);
  for(int i = 2; i <= _nQuad1D; i++) {
    const double value = ((2 * i - 1) * x * result[0] - (i - 1) * value_minus_1) / i;
    result[1] = i * f * (x * value - result[0]);
    value_minus_1 = result[0];
    result[0] = value;
  }
}

void feQuadrature::computeWeightAndRootTri(int degree, int method)
{
  if(degree <= 20) {
    // Get hardcoded quadrature points from Gmsh
    // Rule with negative weights only for order 20
    IntPt *rule = getGQTPts(degree);
    int nQuad = getNGQTPts(degree);

    _w.resize(nQuad);
    _xr.resize(nQuad);
    _yr.resize(nQuad);
    _zr.resize(nQuad);

    for(int i = 0; i < nQuad; ++i) {
      _xr[i] = rule[i].pt[0];
      _yr[i] = rule[i].pt[1];
      _zr[i] = rule[i].pt[2];
      _w[i] = rule[i].weight;
    }

  } else {
    if(method == 1) {
      // Tensor rule with n^2 points
      _w.resize(pow(_nQuad1D, _dim));
      _xr.resize(pow(_nQuad1D, _dim));
      _yr.resize(pow(_nQuad1D, _dim));
      _zr.resize(pow(_nQuad1D, _dim));

      _nQuad1D = (_degQuad + 2) / 2;
      computeWeightAndRoot(_nQuad1D);
      int l = 0;
      for(int i = 0; i < _nQuad1D; ++i) {
        for(int j = 0; j < _nQuad1D; ++j) {
          _yr[l] = (1 + _x1D[i]) / 2.0;
          _xr[l] = (1 - _x1D[i]) * (1 + _x1D[j]) / 4.0;
          _w[l] = _w1D[j] * _w1D[i] * (1 - _x1D[i]) / 8.0;
          ++l;
        }
      }
    } else {
      // Lighter version with (n(n+1)/2 - 1 points) - Source?
      _w.resize(_degQuad * (_degQuad + 1) / 2 - 1);
      _xr.resize(_degQuad * (_degQuad + 1) / 2 - 1);
      _yr.resize(_degQuad * (_degQuad + 1) / 2 - 1);
      _zr.resize(_degQuad * (_degQuad + 1) / 2 - 1);

      computeWeightAndRoot(_degQuad - 1);
      std::vector<double> _xi = _x1D;
      std::vector<double> _wi = _w1D;
      int l = 0;
      for(int i = 0; i < _degQuad - 1; ++i) {
        computeWeightAndRoot(_degQuad - i);
        std::vector<double> _eta = _x1D;
        std::vector<double> _wj = _w1D;
        for(int j = 0; j < _degQuad - i; ++j) {
          _yr[l] = (1 + _xi[i]) / 2.0;
          _xr[l] = (1 - _xi[i]) * (1 + _eta[j]) / 4.0;
          _w[l] = _wj[j] * _wi[i] * (1 - _xi[i]) / 8.0;
          ++l;
        }
      }
    }
  }
}

void feQuadrature::computeWeightAndRootTetra()
{
  _w.resize(pow(_nQuad1D, _dim));
  _xr.resize(pow(_nQuad1D, _dim));
  _yr.resize(pow(_nQuad1D, _dim));
  _zr.resize(pow(_nQuad1D, _dim));
  computeWeightAndRoot(_nQuad1D);
  int l = 0;
  for(int i = 0; i < _nQuad1D; ++i) {
    for(int j = 0; j < _nQuad1D; ++j) {
      for(int k = 0; k < _nQuad1D; ++k) {
        _xr[l] = (1 - _x1D[i]) / 2;
        _yr[l] = (1 + _x1D[i]) * (1 - _x1D[j]) / 4;
        _zr[l] = (1 + _x1D[i]) * (1 + _x1D[j]) * (1 - _x1D[k]) / 8;
        _w[l] = _w1D[j] * _w1D[i] * _w1D[k] * (1 + _x1D[i]) * (1 + _x1D[i]) * (1 + _x1D[j]) / 64;
        ++l;
      }
    }
  }
}

void feQuadrature::computeWeightAndRootSquare()
{
  _w.resize(pow(_nQuad1D, _dim));
  _xr.resize(pow(_nQuad1D, _dim));
  _yr.resize(pow(_nQuad1D, _dim));
  _zr.resize(pow(_nQuad1D, _dim));
  computeWeightAndRoot(_nQuad1D);
  int l = 0;
  for(int i = 0; i < _nQuad1D; ++i) {
    for(int j = 0; j < _nQuad1D; ++j) {
      _xr[l] = _x1D[i];
      _yr[l] = _x1D[j];
      _w[l] = _w1D[j] * _w1D[i];
      ++l;
    }
  }
}

void feQuadrature::computeWeightAndRootCube()
{
  _w.resize(pow(_nQuad1D, _dim));
  _xr.resize(pow(_nQuad1D, _dim));
  _yr.resize(pow(_nQuad1D, _dim));
  _zr.resize(pow(_nQuad1D, _dim));
  computeWeightAndRoot(_nQuad1D);
  int l = 0;
  for(int i = 0; i < _nQuad1D; ++i) {
    for(int j = 0; j < _nQuad1D; ++j) {
      for(int k = 0; k < _nQuad1D; ++k) {
        _xr[l] = _x1D[i];
        _yr[l] = _x1D[j];
        _zr[l] = _x1D[k];
        _w[l] = _w1D[j] * _w1D[i] * _w1D[k];
        ++l;
      }
    }
  }
}