#ifndef _FEFUNCTION_
#define _FEFUNCTION_

#include <vector>
#include <functional>

#include "feMatrixInterface.h"

// Convenience type for a scalar function f(t, [x y z])
typedef std::function<double(const double, const std::vector<double> &,
                             const std::vector<double> &)> ScalarField;

// Convenience type for a scalar function f(u)
typedef std::function<double(const double, const std::vector<double> &)> ScalarSolField;

// Convenience type for a vector function F(t, [x y z])
typedef std::function<void(const double, const std::vector<double> &, const std::vector<double> &,
                           std::vector<double> &)> VectorField;

// Convenience type for a vector function F(u)
// FIXME: Not ideal because it is identical to ScalarField!
// This allows the creation of feVectorFunctions with ScalarField callbacks
typedef std::function<void(const double, const std::vector<double> &, std::vector<double> &)> VectorSolField;

// Convenience type for a function M(t, [x y z]) that returns a 2x2 matrix
typedef std::function<void(const double, const std::vector<double> &, const std::vector<double> &,
                           MetricTensor &)>
  MetricTensorField;

static double constantCallback(const double t, const std::vector<double> &pos,
                               const std::vector<double> &par)
{
  UNUSED(t, pos);
  return par[0];
}

static double constantSolCallback(const double u, const std::vector<double> &par)
{
  UNUSED(u);
  return par[0];
}

static void constantVectorCallback(const double t, const std::vector<double> &pos,
                                   const std::vector<double> &par, std::vector<double> &res)
{
  UNUSED(t, pos);
  res[0] = par[0];
  res[1] = par[1];
}

static void constantVectorSolCallback(const double u, const std::vector<double> &par, std::vector<double> &res)
{
  UNUSED(u);
  res[0] = par[0];
  res[1] = par[1];
}

// Wrappers for scalar and vector fields
class feFunction
{
protected:
  ScalarField _fct1;
  ScalarSolField _fct2;
  std::vector<double> _par;

public:
  feFunction(ScalarField fct, std::vector<double> par = std::vector<double>())
    : _fct1(fct), _par(par){};
  feFunction(ScalarSolField fct, std::vector<double> par = std::vector<double>())
    : _fct2(fct), _par(par){};
  feFunction(ScalarField fct1, ScalarSolField fct2, std::vector<double> par = std::vector<double>())
    : _fct1(fct1), _fct2(fct2), _par(par){};
  ~feFunction(){};

  void resetParameters(const std::vector<double> &newParameters)
  {
    for(size_t i = 0; i < _par.size(); ++i)
      _par[i] = newParameters[i];
  };

  // Evaluate the ScalarField at time t and position x
  // Can use 'eval' or operator() overloading
  double eval(const double t, const std::vector<double> &x) { return _fct1(t, x, _par); };
  double eval(const double u) { return _fct2(u, _par); };
  double operator()(const double t, const std::vector<double> &x) { return _fct1(t, x, _par); };
  double operator()(const double u) { return _fct2(u, _par); };
};

class feConstantFunction : public feFunction
{
public:
  feConstantFunction(double C) : feFunction(constantCallback, constantSolCallback, {C}){};
  ~feConstantFunction(){};
};

class feVectorFunction
{
protected:
  VectorField _fct1;
  VectorSolField _fct2;
  std::vector<double> _par;

public:
  feVectorFunction(VectorField fct, std::vector<double> par = std::vector<double>())
    : _fct1(fct), _par(par){};
  feVectorFunction(VectorSolField fct, std::vector<double> par = std::vector<double>())
    : _fct2(fct), _par(par){};
  feVectorFunction(VectorField fct1, VectorSolField fct2, std::vector<double> par = std::vector<double>())
    : _fct1(fct1), _fct2(fct2), _par(par){};
  ~feVectorFunction(){};

  // Evaluate the VectorField at time t and position x or at solution u
  void eval(const double t, const std::vector<double> &x, std::vector<double> &res)
  {
    _fct1(t, x, _par, res);
  };
  void eval(const double u, std::vector<double> &res)
  {
    _fct2(u, _par, res);
  };
  void operator()(const double t, const std::vector<double> &x, std::vector<double> &res)
  {
    _fct1(t, x, _par, res);
  };
  void operator()(const double u, std::vector<double> &res)
  {
    _fct2(u, _par, res);
  };
  std::vector<double> getParam() { return _par; }
};

class feConstantVectorFunction : public feVectorFunction
{
public:
  feConstantVectorFunction(std::vector<double> C)
   : feVectorFunction(constantVectorCallback, constantVectorSolCallback, C){};
  ~feConstantVectorFunction(){};
};

class feMetricTensorFunction
{
protected:
  MetricTensorField _fct;
  std::vector<double> _par;

public:
  feMetricTensorFunction(MetricTensorField fct, std::vector<double> par = std::vector<double>())
    : _fct(fct), _par(par){};
  ~feMetricTensorFunction(){};

  // Evaluate the MetricTensorField at time t and position x
  void eval(const double t, const std::vector<double> &x, MetricTensor &res)
  {
    _fct(t, x, _par, res);
  };
  void operator()(const double t, const std::vector<double> &x, MetricTensor &res)
  {
    _fct(t, x, _par, res);
  };
  std::vector<double> getParam() { return _par; }
};

#endif