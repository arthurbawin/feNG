#ifndef _FEFUNCTION_
#define _FEFUNCTION_

#include <vector>
#include <functional>

// Convenience type for a scalar function f(t, [x y z])
typedef std::function<double(const double,
  const std::vector<double> &,
  const std::vector<double> &)> ScalarField;

// Convenience type for a scalar function f(u)
typedef std::function<double(const double,
  const std::vector<double> &)> ScalarSolField;

// Convenience type for a vector function F(t, [x y z])
typedef std::function<void(const double,
  const std::vector<double> &,
  const std::vector<double> &,
  std::vector<double> &)> VectorField;

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
  ~feFunction(){};

  // Evaluate the ScalarField at time t and position x
  // Can use 'eval' or operator() overloading
  double eval(const double t, const std::vector<double> &x){ return _fct1(t, x, _par); };
  double eval(const double u){ return _fct2(u, _par); };
  double operator()(const double t, const std::vector<double> &x){ return _fct1(t, x, _par); };
  double operator()(const double u){ return _fct2(u, _par); };
};

class feVectorFunction
{
protected:
  VectorField _fct;
  std::vector<double> _par;

public:
  feVectorFunction(VectorField fct, std::vector<double> par = std::vector<double>())
    : _fct(fct), _par(par){};
  ~feVectorFunction(){};

  // Evaluate the VectorField at time t and position x
  void eval(const double t, const std::vector<double> &x, std::vector<double> &res){ _fct(t, x, _par, res); };
  void operator()(const double t, const std::vector<double> &x, std::vector<double> &res){ _fct(t, x, _par, res); };
  std::vector<double> getParam() { return _par; }
};

#endif