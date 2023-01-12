#ifndef _FEFUNCTION_
#define _FEFUNCTION_

#include <vector>
#include <functional>

typedef std::function<double(const double,
  const std::vector<double> &,
  const std::vector<double> &)> ScalarField;

typedef std::function<void(const double,
  const std::vector<double> &,
  const std::vector<double> &,
  std::vector<double> &)> VectorField;

// Wrappers for scalar and vector fields
class feFunction
{
protected:
  ScalarField _fct;
  std::vector<double> _par;

public:
  feFunction(ScalarField fct, std::vector<double> par = std::vector<double>())
    : _fct(fct), _par(par){};
  ~feFunction() {}

  // Evaluate the ScalarField at time t and position x
  double eval(const double t, const std::vector<double> &x) { return _fct(t, x, _par); }
};

class feVectorFunction
{
protected:
  VectorField _fct;
  std::vector<double> _par;

public:
  feVectorFunction(VectorField fct, std::vector<double> par = std::vector<double>())
    : _fct(fct), _par(par){};
  ~feVectorFunction() {}

  // Evaluate the VectorField at time t and position x
  void eval(const double t, const std::vector<double> &x, std::vector<double> &res)
  {
    return _fct(t, x, _par, res);
  }
  std::vector<double> getParam() { return _par; }
};

#endif