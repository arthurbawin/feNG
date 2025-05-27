#ifndef _FEFUNCTION_
#define _FEFUNCTION_

#include <functional>
#include <vector>

#include "feMatrixInterface.h"

//
// Wrapper for the possible arguments of a function callback
// Other members can be added if necessary (e.g. solution derivatives, etc.)
//
class feFunctionArguments
{
public:
  // The current time
  double t = 0.;
  // Position
  std::vector<double> pos = {0., 0., 0.};

  // FIXME: We should be able to give any number of scalar/vector fields
  // Scalar-valued solution or variable
  double u = 0.;

  //
  std::vector<double> userScalarFields;

  // Vector-valued solution or variable
  std::vector<double> v = {0., 0., 0.};

  feFunctionArguments(const double time = 0.) { t = time; }
  feFunctionArguments(const double time, const std::vector<double> &position)
  {
    t      = time;
    pos[0] = position[0];
    pos[1] = position[1];
    pos[2] = position[2];
  };
};

// Convenience type for a scalar-valued function f(arguments, parameters).
// Most functions are typically of the form f(t,x,y,z),
// but this also allows functions depending on a problem variable, e.g. f(u)
typedef std::function<double(const feFunctionArguments &args,
                             const std::vector<double> &param)>
  ScalarField;

// Convenience type for a vector-valued function f(arguments, parameters)
typedef std::function<void(const feFunctionArguments &args,
                           const std::vector<double> &param,
                           std::vector<double>       &res)>
  VectorField;

// Convenience type for a SPD 2x2 matrix-valued function f(arguments,
// parameters)
typedef std::function<void(const feFunctionArguments &args,
                           const std::vector<double> &param,
                           MetricTensor              &res)>
  MetricTensorField;

static double constantCallback(const feFunctionArguments & /*args*/,
                               const std::vector<double> &par)
{
  return par[0];
}

static void constantVectorCallback(const feFunctionArguments & /*args*/,
                                   const std::vector<double> &par,
                                   std::vector<double>       &res)
{
  res[0] = par[0];
  res[1] = par[1];
}

//
// Wrappers for scalar fields
//
class feFunction
{
protected:
  ScalarField         _fct;
  std::vector<double> _param;

public:
  feFunction(ScalarField fct, std::vector<double> param = std::vector<double>())
    : _fct(fct)
    , _param(param){};

  void resetParameters(const std::vector<double> &newParameters)
  {
    for (size_t i = 0; i < _param.size(); ++i)
      _param[i] = newParameters[i];
  };

  // Evaluate the ScalarField at time t and position x
  // Can use 'eval' or operator() overloading
  double eval(const feFunctionArguments &args) const
  {
    return _fct(args, _param);
  };
  double operator()(const feFunctionArguments &args) const
  {
    return _fct(args, _param);
  };

  // Prevent implicit conversion when calling eval or operator()
  template <class T>
  void eval(const T &) const = delete;
  template <class T>
  void operator()(const T &) const = delete;
};

class feConstantFunction : public feFunction
{
public:
  feConstantFunction(double C)
    : feFunction(constantCallback, {C}){};
};

//
// Wrappers for vector fields
//
class feVectorFunction
{
protected:
  VectorField         _fct;
  std::vector<double> _param;

public:
  feVectorFunction(VectorField         fct,
                   std::vector<double> param = std::vector<double>())
    : _fct(fct)
    , _param(param){};

  // Evaluate the VectorField at time t and position x or at solution u
  void eval(const feFunctionArguments &args, std::vector<double> &res) const
  {
    _fct(args, _param, res);
  };
  void operator()(const feFunctionArguments &args,
                  std::vector<double>       &res) const
  {
    _fct(args, _param, res);
  };

  // Prevent implicit conversion when calling eval or operator()
  template <class T>
  void eval(const T &, std::vector<double> &res) const = delete;
  template <class T>
  void operator()(const T &, std::vector<double> &res) const = delete;

  std::vector<double> getParam() { return _param; }
};

class feConstantVectorFunction : public feVectorFunction
{
public:
  feConstantVectorFunction(std::vector<double> C)
    : feVectorFunction(constantVectorCallback, C){};
};

//
// Wrappers for fields of metric tensors (2x2 SPD matrices in 2D)
//
class feMetricTensorFunction
{
protected:
  MetricTensorField   _fct;
  std::vector<double> _par;

public:
  feMetricTensorFunction(MetricTensorField   fct,
                         std::vector<double> par = std::vector<double>())
    : _fct(fct)
    , _par(par){};
  ~feMetricTensorFunction(){};

  // Evaluate the MetricTensorField at time t and position x
  void eval(const feFunctionArguments &args, MetricTensor &res) const
  {
    _fct(args, _par, res);
  };
  void operator()(const feFunctionArguments &args, MetricTensor &res) const
  {
    _fct(args, _par, res);
  };
  std::vector<double> getParam() { return _par; }
};

#endif