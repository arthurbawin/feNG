#ifndef _FEFUNCTION_
#define _FEFUNCTION_

#include <vector>
#include <functional>

class feFunction {
protected:
  std::function<double(const double, const std::vector<double> &, const std::vector<double> &)>
    _fct;
  std::vector<double> _par;

public:
  feFunction(
    std::function<double(const double, const std::vector<double> &, const std::vector<double> &)>
      fct,
    std::vector<double> par = std::vector<double>())
    : _fct(fct), _par(par){};
  ~feFunction() {}

  double eval(const double t, const std::vector<double> &x) { return _fct(t, x, _par); }
};

class feVectorFunction {
protected:
  std::function<void(const double, const std::vector<double> &, const std::vector<double> &,
                     std::vector<double> &)>
    _fct;
  std::vector<double> _par;

public:
  feVectorFunction(std::function<void(const double, const std::vector<double> &,
                                      const std::vector<double> &, std::vector<double> &)>
                     fct,
                   std::vector<double> par = std::vector<double>())
    : _fct(fct), _par(par){};
  ~feVectorFunction() {}

  void eval(const double t, const std::vector<double> &x, std::vector<double> &res) {
    return _fct(t, x, _par, res);
  }
};

#endif