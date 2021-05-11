#ifndef _FEFUNCTION_
#define _FEFUNCTION_

#include <vector>
#include <functional>

class feFunction{

protected:
	std::function<double(const double, const std::vector<double> &, const std::vector<double>)> _fct;
	std::vector<double> _par;

public:
	feFunction(std::function<double(const double, const std::vector<double> &, const std::vector<double>)> fct,
		std::vector<double> par = std::vector<double>()) : _fct(fct), _par(par)
	{
	};
	~feFunction() {}

  double eval(const double t, const std::vector<double> &x){ return _fct(t,x,_par); }
};

#endif