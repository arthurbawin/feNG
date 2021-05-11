#ifndef _FESYSELM_
#define _FESYSELM_

#include "feSpace.h"

class feSysElm{

protected:
  std::function<double(const double, const std::vector<double> &)> _fct;

  std::vector<int> _iVar;
  std::vector<int> _jVar;

public:
	feSysElm(std::function<double(const double, const std::vector<double> &)> fct)
    : _fct(fct)
	{
    _iVar.resize(1);
    _jVar.resize(1);
	};
	virtual ~feSysElm() {}

  std::vector<int> &getIVar(){ return _iVar; }
  std::vector<int> &getJVar(){ return _jVar; }

  virtual void createElementarySystem(std::vector<feSpace*> &space);
  virtual void computeAe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn, std::vector<double> &Ae) {};
  virtual void computeBe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn, std::vector<double> &Be) {};
};

class feSysElm_1D_Source : public feSysElm{
protected:
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;
public:
  feSysElm_1D_Source(double par, std::function<double(const double, const std::vector<double> &)> fct) : _par(par), feSysElm(fct){};
  virtual ~feSysElm_1D_Source() {}

  virtual void createElementarySystem(std::vector<feSpace*> &space);
  virtual void computeAe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn, std::vector<double> &Ae) {};
  virtual void computeBe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn, std::vector<double> &Be);
};

class feSysElm_1D_Diffusion : public feSysElm{
protected:
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;
public:
  feSysElm_1D_Diffusion(double par, std::function<double(const double, const std::vector<double> &)> fct) : _par(par), feSysElm(fct){};
  virtual ~feSysElm_1D_Diffusion() {}
  
  virtual void createElementarySystem(std::vector<feSpace*> &space);
  virtual void computeAe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn, std::vector<double> &Ae);
  virtual void computeBe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn, std::vector<double> &Be);
};

#endif