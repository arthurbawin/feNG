#ifndef _FESYSELM_
#define _FESYSELM_

#include "feFunction.h"
#include "feSpace.h"

class feSysElm{

protected:
  feFunction *_fct;

  std::vector<int> _iVar;
  std::vector<int> _jVar;

public:
	feSysElm(feFunction *fct)
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
  feSysElm_1D_Source(double par, feFunction *fct) : feSysElm(fct), _par(par){};
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
  feSysElm_1D_Diffusion(double par, feFunction *fct) : feSysElm(fct), _par(par){};
  virtual ~feSysElm_1D_Diffusion() {}
  
  virtual void createElementarySystem(std::vector<feSpace*> &space);
  virtual void computeAe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn, std::vector<double> &Ae);
  virtual void computeBe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn, std::vector<double> &Be);
};

class feSysElm_1D_Masse : public feSysElm{
protected:
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;
public:
  feSysElm_1D_Masse(double par, feFunction *fct) : feSysElm(fct), _par(par){};
  virtual ~feSysElm_1D_Masse() {}
  
  virtual void createElementarySystem(std::vector<feSpace*> &space);
  virtual void computeAe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn, std::vector<double> &Ae);
  virtual void computeBe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn, std::vector<double> &Be);
};

#endif