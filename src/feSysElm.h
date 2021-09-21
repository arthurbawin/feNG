#ifndef _FESYSELM_
#define _FESYSELM_

#include "feFunction.h"
#include "feSpace.h"

class feSysElm {
protected:
  std::vector<int> _iVar;
  std::vector<int> _jVar;

  std::string _ID;

  bool matrixAnalyticalStatus = true;

public:
  feSysElm() {
    _iVar.resize(1);
    _jVar.resize(1);
  };
  virtual ~feSysElm() {}

  std::vector<int> &getIVar() { return _iVar; }
  std::vector<int> &getJVar() { return _jVar; }

  std::string getID() { return _ID; }

  bool getMatrixAnalyticalStatus() { return matrixAnalyticalStatus; };

  virtual void createElementarySystem(std::vector<feSpace *> &space) = 0;
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae) = 0;
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double *Be) = 0;
};

class feSysElm_0D_weakBC : public feSysElm {
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  int _idL;
  std::vector<double> _feU;
  std::vector<double> _feL;

public:
  feSysElm_0D_weakBC(double par, feFunction *fct) : _fct(fct), _par(par){};
  virtual ~feSysElm_0D_weakBC() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double *Be);
};

class feSysElm_0D_weakBC_edo1 : public feSysElm {
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  int _idL;
  int _idV;
  std::vector<double> _feU;
  std::vector<double> _feL;

public:
  feSysElm_0D_weakBC_edo1(double par, feFunction *fct) : _fct(fct), _par(par){};
  virtual ~feSysElm_0D_weakBC_edo1() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double *Be);
};

class feSysElm_0D_weakBC_edo2 : public feSysElm {
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  int _idL;
  int _idV;
  int _idW;
  std::vector<double> _feU;
  std::vector<double> _feL;

public:
  feSysElm_0D_weakBC_edo2(double par, feFunction *fct) : _fct(fct), _par(par){};
  virtual ~feSysElm_0D_weakBC_edo2() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double *Be);
};

class feSysElm_0D_Masse : public feSysElm {
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;

public:
  feSysElm_0D_Masse(double par, feFunction *fct) : _fct(fct), _par(par){};
  virtual ~feSysElm_0D_Masse() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double *Be);
};

class feSysElm_0D_Source : public feSysElm {
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;

public:
  feSysElm_0D_Source(double par, feFunction *fct) : _fct(fct), _par(par){};
  virtual ~feSysElm_0D_Source() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae){};
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double *Be);
};

class feSysElm_0D_Source_crossed : public feSysElm {
protected:
  feVectorFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;

public:
  feSysElm_0D_Source_crossed(double par, feVectorFunction *fct) : _fct(fct), _par(par){};
  virtual ~feSysElm_0D_Source_crossed() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae){};
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double *Be);
};

class feSysElm_1D_Source : public feSysElm {
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;

public:
  feSysElm_1D_Source(double par, feFunction *fct) : _fct(fct), _par(par){};
  virtual ~feSysElm_1D_Source() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae){};
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double *Be);
};

class feSysElm_1D_Diffusion : public feSysElm {
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;

public:
  feSysElm_1D_Diffusion(double par, feFunction *fct) : _fct(fct), _par(par){};
  virtual ~feSysElm_1D_Diffusion() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double *Be);
};

class feSysElm_1D_Masse : public feSysElm {
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;

public:
  feSysElm_1D_Masse(double par, feFunction *fct) : _fct(fct), _par(par){};
  virtual ~feSysElm_1D_Masse() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double *Be);
};

class feSysElm_2D_Source : public feSysElm {
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;
  std::vector<double> _feUdy;

public:
  feSysElm_2D_Source(double par, feFunction *fct) : _fct(fct), _par(par) { _ID = "source2D"; };
  virtual ~feSysElm_2D_Source() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae){};
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double *Be);
};

class feSysElm_2D_Diffusion : public feSysElm {
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;
  std::vector<double> _feUdy;

public:
  feSysElm_2D_Diffusion(double par, feFunction *fct) : _fct(fct), _par(par) { _ID = "diff2D"; };
  virtual ~feSysElm_2D_Diffusion() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double *Be);
};

class feSysElm_2D_Stokes : public feSysElm {
protected:
  feVectorFunction *_fct;
  std::vector<double> &_par; // Parametre
  int _idU;
  int _idV;
  int _idP;
  std::vector<double> _feU;
  std::vector<double> _feUdx;
  std::vector<double> _feUdy;
  std::vector<double> _feV;
  std::vector<double> _feVdx;
  std::vector<double> _feVdy;
  std::vector<double> _feP;

public:
  feSysElm_2D_Stokes(std::vector<double> &par, feVectorFunction *fct) : _fct(fct), _par(par) {
    _ID = "Stokes2D";
    matrixAnalyticalStatus = true;
  };
  virtual ~feSysElm_2D_Stokes() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double *Be);
};

class feSysElm_2D_NavierStokes : public feSysElm {
protected:
  feVectorFunction *_fct;
  std::vector<double> &_par; // Parametre
  int _idU;
  int _idV;
  int _idP;
  std::vector<double> _feU;
  std::vector<double> _feUdx;
  std::vector<double> _feUdy;
  std::vector<double> _feV;
  std::vector<double> _feVdx;
  std::vector<double> _feVdy;
  std::vector<double> _feP;

public:
  feSysElm_2D_NavierStokes(std::vector<double> &par, feVectorFunction *fct) : _fct(fct), _par(par) {
    _ID = "NavierStokes2D";
    matrixAnalyticalStatus = true;
  };
  virtual ~feSysElm_2D_NavierStokes() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double *Be);
};

#endif