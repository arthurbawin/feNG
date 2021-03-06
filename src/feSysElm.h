#ifndef _FESYSELM_
#define _FESYSELM_

#include "feFunction.h"
#include "feSpace.h"

class feSysElm
{
protected:
  std::vector<int> _iVar;
  std::vector<int> _jVar;

  std::string _ID;

  bool matrixAnalyticalStatus = true;

  // Does this form has an elementary matrix to assemble ?
  bool _hasMatrix = false;

public:
  feSysElm(bool hasMatrix) : _hasMatrix(hasMatrix)
  {
    _iVar.resize(1);
    _jVar.resize(1);
  };
  virtual ~feSysElm() {}

  std::vector<int> &getIVar() { return _iVar; }
  std::vector<int> &getJVar() { return _jVar; }

  std::string getID() { return _ID; }

  bool getMatrixAnalyticalStatus() { return matrixAnalyticalStatus; };

  bool hasMatrix() { return _hasMatrix; }

  virtual void createElementarySystem(std::vector<feSpace *> &space) = 0;
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae) = 0;
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be) = 0;
};

class feSysElm_0D_StiffSpring : public feSysElm
{
protected:
  feFunction *_fct;
  std::vector<double> _par; // Parametre
  int _idX; // pos
  int _idV; // vit
  std::vector<double> _feX;
  std::vector<double> _feV;

public:
  feSysElm_0D_StiffSpring(std::vector<double> par, feFunction *fct)
    : feSysElm(true), _fct(fct), _par(par){};
  virtual ~feSysElm_0D_StiffSpring() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be);
};

class feSysElm_0D_Stiff2 : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idX;
  int _idY;
  int _idZ;
  std::vector<double> _feX;
  std::vector<double> _feY;
  std::vector<double> _feZ;

public:
  feSysElm_0D_Stiff2(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par){};
  virtual ~feSysElm_0D_Stiff2() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be);
};

class feSysElm_0D_Stiff3 : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idX;
  int _idY;
  int _idZ;
  std::vector<double> _feX;
  std::vector<double> _feY;
  std::vector<double> _feZ;

public:
  feSysElm_0D_Stiff3(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par){};
  virtual ~feSysElm_0D_Stiff3() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be);
};

class feSysElm_0D_weakBC : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  int _idL;
  std::vector<double> _feU;
  std::vector<double> _feL;

public:
  feSysElm_0D_weakBC(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par){};
  virtual ~feSysElm_0D_weakBC() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be);
};

class feSysElm_0D_weakBC_edo1 : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  int _idL;
  int _idV;
  std::vector<double> _feU;
  std::vector<double> _feL;

public:
  feSysElm_0D_weakBC_edo1(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par){};
  virtual ~feSysElm_0D_weakBC_edo1() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be);
};

class feSysElm_0D_weakBC_edo1_V2 : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  int _idL;
  int _idV;
  std::vector<double> _feU;
  std::vector<double> _feL;

public:
  feSysElm_0D_weakBC_edo1_V2(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par){};
  virtual ~feSysElm_0D_weakBC_edo1_V2() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be);
};

class feSysElm_0D_weakBC_edo2 : public feSysElm
{
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
  feSysElm_0D_weakBC_edo2(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par){};
  virtual ~feSysElm_0D_weakBC_edo2() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be);
};

class feSysElm_0D_Masse : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;

public:
  feSysElm_0D_Masse(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par){};
  virtual ~feSysElm_0D_Masse() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be);
};

class feSysElm_0D_Source : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;

public:
  feSysElm_0D_Source(double par, feFunction *fct) : feSysElm(false), _fct(fct), _par(par){};
  virtual ~feSysElm_0D_Source() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae){};
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be);
};

class feSysElm_0D_Source_crossed : public feSysElm
{
protected:
  feVectorFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;

public:
  feSysElm_0D_Source_crossed(double par, feVectorFunction *fct)
    : feSysElm(false), _fct(fct), _par(par){};
  virtual ~feSysElm_0D_Source_crossed() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae){};
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be);
};

// class feSysElm_1D_weakBC_Vec : public feSysElm {
// protected:
//   feVectorFunction *_fct;
//   std::vector<double> &_par; // Parametre
//   int _idU;
//   int _idV;
//   int _idL;
//   int _idBu;
//   int _idBv;
//   std::vector<double> _feU;
//   std::vector<double> _feV;
//   std::vector<double> _feL;
//   std::vector<double> _feBu;
//   std::vector<double> _feBv;

// public:
//   feSysElm_1D_weakBC_Vec(std::vector<double> &par, feVectorFunction *fct) : _fct(fct),
//   _par(par){}; virtual ~feSysElm_1D_weakBC_Vec() {}

//   virtual void createElementarySystem(std::vector<feSpace *> &space);
//   virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
//                          feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
//                          double **Ae);
//   virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
//                          feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
//                          double dt, double *Be);
// };

class feSysElm_1D_weakBC_edo1 : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  int _idL;
  int _idV;
  std::vector<double> _feU;
  std::vector<double> _feV;
  std::vector<double> _feL;

public:
  feSysElm_1D_weakBC_edo1(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par){};
  virtual ~feSysElm_1D_weakBC_edo1() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be);
};

class feSysElm_1D_Source : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;

public:
  feSysElm_1D_Source(double par, feFunction *fct) : feSysElm(false), _fct(fct), _par(par){};
  virtual ~feSysElm_1D_Source() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae){};
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be);
};

class feSysElm_1D_Diffusion : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;

public:
  feSysElm_1D_Diffusion(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par){};
  virtual ~feSysElm_1D_Diffusion() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be);
};

class feSysElm_1D_Masse : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;

public:
  feSysElm_1D_Masse(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par){};
  virtual ~feSysElm_1D_Masse() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be);
};

// class feSysElm_1D_NeumannBC : public feSysElm {
// protected:
//   feFunction *_fct;
//   double _par; // Parametres
//   int _idU;
//   std::vector<double> _feU;
//   std::vector<double> _feUdx;

// public:
//   feSysElm_1D_Masse(double par, feFunction *fct) : _fct(fct), _par(par){};
//   virtual ~feSysElm_1D_Masse() {}

//   virtual void createElementarySystem(std::vector<feSpace *> &space);
//   virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
//                          feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
//                          double **Ae);
//   virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
//                          feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
//                          double dt, double *Be);
// };

class feSysElm_2D_Masse : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;
  std::vector<double> _feUdy;

public:
  feSysElm_2D_Masse(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par){};
  virtual ~feSysElm_2D_Masse() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be);
};

class feSysElm_2D_Source : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;
  std::vector<double> _feUdy;

public:
  feSysElm_2D_Source(double par, feFunction *fct) : feSysElm(false), _fct(fct), _par(par)
  {
    _ID = "source2D";
  };
  virtual ~feSysElm_2D_Source() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae){};
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be);
};

class feSysElm_2D_Diffusion : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;
  std::vector<double> _feUdy;

public:
  feSysElm_2D_Diffusion(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par)
  {
    _ID = "diff2D";
    matrixAnalyticalStatus = false;
  };
  virtual ~feSysElm_2D_Diffusion() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be);
};

class feSysElm_2D_Advection : public feSysElm
{
protected:
  feVectorFunction *_fct; // The imposed velocity field
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;
  std::vector<double> _feUdy;

public:
  feSysElm_2D_Advection(double par, feVectorFunction *fct) : feSysElm(true), _fct(fct), _par(par)
  {
    _ID = "adv2D";
    matrixAnalyticalStatus = false;
  };
  virtual ~feSysElm_2D_Advection() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be);
};

class feSysElm_2D_Stokes : public feSysElm
{
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
  feSysElm_2D_Stokes(std::vector<double> &par, feVectorFunction *fct)
    : feSysElm(true), _fct(fct), _par(par)
  {
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
                         double dt, double *Be);
};

class feSysElm_2D_NavierStokes : public feSysElm
{
protected:
  feVectorFunction *_fct;
  feFunction *_viscosityFct;
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
  feSysElm_2D_NavierStokes(std::vector<double> &par, feVectorFunction *fct,
                           feFunction *viscosityFct = nullptr)
    : feSysElm(true), _fct(fct), _par(par), _viscosityFct(viscosityFct)
  {
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
                         double dt, double *Be);
};

// Outflow boundary condition
class feSysElm_1D_DirectionalDoNothing : public feSysElm
{
protected:
  feVectorFunction *_fct;
  feFunction *_viscosityFct;
  std::vector<double> &_par; // Parametre
  int _idU;
  int _idV;
  std::vector<double> _feU;
  std::vector<double> _feUdx;
  std::vector<double> _feUdy;
  std::vector<double> _feV;
  std::vector<double> _feVdx;
  std::vector<double> _feVdy;

public:
  feSysElm_1D_DirectionalDoNothing(std::vector<double> &par, feVectorFunction *fct,
                                   feFunction *viscosityFct = nullptr)
    : feSysElm(true), _fct(fct), _par(par), _viscosityFct(viscosityFct)
  {
    _ID = "DirectionalDoNothing1D";
    matrixAnalyticalStatus = false;
  };
  virtual ~feSysElm_1D_DirectionalDoNothing() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  // Matrix not implemented : finite differences only
  virtual void computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae);
  virtual void computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be);
};

#endif
