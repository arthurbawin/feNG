#ifndef _FESYSELM_
#define _FESYSELM_

#include "feFunction.h"
#include "feSpace.h"

typedef enum {
  STIFFSPRING_0D,
  STIFF2_0D,
  STIFF3_0D,
  WEAKBC_0D,
  WEAKBC_EDO1_0D,
  WEAKBC_EDO1_V2_0D,
  WEAKBC_EDO2_0D,
  MASS_0D,
  SOURCE_0D,
  SOURCE_CROSSED_0D,
  WEAKBC_EDO1_1D,
  SOURCE_1D,
  DIFFUSION_1D,
  MASS_1D,
  NEUMANN_1D,
  ADVECTION_1D,
  DG_ADVECTION_1D,
  SUPG_STABILIZATION_1D,
  MASS_2D,
  SOURCE_2D,
  DIFFUSION_2D,
  ADVECTION_2D,
  STOKES_2D,
  NAVIERSTOKES_2D,
  DIRECTIONALDONOTHING_1D
} elementSystemType;

inline const std::string toString(elementSystemType t)
{
  switch(t) {
    case MASS_0D:
      return "MASS_0D";
    case SOURCE_0D:
      return "SOURCE_0D";

    case MASS_1D:
      return "MASS_1D";
    case SOURCE_1D:
      return "SOURCE_1D";
    case DIFFUSION_1D:
      return "DIFFUSION_1D";
    case ADVECTION_1D:
      return "ADVECTION_1D";
    case DG_ADVECTION_1D:
      return "DG_ADVECTION_1D";
    case SUPG_STABILIZATION_1D:
      return "SUPG_STABILIZATION_1D";
    case NEUMANN_1D:
      return "NEUMANN_1D";

    case MASS_2D:
      return "MASS_2D";
    case SOURCE_2D:
      return "SOURCE_2D";
    case DIFFUSION_2D:
      return "DIFFUSION_2D";
    case ADVECTION_2D:
      return "ADVECTION_2D";
    case STOKES_2D:
      return "STOKES_2D";
    case NAVIERSTOKES_2D:
      return "NAVIERSTOKES_2D";
    default:
      return "[Unknown elementSystemType]";
  }
}

class feBilinearForm;

class feSysElm
{
protected:
  int _dim;
  std::vector<int> _iVar;
  std::vector<int> _jVar;
  elementSystemType _ID;

  // Compute the element FE matrix using finite differences

  bool _computeMatrixWithFD = false;
  // Does this form have an elementary matrix to assemble?
  bool _hasMatrix = false;

public:
  feSysElm(bool hasMatrix) : _hasMatrix(hasMatrix)
  {
    _iVar.resize(1);
    _jVar.resize(1);
  };
  virtual ~feSysElm() {}

  int getDim(){ return _dim; }
  std::vector<int> &getIVar() { return _iVar; }
  std::vector<int> &getJVar() { return _jVar; }

  elementSystemType getID() { return _ID; }
  std::string getIDName() { return toString(_ID); }

  bool computeMatrixWithFD() { return _computeMatrixWithFD; };

  bool hasMatrix() { return _hasMatrix; }

  virtual void createElementarySystem(std::vector<feSpace *> &space) = 0;
  virtual void computeAe(feBilinearForm *form) = 0;
  virtual void computeBe(feBilinearForm *form) = 0;
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
    : feSysElm(true), _fct(fct), _par(par)
  {
    _dim = 0;
    _ID = STIFFSPRING_0D;
  };
  virtual ~feSysElm_0D_StiffSpring() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
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
  feSysElm_0D_Stiff2(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par)
  {
    _dim = 0;
    _ID = STIFF2_0D;
  };
  virtual ~feSysElm_0D_Stiff2() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
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
  feSysElm_0D_Stiff3(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par)
  {
    _dim = 0;
    _ID = STIFF3_0D;
  };
  virtual ~feSysElm_0D_Stiff3() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
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
  feSysElm_0D_weakBC(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par)
  {
    _dim = 0;
    _ID = WEAKBC_0D;
  };
  virtual ~feSysElm_0D_weakBC() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
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
  feSysElm_0D_weakBC_edo1(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par)
  {
    _dim = 0;
    _ID = WEAKBC_EDO1_0D;
  };
  virtual ~feSysElm_0D_weakBC_edo1() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
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
  feSysElm_0D_weakBC_edo1_V2(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par)
  {
    _dim = 0;
    _ID = WEAKBC_EDO1_V2_0D;
  };
  virtual ~feSysElm_0D_weakBC_edo1_V2() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
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
  feSysElm_0D_weakBC_edo2(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par)
  {
    _dim = 0;
    _ID = WEAKBC_EDO2_0D;
  };
  virtual ~feSysElm_0D_weakBC_edo2() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
};

class feSysElm_0D_Masse : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;

public:
  feSysElm_0D_Masse(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par)
  {
    _dim = 0;
    _ID = MASS_0D;
  };
  virtual ~feSysElm_0D_Masse() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
};

class feSysElm_0D_Source : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;

public:
  feSysElm_0D_Source(double par, feFunction *fct) : feSysElm(false), _fct(fct), _par(par)
  {
    _dim = 0;
    _ID = SOURCE_0D;
  };
  virtual ~feSysElm_0D_Source() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
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
    : feSysElm(false), _fct(fct), _par(par)
  {
    _dim = 0;
    _ID = SOURCE_CROSSED_0D;
  };
  virtual ~feSysElm_0D_Source_crossed() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
};

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
  feSysElm_1D_weakBC_edo1(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par)
  {
    _dim = 1;
    _ID = WEAKBC_EDO1_1D;
  };
  virtual ~feSysElm_1D_weakBC_edo1() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
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
  feSysElm_1D_Source(double par, feFunction *fct) : feSysElm(false), _fct(fct), _par(par)
  {
    _dim = 1;
    _ID = SOURCE_1D;
  };
  virtual ~feSysElm_1D_Source() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form){};
  virtual void computeBe(feBilinearForm *form);
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
  feSysElm_1D_Diffusion(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par)
  {
    _dim = 1;
    _ID = DIFFUSION_1D;
    _computeMatrixWithFD = true;
  };
  virtual ~feSysElm_1D_Diffusion() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
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
  feSysElm_1D_Masse(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par)
  {
    _dim = 1;
    _ID = MASS_1D;
    _computeMatrixWithFD = false;
  };
  virtual ~feSysElm_1D_Masse() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
};

class feSysElm_1D_NeumannBC : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametres
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;

public:
  feSysElm_1D_NeumannBC(double par, feFunction *fct) : feSysElm(false), _fct(fct), _par(par)
  {
    _dim = 1;
    _ID = NEUMANN_1D;
  };
  virtual ~feSysElm_1D_NeumannBC() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form){};
  virtual void computeBe(feBilinearForm *form);
};

class feSysElm_1D_Advection : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametres
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;

public:
  feSysElm_1D_Advection(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par)
  {
    _dim = 1;
    _ID = ADVECTION_1D;
    _computeMatrixWithFD = false;
  };
  virtual ~feSysElm_1D_Advection() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
};

class feSysElm_1D_SUPGStab : public feSysElm
{
protected:
  feFunction *_fct;
  std::vector<double> _par;
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;

public:
  feSysElm_1D_SUPGStab(std::vector<double> par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par)
  {
    _dim = 1;
    _ID = SUPG_STABILIZATION_1D;
    _computeMatrixWithFD = true;
  };
  virtual ~feSysElm_1D_SUPGStab() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
};

class feSysElm_1D_DG_Advection : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametres
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;

public:
  feSysElm_1D_DG_Advection(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par)
  {
    _dim = 1;
    _ID = DG_ADVECTION_1D;
    _computeMatrixWithFD = false;
  };
  virtual ~feSysElm_1D_DG_Advection() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
};

class feSysElm_2D_Masse : public feSysElm
{
protected:
  feFunction *_fct;
  double _par; // Parametre
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;
  std::vector<double> _feUdy;
  std::vector<double> _dxdr; // [dx/dr, dy/dr, dz/dr]
  std::vector<double> _dxds; // [dx/ds, dy/ds, dz/ds]

public:
  feSysElm_2D_Masse(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par)
  {
    _dim = 2;
    _ID = MASS_2D;
    _computeMatrixWithFD = false;
  };
  virtual ~feSysElm_2D_Masse() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
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
  std::vector<double> _x;

public:
  feSysElm_2D_Source(double par, feFunction *fct) : feSysElm(false), _fct(fct), _par(par)
  {
    _dim = 2;
    _ID = SOURCE_2D;
  };

  virtual ~feSysElm_2D_Source() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form){};
  virtual void computeBe(feBilinearForm *form);
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
  std::vector<double> _dxdr; // [dx/dr, dy/dr, dz/dr]
  std::vector<double> _dxds; // [dx/ds, dy/ds, dz/ds]

public:
  feSysElm_2D_Diffusion(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par)
  {
    _dim = 2;
    _ID = DIFFUSION_2D;
    _computeMatrixWithFD = false;
  };

  virtual ~feSysElm_2D_Diffusion() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
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
  std::vector<double> _dxdr; 
  std::vector<double> _dxds; 

public:
  feSysElm_2D_Advection(double par, feVectorFunction *fct) : feSysElm(true), _fct(fct), _par(par)
  {
    _dim = 2;
    _ID = ADVECTION_2D;
    _computeMatrixWithFD = true;
  };
  virtual ~feSysElm_2D_Advection() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
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
  std::vector<double> _dxdr; 
  std::vector<double> _dxds; 

public:
  feSysElm_2D_Stokes(std::vector<double> &par, feVectorFunction *fct)
    : feSysElm(true), _fct(fct), _par(par)
  {
    _dim = 2;
    _ID = STOKES_2D;
    _computeMatrixWithFD = false;
  };
  virtual ~feSysElm_2D_Stokes() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
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
  std::vector<double> _dxdr;
  std::vector<double> _dxds;

public:
  feSysElm_2D_NavierStokes(std::vector<double> &par, feVectorFunction *fct,
                           feFunction *viscosityFct = nullptr)
    : feSysElm(true), _fct(fct), _par(par), _viscosityFct(viscosityFct)
  {
    _dim = 2;
    _ID = NAVIERSTOKES_2D;
    _computeMatrixWithFD = false;
  };
  virtual ~feSysElm_2D_NavierStokes() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
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
    _dim = 1;
    _ID = DIRECTIONALDONOTHING_1D;
    _computeMatrixWithFD = true;
  };
  virtual ~feSysElm_1D_DirectionalDoNothing() {}

  virtual void createElementarySystem(std::vector<feSpace *> &space);
  virtual void computeAe(feBilinearForm *form);
  virtual void computeBe(feBilinearForm *form);
};

#endif
