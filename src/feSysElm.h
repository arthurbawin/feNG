#ifndef _FESYSELM_
#define _FESYSELM_

#include "feFunction.h"
#include "feSpace.h"

typedef enum {

  // 0D weak forms
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

  // 1D weak forms
  WEAKBC_EDO1_1D,
  SOURCE_1D,
  DIFFUSION_1D,
  MASS_1D,
  NEUMANN_1D,
  ADVECTION_1D,
  DG_ADVECTION_1D,
  SUPG_STABILIZATION_1D,
  BEAM_1D,

  // 2D weak forms
  MASS_2D,
  SOURCE_2D,
  DIFFUSION_2D,
  ADVECTION_2D,
  STOKES_2D,
  NAVIERSTOKES_2D,
} elementSystemType;

inline const std::string toString(elementSystemType t)
{
  switch(t) {

    case MASS_0D:
      return "MASS_0D";
    case SOURCE_0D:
      return "SOURCE_0D";

    case SOURCE_1D:
      return "SOURCE_1D";
    case DIFFUSION_1D:
      return "DIFFUSION_1D";
    case MASS_1D:
      return "MASS_1D";
      case NEUMANN_1D:
      return "NEUMANN_1D";
    case ADVECTION_1D:
      return "ADVECTION_1D";
    case DG_ADVECTION_1D:
      return "DG_ADVECTION_1D";
    case SUPG_STABILIZATION_1D:
      return "SUPG_STABILIZATION_1D";
    case BEAM_1D:
      return "BEAM_1D";
    

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

//
// Abstract class for the computation of weak forms.
//
class feSysElm
{
protected:
  // Weak form identifier
  elementSystemType _ID;
  // Dimension
  int _dim;
  // Number of FE fields in the weak form
  int _nFields;
  // Fields layout, see feBilinearForm.h for details
  std::vector<int> _fieldsLayoutI;
  std::vector<int> _fieldsLayoutJ;
  // Number of shape functions on the element
  int _nFunctions;

  // Compute the element FE matrix using finite differences?
  bool _computeMatrixWithFD = false;
  // Does this form have an elementary matrix to assemble?
  bool _hasMatrix = false;

  // Physical coordinates of a node
  std::vector<double> _pos;

public:
  // Quadrature rule
  int _nQuad;
  std::vector<double> _wQuad;
  std::vector<double> _rQuad;
  std::vector<double> _sQuad;
  std::vector<double> _tQuad;

public:
  // Do not call directly, call derived classes instead
  feSysElm(bool hasMatrix) : _hasMatrix(hasMatrix), _pos(3, 0.)
  {
    _fieldsLayoutI.resize(1);
    _fieldsLayoutJ.resize(1);
  };
  virtual ~feSysElm() {}

  int getDim(){ return _dim; }
  int getNumFields(){ return _nFields; }
  std::vector<int> &getFieldLayoutI() { return _fieldsLayoutI; }
  std::vector<int> &getFieldLayoutJ() { return _fieldsLayoutJ; }
  elementSystemType getID() { return _ID; }
  std::string getWeakFormName() { return toString(_ID); }
  bool computeMatrixWithFD() { return _computeMatrixWithFD; }
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
      _nFields = 2;
      _ID = STIFFSPRING_0D;
    };
    ~feSysElm_0D_StiffSpring() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
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
      _nFields = 3;
      _ID = STIFF2_0D;
    };
    ~feSysElm_0D_Stiff2() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
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
      _nFields = 3;
      _ID = STIFF3_0D;
    };
    ~feSysElm_0D_Stiff3() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
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
      _nFields = 2;
      _ID = WEAKBC_0D;
    };
    ~feSysElm_0D_weakBC() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
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
      _nFields = 3;
      _ID = WEAKBC_EDO1_0D;
    };
    ~feSysElm_0D_weakBC_edo1() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
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
      _nFields = 3;
      _ID = WEAKBC_EDO1_V2_0D;
    };
    ~feSysElm_0D_weakBC_edo1_V2() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
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
      _nFields = 4;
      _ID = WEAKBC_EDO2_0D;
    };
    ~feSysElm_0D_weakBC_edo2() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
};

class feSysElm_0D_Masse : public feSysElm
{
  protected:
    feFunction *_fct;
    double _par;
    int _idU;
    std::vector<double> _feU;

  public:
    feSysElm_0D_Masse(double par, feFunction *fct) : feSysElm(true), _fct(fct), _par(par)
    {
      _dim = 0;
      _nFields = 1;
      _ID = MASS_0D;
    };
    ~feSysElm_0D_Masse() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
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
      _nFields = 1;
      _ID = SOURCE_0D;
    };
    ~feSysElm_0D_Source() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
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
      _nFields = 1;
      _ID = SOURCE_CROSSED_0D;
    };
    ~feSysElm_0D_Source_crossed() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
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
      _nFields = 3;
      _ID = WEAKBC_EDO1_1D;
    };
    ~feSysElm_1D_weakBC_edo1() {}
    void createElementarySystem(std::vector<feSpace *> &space);
    void computeAe(feBilinearForm *form);
    void computeBe(feBilinearForm *form);
};

//
// 1D Source term weak form.
// Residual only
//
//  /
//  | fv dx
//  /
//
// # fields: 1 (test functions)
//
// Fields layout: phi_U []
//
class feSysElm_1D_Source : public feSysElm
{
protected:
  feFunction *_sourceFun;
  int _idU;
  std::vector<double> _feU;
public:
  feSysElm_1D_Source(feFunction *sourceFun) : feSysElm(false), _sourceFun(sourceFun)
  {
    _dim = 1;
    _nFields = 1;
    _ID = SOURCE_1D;
  };
  ~feSysElm_1D_Source() {}
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form){};
  void computeBe(feBilinearForm *form);
};

//
// 1D Diffusion term weak form (stiffness matrix).
// Matrix and residual
//
//  /
//  | kD * grad(u) dot grad(v) dx 
//  /
//
// # fields: 1 (FE solution and test functions)
//                        U
// Fields layout: phi_U [   ]
//
class feSysElm_1D_Diffusion : public feSysElm
{
protected:
  feFunction *_diffusivityFun;
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;
public:
  feSysElm_1D_Diffusion(feFunction *diffusivityFun) : feSysElm(true), _diffusivityFun(diffusivityFun)
  {
    _dim = 1;
    _nFields = 1;
    _ID = DIFFUSION_1D;
     _computeMatrixWithFD = false;
  };
  ~feSysElm_1D_Diffusion() {}
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
};

//
// 1D Transient term weak form (mass matrix).
// Matrix and residual
//
//  /
//  | coeff * dudt * v dx 
//  /
//
// # fields: 1 (FE solution and test functions)
//                        U
// Fields layout: phi_U [   ]
//
class feSysElm_1D_Masse : public feSysElm
{
protected:
  feFunction *_fct;
  double _coeff;
  int _idU;
  std::vector<double> _feU;

public:
  feSysElm_1D_Masse(double coeff, feFunction *fct) : feSysElm(true), _fct(fct), _coeff(coeff)
  {
    _dim = 1;
    _nFields = 1;
    _ID = MASS_1D;
    _computeMatrixWithFD = false;
  };
  ~feSysElm_1D_Masse() {}
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
};

//
// 1D Neumann boundary condition weak form.
// Residual only
//
//  /                         /
//  | (grad(u) dot n) v dx =  | h v dx
//  /Boundary                 /Boundary
//
// where (grad(u) dot n) = h on this boundary.
//
// # fields: 1 (test functions)
//                        
// Fields layout: phi_U [   ]
//
class feSysElm_1D_NeumannBC : public feSysElm
{
protected:
  feFunction *_neumannBC;
  int _idU;
  std::vector<double> _feU;

public:
  feSysElm_1D_NeumannBC(feFunction *neumannBC) : feSysElm(false), _neumannBC(neumannBC)
  {
    _dim = 1;
    _nFields = 1;
    _ID = NEUMANN_1D;
  };
  ~feSysElm_1D_NeumannBC() {}
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form){};
  void computeBe(feBilinearForm *form);
};

//
// 1D Advection term weak form
// Matrix and residual
//
//    /
//  - | c * u * dvdx dx
//    /
//
// with c an imposed velocity field
//
// # fields: 1 (FE solution and test functions)
//                        U
// Fields layout: phi_U [   ]
//
class feSysElm_1D_Advection : public feSysElm
{
protected:
  feFunction *_cVelocity;
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;

public:
  feSysElm_1D_Advection(feFunction *cVelocity) : feSysElm(true), _cVelocity(cVelocity)
  {
    _dim = 1;
    _nFields = 1;
    _ID = ADVECTION_1D;
    _computeMatrixWithFD = false;
  };
  ~feSysElm_1D_Advection() {}
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
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
    _nFields = 1;
    _ID = SUPG_STABILIZATION_1D;
    _computeMatrixWithFD = true;
  };
  ~feSysElm_1D_SUPGStab() {}

  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
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
    _nFields = 1;
    _ID = DG_ADVECTION_1D;
    _computeMatrixWithFD = false;
  };
  ~feSysElm_1D_DG_Advection() {}

  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
};

class feSysElm_1D_Beam : public feSysElm
{
protected:
  feFunction *_fct;
  double _par;
  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx2;

public:
  feSysElm_1D_Beam(double EI) : feSysElm(true), _par(EI)
  {
    _dim = 1;
    _nFields = 1;
    _ID = BEAM_1D;
    _computeMatrixWithFD = true;
  };
  ~feSysElm_1D_Beam() {}

  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
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
    _nFields = 1;
    _ID = MASS_2D;
    _computeMatrixWithFD = false;
  };
  ~feSysElm_2D_Masse() {}

  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
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
    _nFields = 1;
    _ID = SOURCE_2D;
  };

  ~feSysElm_2D_Source() {}

  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form){};
  void computeBe(feBilinearForm *form);
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
    _nFields = 1;
    _ID = DIFFUSION_2D;
    _computeMatrixWithFD = false;
  };

  ~feSysElm_2D_Diffusion() {}

  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
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
    _nFields = 1;
    _ID = ADVECTION_2D;
    _computeMatrixWithFD = true;
  };
  ~feSysElm_2D_Advection() {}

  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
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
    _nFields = 3;
    _ID = STOKES_2D;
    _computeMatrixWithFD = false;
  };
  ~feSysElm_2D_Stokes() {}

  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
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
    _nFields = 1;
    _ID = NAVIERSTOKES_2D;
    _computeMatrixWithFD = false;
  };
  ~feSysElm_2D_NavierStokes() {}

  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
};

#endif
