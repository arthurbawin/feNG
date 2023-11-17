#ifndef _FESYSELM_
#define _FESYSELM_

#include "feFunction.h"
#include "feSpace.h"

#define CLONEABLE(Type)                                                                            \
  virtual Type *clone() const { return new Type(*this); }

typedef enum {
  SOURCE,
  SOURCE_DIRAC,
  VECTOR_SOURCE,
  MIXED_SCALAR_VECTOR_SOURCE,
  GRAD_SOURCE,
  MASS,
  MIXED_MASS,
  VECTOR_MASS,
  MIXED_VECTOR_MASS,
  MIXED_SCALAR_VECTOR_MASS,
  TRANSIENT_MASS,
  TRANSIENT_VECTOR_MASS,
  DIFFUSION,
  VECTOR_DIFFUSION,
  NONLINEAR_DIFFUSION,
  ADVECTION,
  NONLINEAR_ADVECTION,
  VECTOR_CONVECTIVE_ACCELERATION,
  VECTOR_ADJOINT_CONVECTIVE_ACCELERATION,
  DIV_NEWTONIAN_STRESS,
  MIXED_GRADIENT,
  MIXED_DIVERGENCE,
  MIXED_CURL,
  ZERO_BLOCK,

  SUPG_PSPG_STOKES,
  SUPG_PSPG_NAVIERSTOKES,

  GLS_STOKES,
  GLS_NAVIERSTOKES,

  // 1D weak forms
  NEUMANN_1D,
  SUPG_STABILIZATION_1D,

  // 2D weak forms
  SUPG_STABILIZATION_2D,

} elementSystemType;

inline const std::string toString(elementSystemType t)
{
  switch(t) {
    case SOURCE:
      return "SOURCE";
    case SOURCE_DIRAC:
      return "SOURCE_DIRAC";
    case VECTOR_SOURCE:
      return "VECTOR_SOURCE";
    case MIXED_SCALAR_VECTOR_SOURCE:
      return "MIXED_SCALAR_VECTOR_SOURCE";
    case GRAD_SOURCE:
      return "GRAD_SOURCE";
    case MASS:
      return "MASS";
    case MIXED_MASS:
      return "MIXED_MASS";
    case VECTOR_MASS:
      return "VECTOR_MASS";
    case MIXED_VECTOR_MASS:
      return "MIXED_VECTOR_MASS";
    case MIXED_SCALAR_VECTOR_MASS:
      return "MIXED_SCALAR_VECTOR_MASS";
    case TRANSIENT_MASS:
      return "TRANSIENT_MASS";
    case TRANSIENT_VECTOR_MASS:
      return "TRANSIENT_VECTOR_MASS";
    case DIFFUSION:
      return "DIFFUSION";
    case VECTOR_DIFFUSION:
      return "VECTOR_DIFFUSION";
    case NONLINEAR_DIFFUSION:
      return "NONLINEAR_DIFFUSION";
    case ADVECTION:
      return "ADVECTION";
    case NONLINEAR_ADVECTION:
      return "NONLINEAR_ADVECTION";
    case VECTOR_CONVECTIVE_ACCELERATION:
      return "VECTOR_CONVECTIVE_ACCELERATION";
    case VECTOR_ADJOINT_CONVECTIVE_ACCELERATION:
      return "VECTOR_ADJOINT_CONVECTIVE_ACCELERATION";
    case DIV_NEWTONIAN_STRESS:
      return "DIV_NEWTONIAN_STRESS";
    case MIXED_GRADIENT:
      return "MIXED_GRADIENT";
    case MIXED_DIVERGENCE:
      return "MIXED_DIVERGENCE";
    case MIXED_CURL:
      return "MIXED_CURL";
    case ZERO_BLOCK:
      return "ZERO_BLOCK";
    case NEUMANN_1D:
      return "NEUMANN_1D";
    case SUPG_STABILIZATION_1D:
      return "SUPG_STABILIZATION_1D";
    case SUPG_STABILIZATION_2D:
      return "SUPG_STABILIZATION_2D";
    case SUPG_PSPG_STOKES:
      return "SUPG_PSPG_STOKES";
    case SUPG_PSPG_NAVIERSTOKES:
      return "SUPG_PSPG_NAVIERSTOKES";
    case GLS_STOKES:
      return "GLS_STOKES";
    case GLS_NAVIERSTOKES:
      return "GLS_NAVIERSTOKES";
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
  // Number of vector components on the FE space
  int _nComponents;

  // Compute the element FE matrix using finite differences?
  bool _computeMatrixWithFD = false;
  // Does this form have an elementary matrix to assemble?
  bool _hasMatrix = false;
  // Is the computation thread safe?
  bool _isThreadSafe = true;

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
  feSysElm(int dim, int nFields, elementSystemType ID, bool hasMatrix)
    : _dim(dim), _nFields(nFields), _ID(ID), _hasMatrix(hasMatrix), _pos(3, 0.)
  {
    _fieldsLayoutI.resize(1);
    _fieldsLayoutJ.resize(1);
  };
  virtual ~feSysElm() {}

  virtual feSysElm *clone() const = 0;

  int getDim() { return _dim; }
  int getNumFields() { return _nFields; }
  std::vector<int> &getFieldLayoutI() { return _fieldsLayoutI; }
  std::vector<int> &getFieldLayoutJ() { return _fieldsLayoutJ; }
  elementSystemType getID() { return _ID; }
  std::string getWeakFormName() { return toString(_ID); }
  bool computeMatrixWithFD() { return _computeMatrixWithFD; }
  bool hasMatrix() { return _hasMatrix; }
  bool isThreadSafe() { return _isThreadSafe; }

  virtual void createElementarySystem(std::vector<feSpace *> &spaces){};
  virtual void computeAe(feBilinearForm *form){};
  virtual void computeBe(feBilinearForm *form) = 0;
};

//
// Source term weak form.
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
class feSysElm_Source : public feSysElm
{
protected:
  feFunction *_source;
  int _idU;
  std::vector<double> _phiU;

public:
  feSysElm_Source(feFunction *source) : feSysElm(-1, 1, SOURCE, false), _source(source){};
  ~feSysElm_Source(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_Source)
};

//
// A particular case where the source term is the Dirac delta.
// Used to compute the adjoint solution when the output functional
// is the solution evaluated at x_0.
//
// !!! NOT THREAD SAFE BECAUSE OF THE LOCATE VERTEX CALLS
//
class feSysElm_SourceDirac : public feSysElm
{
protected:
  std::vector<double> &_x0;
  int _idU;
  std::vector<double> _phiAtX0;

public:
  feSysElm_SourceDirac(std::vector<double> &x0) : feSysElm(-1, 1, SOURCE, false), _x0(x0)
  { _isThreadSafe = false; };
  ~feSysElm_SourceDirac(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_SourceDirac)
};

class feSysElm_VectorSource : public feSysElm
{
protected:
  feVectorFunction *_source;
  std::vector<double> _S;
  int _idU;
  std::vector<double> _phiU;

public:
  feSysElm_VectorSource(feVectorFunction *source)
    : feSysElm(-1, 1, VECTOR_SOURCE, false), _source(source), _S(3, 0.){};
  ~feSysElm_VectorSource(){};
  void createElementarySystem(std::vector<feSpace *> &spaces);
  // void createElementarySystem(std::vector<feVectorSpace *> &vectorSpaces){};
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_VectorSource)
};

//
// Reaction weak form (mass matrix).
// Matrix and residual
//
//  /
//  | coeff * u * v dx
//  /
//
// where v is the test function of u.
//
// # fields: 1 (FE solution and test functions)
//                        U
// Fields layout: phi_U [   ]
//
class feSysElm_Mass : public feSysElm
{
protected:
  feFunction *_coeff;
  int _idU;
  std::vector<double> _phiU;

public:
  feSysElm_Mass(feFunction *coeff) : feSysElm(-1, 1, MASS, true), _coeff(coeff){};
  ~feSysElm_Mass(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_Mass)
};

//
// Mass matrix for one variable tested with the test functions of another variable
// Matrix and residual
//
//  /
//  | coeff * u * v dx
//  /
//
// where v is NOT the test function of u, but of another variable.
//
// # fields: 2 (FE solution and distinct test functions)
//                    U
// Fields layout: V [   ]
//
class feSysElm_MixedMass : public feSysElm
{
protected:
  feFunction *_coeff;
  int _idU, _idV, _nFunctionsU, _nFunctionsV;
  std::vector<double> _phiU, _phiV;

public:
  feSysElm_MixedMass(feFunction *coeff) : feSysElm(-1, 2, MIXED_MASS, true), _coeff(coeff){};
  ~feSysElm_MixedMass(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_MixedMass)
};

class feSysElm_VectorMass : public feSysElm
{
protected:
  feFunction *_coeff;
  std::vector<double> _u;
  int _idU;
  std::vector<double> _phiU;

public:
  feSysElm_VectorMass(feFunction *coeff) : feSysElm(-1, 1, VECTOR_MASS, true), _coeff(coeff){};
  ~feSysElm_VectorMass(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_VectorMass)
};

//
// Mass matrix for one vector-valued variable tested with the test functions 
// of another vector-valued variable"
// Matrix and residual
//
//  /
//  | coeff * u dot v dx
//  /
//
// where v is NOT the test function of u, but of another variable.
//
// # fields: 2 (FE solution and distinct test functions)
//                    U
// Fields layout: V [   ]
//
class feSysElm_MixedVectorMass : public feSysElm
{
protected:
  feFunction *_coeff;
  std::vector<double> _u;
  int _idU, _idV, _nFunctionsU, _nFunctionsV, _nComponentsU, _nComponentsV;
  std::vector<double> _phiU, _phiV;

public:
  feSysElm_MixedVectorMass(feFunction *coeff) : feSysElm(-1, 2, MIXED_VECTOR_MASS, true), _coeff(coeff){};
  ~feSysElm_MixedVectorMass(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_MixedVectorMass)
};

//
// Mass matrix for one vector-valued variable tested with the test functions of 
// a scalar-valued variable.
//
// Matrix and residual
//
//  /
//  | coeff * u * v dx with u vector-valued, v scalar 
//  /
//
// # fields: 2 (FE solution and distinct test functions)
//                    U
// Fields layout: V [   ]
//
class feSysElm_MixedScalarVectorMass : public feSysElm
{
protected:
  feFunction *_coeff;
  std::vector<double> _u;
  int _idU, _idV, _nFunctionsU, _nFunctionsV, _nComponentsU;
  std::vector<double> _phiU, _phiV;

public:
  feSysElm_MixedScalarVectorMass(feFunction *coeff) : feSysElm(-1, 2, MIXED_SCALAR_VECTOR_MASS, true), _coeff(coeff)
  {
    _computeMatrixWithFD = true;
  };
  ~feSysElm_MixedScalarVectorMass(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_MixedScalarVectorMass)
};

//
// Transient term weak form (transient mass matrix).
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
class feSysElm_TransientMass : public feSysElm
{
protected:
  double _coeff;
  int _idU;
  std::vector<double> _feU;

public:
  feSysElm_TransientMass(double coeff) : feSysElm(-1, 1, TRANSIENT_MASS, true), _coeff(coeff){
    _computeMatrixWithFD = true;
  };
  ~feSysElm_TransientMass(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_TransientMass)
};

//
// Transient term weak form (transient mass matrix) for vector-valued field
// Matrix and residual
//
//  /
//  | coeff * dudt dot v dx
//  /
//
// # fields: 1 (FE solution and test functions)
//                        U
// Fields layout: phi_U [   ]
//
class feSysElm_TransientVectorMass : public feSysElm
{
protected:
  double _coeff;
  int _idU;
  std::vector<double> _dudt, _phiU;

public:
  feSysElm_TransientVectorMass(double coeff) : feSysElm(-1, 1, TRANSIENT_VECTOR_MASS, true), _coeff(coeff)
  {
    _computeMatrixWithFD = true;
  };
  ~feSysElm_TransientVectorMass(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_TransientVectorMass)
};


//
// Diffusion term weak form (stiffness matrix).
// Matrix and residual
//
//  /
//  | k * grad(u) dot grad(v) dx
//  /
//
// # fields: 1 (FE solution and test functions)
//
// Parameters:
//  diffusivity: scalar field returning the diffusivity k(x)
//
//                        U
// Fields layout: phi_U [   ]
//
template <int dim> class feSysElm_Diffusion : public feSysElm
{
protected:
  feFunction *_coeff;
  int _idU;
  std::vector<double> _gradPhi;

public:
  feSysElm_Diffusion(feFunction *coeff) : feSysElm(dim, 1, DIFFUSION, true), _coeff(coeff){};
  ~feSysElm_Diffusion(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_Diffusion)
};

class feSysElm_VectorDiffusion : public feSysElm
{
protected:
  feFunction *_coeff;
  feFunction *_diffusivity;
  int _idU;
  std::vector<double> _gradu;
  std::vector<double> _gradPhi;

public:
  feSysElm_VectorDiffusion(feFunction *coeff, feFunction *diffusivity)
    : feSysElm(-1, 1, VECTOR_DIFFUSION, true), _coeff(coeff), _diffusivity(diffusivity)
  {
    _computeMatrixWithFD = false;
  };
  ~feSysElm_VectorDiffusion(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_VectorDiffusion)
};

template <int dim> class feSysElm_NonlinearDiffusion : public feSysElm
{
protected:
  feFunction *_diffusivity;
  feFunction *_ddiffdu;
  int _idU;
  std::vector<double> _phiU;
  std::vector<double> _gradPhi;

public:
  feSysElm_NonlinearDiffusion(feFunction *diffusivity, feFunction *ddiffdu)
    : feSysElm(dim, 1, NONLINEAR_DIFFUSION, true), _diffusivity(diffusivity), _ddiffdu(ddiffdu)
  {
    _computeMatrixWithFD = false;
  };
  ~feSysElm_NonlinearDiffusion(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_NonlinearDiffusion)
};

//
// Advection term weak form
// Matrix and residual
//
//    /
//  - | (c dot grad(u)) * v dx
//    /
//
// with c an imposed velocity field
//
// # fields: 1 (FE solution and test functions)
//                        U
// Fields layout: phi_U [   ]
//
template <int dim> class feSysElm_Advection : public feSysElm
{
protected:
  // Imposed velocity field
  feVectorFunction *_velocity;
  int _idU;
  std::vector<double> _phiU;
  std::vector<double> _gradPhi;

public:
  feSysElm_Advection(feVectorFunction *velocity)
    : feSysElm(-1, 1, ADVECTION, true), _velocity(velocity)
  {
    _computeMatrixWithFD = true;
  };
  ~feSysElm_Advection(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_Advection)
};

//
// Nonlinear advection term weak form
// Matrix and residual
//
// PDE term: div(f(u)) * v
//
// Term in weak form:
//    /
//  - | f(u) cdot grad(v) dx
//    /
//
// with f(u) a possibly nonlinear vector flux (e.g. u*u/2)
//
// # fields: 1 (FE solution and test functions)
//                        U
// Fields layout: phi_U [   ]
//
template <int dim> class feSysElm_NonlinearAdvection : public feSysElm
{
protected:
  // Nonlinear flux
  feVectorFunction *_flux;
  int _idU;
  std::vector<double> _phiU;
  std::vector<double> _gradPhi;

public:
  feSysElm_NonlinearAdvection(feVectorFunction *flux)
    : feSysElm(-1, 1, NONLINEAR_ADVECTION, true), _flux(flux)
  {
    _computeMatrixWithFD = true;
  };
  ~feSysElm_NonlinearAdvection(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_NonlinearAdvection)
};

//
// Integral of coeff * (u dot grad) u dot v
//
class feSysElm_VectorConvectiveAcceleration : public feSysElm
{
protected:
  feFunction *_coeff;
  int _idU;
  std::vector<double> _u;
  std::vector<double> _gradu;
  std::vector<double> _uDotGradu;
  std::vector<double> _phiU;
  std::vector<double> _gradPhiU;

public:
  feSysElm_VectorConvectiveAcceleration(feFunction *coeff)
    : feSysElm(-1, 1, VECTOR_CONVECTIVE_ACCELERATION, true), _coeff(coeff)
  {
    _computeMatrixWithFD = false;
  };
  ~feSysElm_VectorConvectiveAcceleration(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_VectorConvectiveAcceleration)
};

//
// Integral of coeff * (u dot grad uAd) dot phi_uAd
//
// where u is the velocity field
//       uAd is the adjoint velocity field
//       phi_uAd is the test function for the uAd
//
// Convective term in the adjoint Navier-Stokes equations
//
class feSysElm_VectorAdjointConvectiveAcceleration : public feSysElm
{
protected:
  feFunction *_coeff;
  int _idU, _idUad;
  int _nFunctionsU, _nFunctionsUad;
  std::vector<double> _u;
  std::vector<double> _grad_uAd;
  std::vector<double> _uDotGraduAd;
  std::vector<double> _phi_uAd;
  std::vector<double> _gradPhiU;

public:
  feSysElm_VectorAdjointConvectiveAcceleration(feFunction *coeff)
    : feSysElm(-1, 2, VECTOR_ADJOINT_CONVECTIVE_ACCELERATION, true), _coeff(coeff)
  {
    _computeMatrixWithFD = true;
  };
  ~feSysElm_VectorAdjointConvectiveAcceleration(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_VectorAdjointConvectiveAcceleration)
};

//
// Integral of - sigma : grad(v), with sigma = -pI + 2*mu*d(u)
//                                with d(u) = (grad(u) + grad(u)^T)/2
//
// Hence: integral of -p * div(v) + mu * d : grad(v) dx
//
class feSysElm_DivergenceNewtonianStress : public feSysElm
{
protected:
  feFunction *_coeff;
  feFunction *_viscosity;
  int _idU;
  int _idP;
  std::vector<double> _gradu;
  std::vector<double> _symmetricGradu;
  std::vector<double> _gradPhiU;

public:
  feSysElm_DivergenceNewtonianStress(feFunction *coeff, feFunction *viscosity)
    : feSysElm(-1, 2, DIV_NEWTONIAN_STRESS, true), _coeff(coeff), _viscosity(viscosity)
  {
    _computeMatrixWithFD = true;
  };
  ~feSysElm_DivergenceNewtonianStress(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_DivergenceNewtonianStress)
};

// Weak form of grad(f) cdot v, with f scalar and v vector-valued
//
// Integral of -f div(v)

class feSysElm_GradSource : public feSysElm
{
protected:
  feFunction *_source;
  int _idU;
  std::vector<double> _gradPhi;

public:
  feSysElm_GradSource(feFunction *source) : feSysElm(-1, 1, GRAD_SOURCE, false), _source(source){};
  ~feSysElm_GradSource(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_GradSource)
};

//
// Mixed gradient : (coeff * grad(u)) cdot v, with v vector-valued test functions
//                                                 coeff a scalar coefficient
// Integral of - coeff * u * div(v)
//
class feSysElm_MixedGradient : public feSysElm
{
protected:
  feFunction *_coeff;
  int _idU;
  int _idV;
  int _nFunctionsU;
  int _nFunctionsV;
  std::vector<double> _phiU;
  std::vector<double> _gradPhiV;

public:
  feSysElm_MixedGradient(feFunction *coeff) : feSysElm(-1, 2, MIXED_GRADIENT, true), _coeff(coeff)
  {
    _computeMatrixWithFD = false;
  };
  ~feSysElm_MixedGradient(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_MixedGradient)
};

//
// Mixed divergence : (coeff * div(u)) * v, with v scalar test functions
//                                                 coeff a scalar coefficient
// Integral of (coeff * div(u)) * v
//
class feSysElm_MixedDivergence : public feSysElm
{
protected:
  feFunction *_coeff;
  int _idU;
  int _idV;
  int _nFunctionsU;
  int _nFunctionsV;
  std::vector<double> _gradu;
  std::vector<double> _gradPhiU;
  std::vector<double> _phiV;

public:
  feSysElm_MixedDivergence(feFunction *coeff)
    : feSysElm(-1, 2, MIXED_DIVERGENCE, true), _coeff(coeff)
  {
    _computeMatrixWithFD = false;
  };
  ~feSysElm_MixedDivergence(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_MixedDivergence)
};

class feSysElm_MixedCurl : public feSysElm
{
protected:
  feFunction *_coeff;
  int _idU;
  int _idV;
  int _nFunctionsU;
  int _nFunctionsV;
  std::vector<double> _gradu;
  std::vector<double> _gradPhiU;
  std::vector<double> _phiV;

public:
  feSysElm_MixedCurl(feFunction *coeff) : feSysElm(2, 2, MIXED_CURL, true), _coeff(coeff)
  {
    _computeMatrixWithFD = true;
  };
  ~feSysElm_MixedCurl(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_MixedCurl)
};

class feSysElm_ZeroBlock : public feSysElm
{
protected:
  int _idU;
  int _idV;

public:
  feSysElm_ZeroBlock() : feSysElm(-1, 2, ZERO_BLOCK, true) { _computeMatrixWithFD = false; };
  ~feSysElm_ZeroBlock(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeBe(feBilinearForm *form){};
  CLONEABLE(feSysElm_ZeroBlock)
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
  feSysElm_1D_NeumannBC(feFunction *neumannBC)
    : feSysElm(1, 1, NEUMANN_1D, false), _neumannBC(neumannBC){};
  ~feSysElm_1D_NeumannBC() {}
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form){};
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_1D_NeumannBC)
};

//
// 1D SUPG stabilization for the advection-diffusion equation.
// Matrix and residual
//
//  /
//  | P(v) * tau * R(u) dx
//  /
//
// with:
//  P(v) = c dot v
//  tau: the SUPG parameter defined by nuNum/c^2, nuNum the numerical diffusivity
//  R(u) = c*dudx + kd2u/dx2 + f the residual
//
// # fields: 1 (FE solution and test functions)
//
// Parameters:
//     velocity: imposed velocity (scalar field)
//  diffusivity: physical diffusivity (double)
//       source: source term (scalar field)
//                        U
// Fields layout: phi_U [   ]
//
class feSysElm_1D_SUPGStab : public feSysElm
{
protected:
  feFunction *_velocity;
  double _diffusivity;
  feFunction *_source;

  int _idU;
  std::vector<double> _feU;
  std::vector<double> _feUdx;

public:
  feSysElm_1D_SUPGStab(double diffusivity, feFunction *velocity, feFunction *source)
    : feSysElm(1, 1, SUPG_STABILIZATION_1D, true), _velocity(velocity), _diffusivity(diffusivity),
      _source(source)
  {
    _computeMatrixWithFD = true;
  };
  ~feSysElm_1D_SUPGStab() {}

  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_1D_SUPGStab)
};

class feSysElm_2D_SUPGStab : public feSysElm
{
protected:
  feVectorFunction *_velocity;
  feFunction *_diffusivity;
  feFunction *_reactionCoeff;
  feFunction *_source;

  int _idU;
  std::vector<double> _phiU;
  std::vector<double> _gradPhiU;

public:
  feSysElm_2D_SUPGStab(feFunction *reactionCoeff, feFunction *diffusivity,
                       feVectorFunction *velocity, feFunction *source)
    : feSysElm(2, 1, SUPG_STABILIZATION_2D, true), _reactionCoeff(reactionCoeff),
      _velocity(velocity), _diffusivity(diffusivity), _source(source)
  {
    _computeMatrixWithFD = true;
  };
  ~feSysElm_2D_SUPGStab() {}

  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_2D_SUPGStab)
};

//
// SUPG + PSPG stabilization for the (Navier-)Stokes equations
//
class feSysElm_Stokes_SUPG_PSPG : public feSysElm
{
protected:
  feFunction *_coeff;
  feFunction *_density;
  feFunction *_viscosity;
  feVectorFunction *_volumeForce;
  int _idU;
  int _idP;
  int _nFunctionsU;
  int _nFunctionsP;
  std::vector<double> _f;
  std::vector<double> _residual;
  std::vector<double> _u;
  std::vector<double> _gradu;
  std::vector<double> _hessu;
  std::vector<double> _gradp;
  std::vector<double> _gradPhiU;
  std::vector<double> _gradPhiP;
  std::vector<double> _hessPhiU;

public:
  feSysElm_Stokes_SUPG_PSPG(feFunction *coeff, feFunction *density, feFunction *viscosity,
                           feVectorFunction *volumeForce)
    : feSysElm(-1, 2, SUPG_PSPG_STOKES, true), _coeff(coeff), _density(density),
      _viscosity(viscosity), _volumeForce(volumeForce)
  {
    _computeMatrixWithFD = true;
  };
  ~feSysElm_Stokes_SUPG_PSPG(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_Stokes_SUPG_PSPG)
};

class feSysElm_NS_SUPG_PSPG : public feSysElm
{
protected:
  feFunction *_coeff;
  feFunction *_density;
  feFunction *_viscosity;
  feVectorFunction *_volumeForce;
  int _idU;
  int _idP;
  int _nFunctionsU;
  int _nFunctionsP;
  std::vector<double> _f;
  std::vector<double> _residual;
  std::vector<double> _u;
  std::vector<double> _gradu;
  std::vector<double> _hessu;
  std::vector<double> _gradp;
  std::vector<double> _gradPhiU;
  std::vector<double> _gradPhiP;
  std::vector<double> _hessPhiU;
  std::vector<double> _uDotGradu;
  std::vector<double> _uDotGradPhiu;

public:
  feSysElm_NS_SUPG_PSPG(feFunction *coeff, feFunction *density, feFunction *viscosity,
                           feVectorFunction *volumeForce)
    : feSysElm(-1, 2, SUPG_PSPG_NAVIERSTOKES, true), _coeff(coeff), _density(density),
      _viscosity(viscosity), _volumeForce(volumeForce)
  {
    _computeMatrixWithFD = true;
  };
  ~feSysElm_NS_SUPG_PSPG(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_NS_SUPG_PSPG)
};

//
// Galerkin least-squares stabilization for the (Navier-)Stokes equations
//
class feSysElm_GLS_Stokes : public feSysElm
{
protected:
  feFunction *_coeff;
  feFunction *_density;
  feFunction *_viscosity;
  feVectorFunction *_volumeForce;
  int _idU;
  int _idP;
  int _nFunctionsU;
  int _nFunctionsP;
  std::vector<double> _f;
  std::vector<double> _residual;
  std::vector<double> _u;
  std::vector<double> _gradu;
  std::vector<double> _hessu;
  std::vector<double> _gradp;
  std::vector<double> _gradPhiU;
  std::vector<double> _hessPhiU;
  std::vector<double> _gradPhiP;

public:
  feSysElm_GLS_Stokes(feFunction *coeff, feFunction *density, feFunction *viscosity,
                           feVectorFunction *volumeForce)
    : feSysElm(-1, 2, GLS_STOKES, true), _coeff(coeff), _density(density),
      _viscosity(viscosity), _volumeForce(volumeForce)
  {
    _computeMatrixWithFD = true;
  };
  ~feSysElm_GLS_Stokes(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_GLS_Stokes)
};

class feSysElm_GLS_NavierStokes : public feSysElm
{
protected:
  feFunction *_coeff;
  feFunction *_density;
  feFunction *_viscosity;
  feVectorFunction *_volumeForce;
  int _idU;
  int _idP;
  int _nFunctionsU;
  int _nFunctionsP;
  std::vector<double> _f;
  std::vector<double> _residual;
  std::vector<double> _u;
  std::vector<double> _gradu;
  std::vector<double> _hessu;
  std::vector<double> _gradp;
  std::vector<double> _gradPhiU;
  std::vector<double> _hessPhiU;
  std::vector<double> _gradPhiP;
  std::vector<double> _uDotGradu;
  std::vector<double> _uDotGradPhiu;

public:
  feSysElm_GLS_NavierStokes(feFunction *coeff, feFunction *density, feFunction *viscosity,
                           feVectorFunction *volumeForce)
    : feSysElm(-1, 2, GLS_NAVIERSTOKES, true), _coeff(coeff), _density(density),
      _viscosity(viscosity), _volumeForce(volumeForce)
  {
    _computeMatrixWithFD = true;
  };
  ~feSysElm_GLS_NavierStokes(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_GLS_NavierStokes)
};

#endif
