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
  MASS_POWER,
  MIXED_MASS,
  MIXED_MASS_POWER,
  VECTOR_MASS,
  MIXED_VECTOR_MASS,
  MIXED_SCALAR_VECTOR_MASS,
  TRANSIENT_MASS,
  MIXED_TRANSIENT_MASS,
  TRANSIENT_VECTOR_MASS,

  DIFFUSION,
  VECTOR_DIFFUSION,
  NONLINEAR_DIFFUSION,

  ADVECTION,
  NONLINEAR_ADVECTION,
  VECTOR_CONVECTIVE_ACCELERATION,
  VECTOR_ADJOINT_CONVECTIVE_ACCELERATION,
  TRACER_CONVECTION,

  DIV_NEWTONIAN_STRESS,

  MIXED_GRADIENT,
  MIXED_GRADIENT_FIELD_DEPENDENT_COEFF,
  MIXED_GRADGRAD,
  TRIPLE_MIXED_GRADIENT,

  MIXED_DIVERGENCE,
  MIXED_CURL,
  MIXED_DOTPRODUCT,

  CHNS_MOMENTUM,
  MIXED_DIVERGENCE_CHNS,

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
    case MASS_POWER:
      return "MASS_POWER";
    case MIXED_MASS:
      return "MIXED_MASS";
    case MIXED_MASS_POWER:
      return "MIXED_MASS_POWER";
    case VECTOR_MASS:
      return "VECTOR_MASS";
    case MIXED_VECTOR_MASS:
      return "MIXED_VECTOR_MASS";
    case MIXED_SCALAR_VECTOR_MASS:
      return "MIXED_SCALAR_VECTOR_MASS";
    case TRANSIENT_MASS:
      return "TRANSIENT_MASS";
    case MIXED_TRANSIENT_MASS:
      return "MIXED_TRANSIENT_MASS";
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
    case TRACER_CONVECTION:
      return "TRACER_CONVECTION";
    case DIV_NEWTONIAN_STRESS:
      return "DIV_NEWTONIAN_STRESS";
    case MIXED_GRADIENT:
      return "MIXED_GRADIENT";
    case MIXED_GRADIENT_FIELD_DEPENDENT_COEFF:
      return "MIXED_GRADIENT_FIELD_DEPENDENT_COEFF";
    case MIXED_GRADGRAD:
      return "MIXED_GRADGRAD";
    case MIXED_DIVERGENCE:
      return "MIXED_DIVERGENCE";
    case MIXED_CURL:
      return "MIXED_CURL";
    case MIXED_DOTPRODUCT:
      return "MIXED_DOTPRODUCT";
    case TRIPLE_MIXED_GRADIENT:
      return "TRIPLE_MIXED_GRADIENT";
    case CHNS_MOMENTUM:
      return "CHNS_MOMENTUM";
    case MIXED_DIVERGENCE_CHNS:
      return "MIXED_DIVERGENCE_CHNS";
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

  // Number of vector components on the FE space
  int _nComponents;

  // Compute the element FE matrix using finite differences?
  bool _computeMatrixWithFD = false;
  // Does this form have an elementary matrix to assemble?
  bool _hasMatrix = false;
  // Is this form associated to a d/dt operator?
  bool _isTransientMatrix = false;
  // Is the computation thread safe?
  bool _isThreadSafe = true;

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
    : _ID(ID), _dim(dim), _nFields(nFields), _hasMatrix(hasMatrix)
  {
    _fieldsLayoutI.resize(1);
    _fieldsLayoutJ.resize(1);
  };
  virtual ~feSysElm() {}

  virtual feSysElm *clone() const = 0;

  int getDim() const { return _dim; }
  int getNumFields() const { return _nFields; }
  const std::vector<int> &getFieldLayoutI() const { return _fieldsLayoutI; }
  const std::vector<int> &getFieldLayoutJ() const { return _fieldsLayoutJ; }
  elementSystemType getID() const { return _ID; }
  std::string getWeakFormName() const { return toString(_ID); }
  bool computeMatrixWithFD() const { return _computeMatrixWithFD; }
  void setComputeMatrixWithFD(bool flag) { _computeMatrixWithFD = flag; }
  bool hasMatrix() const { return _hasMatrix; }
  bool isTransientMatrix() const { return _isTransientMatrix; }
  bool isThreadSafe() const { return _isThreadSafe; }

  virtual void createElementarySystem(std::vector<feSpace *> &spaces) { UNUSED(spaces); };
  virtual void computeAe(feBilinearForm *form) { UNUSED(form); };
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
  int _idU, _nFunctions;
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
  int _idU, _nFunctions;
  std::vector<double> _phiAtX0;

public:
  feSysElm_SourceDirac(std::vector<double> &x0) : feSysElm(-1, 1, SOURCE, false), _x0(x0)
  {
    _isThreadSafe = false;
  };
  ~feSysElm_SourceDirac(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_SourceDirac)
};

class feSysElm_VectorSource : public feSysElm
{
protected:
  feVectorFunction *_source;
  int _idU, _nFunctions;
  std::vector<double> _S, _SdotPhi;

public:
  feSysElm_VectorSource(feVectorFunction *source)
    : feSysElm(-1, 1, VECTOR_SOURCE, false), _source(source), _S(3, 0.){};
  ~feSysElm_VectorSource(){};
  void createElementarySystem(std::vector<feSpace *> &spaces);
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
  int _idU, _nFunctions;
  std::vector<double> _phiU;

public:
  feSysElm_Mass(feFunction *coeff) : feSysElm(-1, 1, MASS, true), _coeff(coeff)
  {
    _computeMatrixWithFD = false;
  };
  ~feSysElm_Mass(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_Mass)
};

//
// Mass matrix with field to a given power p.
// Matrix and residual
//
//  /          p
//  | coeff * u * v dx
//  /
//
// where v is the test function of u.
//
// # fields: 1 (FE solution and test functions)
//                        U
// Fields layout: phi_U [   ]
//
class feSysElm_MassPower : public feSysElm
{
protected:
  feFunction *_coeff;
  double _p;
  int _idU, _nFunctions;
  std::vector<double> _phiU;

public:
  feSysElm_MassPower(feFunction *coeff, double exponent)
    : feSysElm(-1, 1, MASS, true), _coeff(coeff), _p(exponent)
  {
    _computeMatrixWithFD = true;
  };
  ~feSysElm_MassPower(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_MassPower)
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
  feSysElm_MixedMass(feFunction *coeff) : feSysElm(-1, 2, MIXED_MASS, true), _coeff(coeff)
  {
    _computeMatrixWithFD = false;
  };
  ~feSysElm_MixedMass(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_MixedMass)
};

//
// Same as MIXED_MASS, but instead computes the integral
//
//  /          p
//  | coeff * u * v dx, with p a real number exponent.
//  /
//
class feSysElm_MixedMassPower : public feSysElm
{
protected:
  feFunction *_coeff;
  double _p;
  int _idU, _idV, _nFunctionsU, _nFunctionsV;
  std::vector<double> _phiU, _phiV;

public:
  feSysElm_MixedMassPower(feFunction *coeff, int exponent)
    : feSysElm(-1, 2, MIXED_MASS_POWER, true), _coeff(coeff), _p(exponent)
  {
    _computeMatrixWithFD = false;
  };
  ~feSysElm_MixedMassPower(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_MixedMassPower)
};

class feSysElm_VectorMass : public feSysElm
{
protected:
  feFunction *_coeff;
  int _idU, _nFunctions;
  std::vector<double> _u, _phi_idotphi_j, _udotphi_i;

public:
  feSysElm_VectorMass(feFunction *coeff) : feSysElm(-1, 1, VECTOR_MASS, true), _coeff(coeff)
  {
    _computeMatrixWithFD = false;
  };
  ~feSysElm_VectorMass(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_VectorMass)
};

//
// Mass matrix for one vector-valued variable tested with the test functions
// of another vector-valued variable
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
  int _idU, _idV, _nFunctionsU, _nFunctionsV, _nComponentsU;
  std::vector<double> _udotphi_i, _phi_idotphi_j;

public:
  feSysElm_MixedVectorMass(feFunction *coeff)
    : feSysElm(-1, 2, MIXED_VECTOR_MASS, true), _coeff(coeff){};
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
  feSysElm_MixedScalarVectorMass(feFunction *coeff)
    : feSysElm(-1, 2, MIXED_SCALAR_VECTOR_MASS, true), _coeff(coeff)
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
  feFunction *_coeff;
  int _idU, _nFunctions;
  std::vector<double> _feU;

public:
  feSysElm_TransientMass(feFunction *coeff) : feSysElm(-1, 1, TRANSIENT_MASS, true), _coeff(coeff)
  {
    _computeMatrixWithFD = false;
    _isTransientMatrix = true;
  };
  ~feSysElm_TransientMass(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_TransientMass)
};

//
// Transient term weak form (transient mass matrix) tested with the test functions of another scalar
// variable. Matrix and residual
//
//  /
//  | coeff * dudt * phi_v dx
//  /
//
// # fields: 2 (test functions and field whose time derivative is computed)
//                        U
// Fields layout: phi_V [   ]
//
class feSysElm_MixedTransientMass : public feSysElm
{
protected:
  feFunction *_coeff;
  int _idU, _idV, _nFunctionsU, _nFunctionsV;
  std::vector<double> _phiU, _phiV;

public:
  feSysElm_MixedTransientMass(feFunction *coeff)
    : feSysElm(-1, 2, MIXED_TRANSIENT_MASS, true), _coeff(coeff)
  {
    _computeMatrixWithFD = true;
    _isTransientMatrix = true;
  };
  ~feSysElm_MixedTransientMass(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_MixedTransientMass)
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
  feFunction *_coeff;
  int _idU, _nFunctions;
  std::vector<double> _dudt, _dudtdotphi_i, _phi_idotphi_j;

public:
  feSysElm_TransientVectorMass(feFunction *coeff)
    : feSysElm(-1, 1, TRANSIENT_VECTOR_MASS, true), _coeff(coeff)
  {
    _computeMatrixWithFD = false;
    _isTransientMatrix = true;
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
// Strong form: - div(coeff * grad(u))
//
//            /
// Weak form: | coeff * grad(u) dot grad(v) dx
//            /
//
// # fields: 1 (FE solution and test functions)
//
// Parameters:
//  coeff: scalar field returning the diffusivity k(x)
//
//                        U
// Fields layout: phi_U [   ]
//
template <int dim> class feSysElm_Diffusion : public feSysElm
{
protected:
  feFunction *_coeff;
  int _idU, _nFunctions;
  std::vector<double> _gradPhi;

public:
  feSysElm_Diffusion(feFunction *coeff) : feSysElm(dim, 1, DIFFUSION, true), _coeff(coeff){};
  ~feSysElm_Diffusion(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_Diffusion)
};

//
// Diffusion of vector-valued field
// Matrix and residual
//
// Strong form: - div(coeff * grad(u)), u vector
//
//             /
// Weak form:  | coeff * grad(u) : grad(v) dx
//             /
//
class feSysElm_VectorDiffusion : public feSysElm
{
protected:
  feFunction *_coeff;
  feFunction *_diffusivity;
  int _idU, _nFunctions;
  std::vector<double> _gradu, _gradPhi, _doubleContraction, _doubleContraction_u;

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
  feFunction *_diffusivity, *_ddiffdu;
  int _idU, _nFunctions;
  std::vector<double> _phiU, _gradPhi;

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
  int _idU, _nFunctions;
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
  int _idU, _nFunctions;
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
// This is the directional derivative along u:
//
//     coeff * (u_i * d/dx_i)u_j * v_j
//   = coeff *  u_i *  du_j/dx_i * v_j
//
// Note that if the velocity gradient is defined as (grad u)_ij = du_i/dx_j,
// then what is computed is (in cartesian coordinates):
//
//             coeff * (u dot (grad u)^T) dot v
// instead of
//             coeff * (u dot (grad u)  ) dot v,
//
// which is associated to (grad u)_ij = du_j/dx_i
//
class feSysElm_VectorConvectiveAcceleration : public feSysElm
{
protected:
  feFunction *_coeff;
  int _idU, _nFunctions;
  std::vector<double> _u, _gradu, _uDotGradu, _phiU, _gradPhiU, _uDotGraduDotPhiU,
    _u0DotGradPhiUDotPhiU, _phiUDotGradu0DotPhiU;

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
// Mixed weak form for: u dot grad c (c scalar).
//
// Integral of coeff * (u dot grad c) dot phi_c,
// where:
//  - u is the solved vector-valued velocity field,
//  - c is a different solved scalar field (active or passive tracer),
//  - phi_c are the test functions of c.
//
class feSysElm_TracerConvection : public feSysElm
{
protected:
  feFunction *_coeff;
  int _idU, _idC;
  int _nFunctionsU, _nFunctionsC;
  std::vector<double> _u, _gradC, _phiC, _gradPhiC, _phiUdotGradc0;

public:
  feSysElm_TracerConvection(feFunction *coeff)
    : feSysElm(-1, 2, TRACER_CONVECTION, true), _coeff(coeff)
  {
    _computeMatrixWithFD = false;
  };
  ~feSysElm_TracerConvection(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_TracerConvection)
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
// Hence: integral of +p * div(v) - mu * d : grad(v) dx
// Residual is the negative of this integral.
//
class feSysElm_DivergenceNewtonianStress : public feSysElm
{
protected:
  feFunction *_coeff;
  feFunction *_viscosity;
  int _idU, _idP, _nFunctionsU, _nFunctionsP;
  std::vector<double> _gradu, _symmetricGradu, _gradPhiU, _phiP, _divPhiU, _doubleContraction,
    _doubleContractionPhiPhi, _doubleContractionPhiPhiT;

public:
  feSysElm_DivergenceNewtonianStress(feFunction *coeff, feFunction *viscosity)
    : feSysElm(-1, 2, DIV_NEWTONIAN_STRESS, true), _coeff(coeff), _viscosity(viscosity)
  {
    _computeMatrixWithFD = false;
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
  int _idU, _nFunctions;
  std::vector<double> _gradPhi, _divPhi;

public:
  feSysElm_GradSource(feFunction *source) : feSysElm(-1, 1, GRAD_SOURCE, false), _source(source){};
  ~feSysElm_GradSource(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_GradSource)
};

//
// Mixed gradient : (coeff * grad(u)) cdot v, with v: vector-valued test functions
//                                                 coeff: a scalar coefficient
// Integral of - coeff * u * div(v)
//
class feSysElm_MixedGradient : public feSysElm
{
protected:
  feFunction *_coeff;
  int _idU, _idV, _nFunctionsU, _nFunctionsV;
  std::vector<double> _phiU, _gradPhiV, _divPhiV;

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
// Mixed gradient with a coefficient depending on another field :
//
// Strong form : coeff(w) * grad(p) cdot v
//
// with : w a scalar-valued variable
//        coeff a scalar coefficient depending on w
//        p a scalar-valued variable
//        v the test functions of a vector-valued variable
//
// Weak forn : integral of - coeff(w) * p * div(v)
//
class feSysElm_MixedGradientFieldDependentCoeff : public feSysElm
{
protected:
  feFunction *_coeff;
  int _idV, _idP, _idW, _nFunctionsV, _nFunctionsP, _nFunctionsW;
  std::vector<double> _phiV, _phiP;
  std::vector<double> _gradPhiV;

public:
  feSysElm_MixedGradientFieldDependentCoeff(feFunction *coeff)
    : feSysElm(-1, 3, MIXED_GRADIENT_FIELD_DEPENDENT_COEFF, true), _coeff(coeff)
  {
    _computeMatrixWithFD = true;
  };
  ~feSysElm_MixedGradientFieldDependentCoeff(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_MixedGradientFieldDependentCoeff)
};

//
// Mixed gradient with test functions of a different variable.
//
// Strong form:       - div (coeff * grad(u))
//   Weak form: integral of (coeff * grad(u)) dot grad(phi_v)
//
// with:
//  - coeff: a scalar coefficient
//  -     u: a resolved SCALAR field
//  - phi_v: the test functions of another resolved SCALAR field
//
template <int dim> class feSysElm_MixedGradGrad : public feSysElm
{
protected:
  feFunction *_coeff;
  int _idU, _idV, _nFunctionsU, _nFunctionsV;
  std::vector<double> _gradPhiU, _gradPhiV;

public:
  feSysElm_MixedGradGrad(feFunction *coeff) : feSysElm(dim, 2, MIXED_GRADGRAD, true), _coeff(coeff)
  {
    _computeMatrixWithFD = false;
  };
  ~feSysElm_MixedGradGrad(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_MixedGradGrad)
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
  int _idU, _idV, _nFunctionsU, _nFunctionsV;
  std::vector<double> _gradu, _gradPhiU, _divPhiU, _phiV;

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

//
// Mixed divergence with density function of the phase marker phi:
//
//    coeff * (rho(phi) * div(u)) * q, with      q (= phiP) : pressure test functions

// Integral of coeff * (div(rho(phi) * u)) * v
//
class feSysElm_MixedDivergenceCHNS : public feSysElm
{
protected:
  feFunction *_coeff;
  feFunction *_density;
  feFunction *_drhodphi;
  int _idP, _idU, _idPhi;
  int _nFunctionsP, _nFunctionsU;
  std::vector<double> _u, _gradu, _gradphi, _gradPhiU, _phiU, _phiP;

public:
  feSysElm_MixedDivergenceCHNS(feFunction *coeff, feFunction *density, feFunction *drhodphi)
    : feSysElm(-1, 3, MIXED_DIVERGENCE_CHNS, true), _coeff(coeff), _density(density),
      _drhodphi(drhodphi)
  {
    _computeMatrixWithFD = false;
  };
  ~feSysElm_MixedDivergenceCHNS(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_MixedDivergenceCHNS)
};

class feSysElm_MixedCurl : public feSysElm
{
protected:
  feFunction *_coeff;
  int _idU, _idV, _nFunctionsU, _nFunctionsV;
  std::vector<double> _gradu, _gradPhiU, _phiV;

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

// Mixed dot product : u cdot f, with u an unknown vector field and
//                                    f a user-defined vector field.
// Integral of (u cdot f) dot v, with v the test function of another scalar unknown field.
class feSysElm_MixedDotProduct : public feSysElm
{
protected:
  feVectorFunction *_f;
  int _idU, _idV, _nFunctionsU, _nFunctionsV;
  std::vector<double> _u, _fVal, _phiU, _phiV;

public:
  feSysElm_MixedDotProduct(feVectorFunction *f) : feSysElm(2, 2, MIXED_DOTPRODUCT, true), _f(f)
  {
    _computeMatrixWithFD = false;
  };
  ~feSysElm_MixedDotProduct(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_MixedDotProduct)
};

//
// Mixed gradient with 3 fields.
//
// Strong form:              coeff * u * grad(v)
//   Weak form: integral of (coeff * u * grad(v)) cdot phi_w
//
// with:
//  - coeff: a scalar coefficient
//  -     u: a resolved SCALAR field
//  -     v: another resolved SCALAR field
//  - phi_w: the test functions of yet another resolved VECTOR field
//
class feSysElm_TripleMixedGradient : public feSysElm
{
protected:
  feFunction *_coeff;
  int _idU, _idV, _idW;
  int _nFunctionsU, _nFunctionsV, _nFunctionsW;
  std::vector<double> _gradV;
  std::vector<double> _phiW;

public:
  feSysElm_TripleMixedGradient(feFunction *coeff)
    : feSysElm(-1, 3, TRIPLE_MIXED_GRADIENT, true), _coeff(coeff)
  {
    _computeMatrixWithFD = true;
  };
  ~feSysElm_TripleMixedGradient(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_TripleMixedGradient)
};

//
// Momentum equation for the Cahn-Hilliard Navier-Stokes model:
//
// Strong form:
//  rho(phi) * (du/dt + (u dot grad) u - f) - coeff * mu * grad(phi) - div(sigma)
//
// Weak form:
//  integral of rho(phi) * (du/dt + (u dot grad) u - f) dot phi_u + sigma : grad(phi_u)
//
// with:
//  -   rho: the density which depends on phi
//  -   phi: the phase marker
//  -     u: the velocity field
//  - coeff: a coefficient, e.g. gamma/eps (surface tension/interface thickness)
//  -    mu: the chemical potential
//  - sigma: the Newtonian stress tensor whose viscosity also depends on phi
//  -     f: volume source term
//  - phi_u: the velocity test functions
//
// Requires 4 fields: u, p (through sigma), phi, mu.
// The density is a function depending on the solution.
//
class feSysElm_CHNS_Momentum : public feSysElm
{
protected:
  feFunction *_density;
  feFunction *_drhodphi;
  feFunction *_viscosity;
  feFunction *_dviscdphi;
  feFunction *_mobility;
  feFunction *_coeffKorteweg;
  feVectorFunction *_volumeForce;
  int _idU, _idP, _idPhi, _idMu;
  int _nFunctionsU, _nFunctionsP, _nFunctionsPhi, _nFunctionsMu;
  // A lot of vectors for the dot products and contractions
  std::vector<double> _f, _u, _dudt, _gradu, _uDotGradu, _gradphi, _gradmu, _gradmuDotGradu,
    _symmetricGradu, _phiU, _gradPhiU, _phiP, _phiPhi, _gradPhiPhi, _phiMu, _gradPhiMu,
    _dudtDotPhiU, _uDotGraduDotPhiU, _fDotPhiU, _gradPhiDotphiU, _gradMuDotgradUdotphiU, _divPhiU,
    _doubleContraction, _phi_idotphi_j, _u0DotGradPhiUDotPhiU, _phiUDotGradu0DotPhiU,
    _doubleContractionPhiPhi, _doubleContractionPhiPhiT, _gradMu0DotgradUdotphiU,
    _gradMuDotgradU0DotphiU, _gradPhi0dotPhiU, _gradPhiPhiDotPhiU, _symGraduDDotGradPhiU;

public:
  feSysElm_CHNS_Momentum(feFunction *density, feFunction *drhodphi, feFunction *viscosity,
                         feFunction *dviscdphi, feFunction *mobility, feFunction *coeffKorteweg,
                         feVectorFunction *volumeForce)
    : feSysElm(-1, 4, CHNS_MOMENTUM, true), _density(density), _drhodphi(drhodphi),
      _viscosity(viscosity), _dviscdphi(dviscdphi), _mobility(mobility),
      _coeffKorteweg(coeffKorteweg), _volumeForce(volumeForce)
  {
    _computeMatrixWithFD = false;
  };
  ~feSysElm_CHNS_Momentum(){};
  void createElementarySystem(std::vector<feSpace *> &space);
  void computeAe(feBilinearForm *form);
  void computeBe(feBilinearForm *form);
  CLONEABLE(feSysElm_CHNS_Momentum)
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
  int _idU, _nFunctions;
  std::vector<double> _feU;

public:
  feSysElm_1D_NeumannBC(feFunction *neumannBC)
    : feSysElm(1, 1, NEUMANN_1D, false), _neumannBC(neumannBC){};
  ~feSysElm_1D_NeumannBC() {}
  void createElementarySystem(std::vector<feSpace *> &space);
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

  int _idU, _nFunctions;
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
  int _idU, _nFunctions;
  std::vector<double> _phiU;
  std::vector<double> _gradPhiU;

public:
  feSysElm_2D_SUPGStab(feFunction *reactionCoeff, feFunction *diffusivity,
                       feVectorFunction *velocity, feFunction *source)
    : feSysElm(2, 1, SUPG_STABILIZATION_2D, true), _velocity(velocity), _diffusivity(diffusivity),
      _reactionCoeff(reactionCoeff), _source(source)
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
  int _idU, _idP, _nFunctionsU, _nFunctionsP;
  std::vector<double> _f, _residual, _u, _gradu, _hessu, _gradp, _gradPhiU, _gradPhiP, _hessPhiU,
    _uDotGradu, _uDotGradPhiu;

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
  int _idU, _idP, _nFunctionsU, _nFunctionsP;
  std::vector<double> _f, _residual, _u, _gradu, _hessu, _gradp, _gradPhiU, _hessPhiU, _gradPhiP;

public:
  feSysElm_GLS_Stokes(feFunction *coeff, feFunction *density, feFunction *viscosity,
                      feVectorFunction *volumeForce)
    : feSysElm(-1, 2, GLS_STOKES, true), _coeff(coeff), _density(density), _viscosity(viscosity),
      _volumeForce(volumeForce)
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
  int _idU, _idP, _nFunctionsU, _nFunctionsP;
  std::vector<double> _f, _residual, _u, _gradu, _hessu, _gradp, _gradPhiU, _hessPhiU, _gradPhiP,
    _uDotGradu, _uDotGradPhiu;

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
