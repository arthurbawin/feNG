#ifndef _FEBILINEARFORM_
#define _FEBILINEARFORM_

#include "feNG.h"
#include "feMesh.h"
#include "feSpace.h"
#include "feSysElm.h"
#include "feQuadrature.h"
#include "feNumber.h"
#include "feSolution.h"
#include "feCncGeo.h"

// Create a linear or bilinear form and perform safety checks.
// This is the recommended way of creating a linear form.
// Call within feCheck( ... ) to safely exit if a problem arises.
//
//                   form: pointer to the form, initially undefined and assigned
//                         during the call.
//                 spaces: vector of FE spaces (interpolation functions) used to
//                         compute the form
//       elementarySystem: the elementary system (discretized weak form) to be
//                         solved on each element
feStatus createBilinearForm(feBilinearForm *&form, const std::vector<feSpace *> &spaces,
                            feSysElm *elementarySystem);

//
// A linear or bilinear form.
//
// Stores the elementary matrix and residual and handles
// the assembly from local to global matrix. Initializes all
// necessary structures before calling the feSysElm class,
// which performs the actual computation of the weak forms.
//
class feBilinearForm
{
protected:
  // Elementary system to be solved on each element.
  // Compute the element matrix and residual.
  feSysElm *_sysElm;

public:
  // These members are used to compute the elementary system.
  // Since friendship is not inherited by feSysElm derived classes,
  // they are public for now, even if it's probably not ideal.

  // Finite element spaces used to discretize each field present
  // in the weak form associated to the chosen feSysElm.
  // For example, for a diffusion weak form, there is only a single
  // field e.g. U, hence _intSpaces[0] = U;
  // For the Navier-Stokes weak forms, we have
  // _intSpaces[0] = U; (vector FE space)
  // _intSpaces[1] = P; (scalar FE space)
  //
  // The fields order must match the one given in feSysElm.
  std::vector<feSpace *> _intSpaces;
  // Ptr to the geometric connectivity (can probably be removed)
  const feCncGeo *_cnc;
  // Ptr to the FE space used to interpolate the geometry
  feSpace *_geoSpace;
  // Physical coordinates of a point in the reference space
  std::vector<double> _geoCoord;

  // Elements jacobians
  const std::vector<double> &_J;

  // First coefficient of the BDF expansion
  double _c0;
  // Current time
  double _tn;
  // Time step
  double _dt;

  // The local index of the current element in the geometric connectivity
  int _numElem;

  // The local matrix and residual
  double **_Ae;
  double *_Be;

  // Jacobian matrix and inverse of reference-to-physical element transformation
  ElementTransformation _transformation;

  // The arguments to give the various callbacks
  feFunctionArguments _args;

protected:
  // Tag and name of the geometric connectivity on which the form is defined
  int _cncGeoTag;
  std::string _cncGeoID;

  // Fields layout of the weak form in I,J in the local matrix.
  // It is a trivial vector _fieldsLayoutI[i] = i for most weak forms,
  // since we test the fields with their own basis functions.
  //
  // Example for Navier-Stokes equation: local matrix has the form:
  //
  //         U P -> unknown fields
  // phi_U [       ] _fieldsLayoutI = {0 1}
  // phi_P [       ] _fieldsLayoutJ = {0 1}
  //   |
  //   |------------> test functions
  //
  // For mixed weak forms, however, the local matrix is rectangular
  // and the layouts differ.
  std::vector<int> _fieldsLayoutI;
  std::vector<int> _fieldsLayoutJ;

public:
  // Dimensions of the element-wise linear system
  feInt _M;
  feInt _N;

protected:
  // Addressing vector (local to global mapping) in I and J
  // Continuous vector that spans all FE spaces, such that
  // the local matrix is written at (_adrI[i], _adrJ[j]) in
  // the global matrix.
  std::vector<feInt> _adrI;
  std::vector<feInt> _adrJ;

  // Residuals to compute the elementary matrix with finite differences
  double *_R0;
  double *_Rh;
  double _h0;
  // Function pointer to the matrix computation method
  void (feBilinearForm::*ptrComputeMatrix)(const feSolution *sol, const int numElem);

public:
  // Addressing vector of the current element for each FE space
  std::vector<std::vector<feInt> > _adr;

  // Solution at DOFs on the current element for each FE space
  std::vector<std::vector<double>> _sol;
  std::vector<std::vector<double>> _solDot;

  // Solution on the previous and next element (for e.g. DG fluxes)
  // Next and previous only make sense in 1D for now.
  std::vector<std::vector<double>> _solPrev;
  std::vector<std::vector<double>> _solNext;

public:
  // Create a (bi-)linear form to compute the element-wise weak form
  // defined by elementarySystem with interpolation and test functions
  // defined in spaces/vectorSpaces.
  feBilinearForm(const std::vector<feSpace*> spaces,
                 feSysElm *elementarySystem);
  feBilinearForm(const feBilinearForm &f);
  ~feBilinearForm();

  const feCncGeo *getCncGeo() const { return _cnc; }
  int getCncGeoTag() const { return _cncGeoTag; }
  feInt getLocalMatrixM() const { return _M; }
  feInt getLocalMatrixN() const { return _N; }
  const std::vector<feInt> &getAdrI() const { return _adrI; }
  const std::vector<feInt> &getAdrJ() const { return _adrJ; }
  const double* const *getAe() const { return _Ae; }
  const double* getBe() const { return _Be; }
  elementSystemType getID() const { return _sysElm->getID(); }
  std::string getWeakFormName() const { return _sysElm->getWeakFormName(); }

  // Return true if there is a local matrix associated to the weak form
  // (false if there is only a residual).
  bool hasMatrix() const { return _sysElm->hasMatrix(); }
  bool isTransientMatrix() const { return _sysElm->isTransientMatrix(); }

  // Sets the Jacobian matrix to be evaluated numerically
  // using finite differences, allocates the necessary arrays.
  void setComputeMatrixWithFD(bool flag);
  double compareAnalyticalAndFDMatrices(const feSolution *sol, const int numElem);

  // Compute element-wise matrix and residual on element numElem
  void computeMatrix(const feSolution *sol, const int numElem);
  void computeResidual(const feSolution *sol, const int numElem);

  // Initialize the addressing vectors _adrI and _adrJ on element numElem
  void initializeAddressingVectors(int numElem);

  double getMatrixNorm() const;
  double getResidualNorm() const;

  void viewLocalMatrix() const;
  void viewLocalResidual() const;

private:
  // Initialize the form on element numElem.
  // Initializes _adr and the local solutions _sol, _solDot, _solPrev, _solNext
  void initialize(const feSolution *sol, const int numElem);

  void computeMatrixAnalytical(const feSolution *sol, const int numElem);
  void computeMatrixFiniteDifference(const feSolution *sol, const int numElem);
};

#endif