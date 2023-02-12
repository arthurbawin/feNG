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
feStatus createBilinearForm(feBilinearForm *&form,
  const std::vector<feSpace *> &spaces, feSysElm *elementarySystem);

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
  // Since friendship is not inherited by feSysElm derived class,
  // they are public for now, even if it's probably not optimal.

  // Finite element spaces used to discretize each field present
  // in the weak form associated to the chosen feSysElm.
  // For example, for a diffusion weak form, there is only a single
  // field e.g. U, hence _intSpaces[0] = U;
  // For the Navier-Stokes weak forms, we have
  // _intSpaces[0] = U;
  // _intSpaces[1] = V;
  // _intSpaces[2] = P;
  //
  // The fields order must match the one given in feSysElm.
  std::vector<feSpace *> _intSpaces;
  // Ptr to the geometric connectivity to access jacobians
  feCncGeo *_cnc;
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

  double _dxdr[3];
  double _dxds[3];
  double _dxdt[3];
  double _drdx[3];
  double _drdy[3];
  double _drdz[3];

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
  //         U V P -> unknown fields
  // phi_U [       ]    
  // phi_V [       ]    _fieldsLayoutI = {0 1 2}
  // phi_P [       ]    _fieldsLayoutJ = {0 1 2}
  //   |
  //   |------------> test functions
  //
  std::vector<int> _fieldsLayoutI;
  std::vector<int> _fieldsLayoutJ;

  // Dimensions of the element-wise linear system
  feInt _M;
  feInt _N;

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
  void (feBilinearForm::*ptrComputeMatrix)(feSolution *sol,
                                           int numElem);

public:
  // Addressing vector of the current element for each FE space
  std::vector<std::vector<feInt>> _adr;

  // Solution at DOFs on the current element for each FE space
  std::vector<std::vector<double>> _sol;
  std::vector<std::vector<double>> _solDot;

  // Solution on the previous and next element (for e.g. DG fluxes)
  std::vector<std::vector<double>> _solPrev;
  std::vector<std::vector<double>> _solNext;

public:
  // Create a (bi-)linear form to compute the element-wise weak form
  // defined by elementarySystem with interpolation and test functions
  // defined in spaces/vectorSpaces.
  feBilinearForm(std::vector<feSpace*> spaces, feSysElm *elementarySystem);

  feBilinearForm(const feBilinearForm &f);
  ~feBilinearForm();

  feCncGeo *getCncGeo() { return _cnc; }
  int getCncGeoTag() { return _cncGeoTag; }
  feInt getLocalMatrixM() { return _M; }
  feInt getLocalMatrixN() { return _N; }
  std::vector<feInt> &getAdrI() { return _adrI; }
  std::vector<feInt> &getAdrJ() { return _adrJ; }
  double **getAe() { return _Ae; }
  double *getBe() { return _Be; }
  elementSystemType getID() { return _sysElm->getID(); }
  std::string getWeakFormName() { return _sysElm->getWeakFormName(); }

  // Return true if there is a local matrix associated to the weak form
  // (false if there is only a residual).
  bool hasMatrix() { return _sysElm->hasMatrix(); }

  // Compute element-wise matrix and residual on element numElem
  void computeMatrix(feSolution *sol, int numElem);
  void computeResidual(feSolution *sol, int numElem);

  // Initialize the addressing vectors _adrI and _adrJ on element numElem
  void initializeAddressingVectors(int numElem);

  double getMatrixNorm();
  double getResidualNorm();

private:
  // Initialize the form on element numElem.
  // Initializes _adr and the local solutions _sol, _solDot, _solPrev, _solNext
  void initialize(feSolution *sol, int numElem);

  void computeMatrixAnalytical(feSolution *sol,
                               int numElem);
  void computeMatrixFiniteDifference(feSolution *sol,
                                     int numElem);
};

#endif