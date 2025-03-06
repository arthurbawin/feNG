#include "feBilinearForm.h"
#include "feQuadrature.h"

extern int FE_VERBOSE;

feStatus createBilinearForm(feBilinearForm *&form, const std::vector<feSpace *> &spaces,
                            feSysElm *elementarySystem)
{
  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "(BI)LINEAR FORM:");

  if(spaces.size() != (size_t) elementarySystem->getNumFields()) {
    return feErrorMsg(FE_STATUS_ERROR, "(Bi-)linear form %s expects %d FE space(s) (given: %d)\n",
                      elementarySystem->getWeakFormName().data(), elementarySystem->getNumFields(),
                      spaces.size());
  }

  int cncGeoTag = spaces[0]->getCncGeoTag();

  // Perform some checks and exit on error:
  for(feSpace *s : spaces)
  {
    if(s == nullptr) {
      return feErrorMsg(FE_STATUS_ERROR,
                        "Null pointer in vector of FE spaces for (bi)linear form \"%s\", maybe you forgot to initialize it.",
                        elementarySystem->getWeakFormName().data());
    }

    // Check that all spaces are defined (i.e. their DOF were numbered)
    if(!s->wasNumbered()) {
      return feErrorMsg(FE_STATUS_ERROR,
        "(Bi-)linear form %s: FE space %s - %s was defined but its degrees of freedom were not numbered.\n"
        " It was most likely not included in the vector of all FE spaces.\n",
        elementarySystem->getWeakFormName().data(),
        s->getFieldID().data(), s->getCncGeoID().data());
    }

    // Check that all interpolation spaces are defined on the same connectivity
    if(s->getCncGeoTag() != cncGeoTag) {
      return feErrorMsg(
        FE_STATUS_ERROR,
        "(Bi-)linear form %s is defined on more than one geometric connectivity.\n"
        "All finite element spaces in \"spaces\" should be defined on the same geometric "
        "connectivity.\n",
        elementarySystem->getWeakFormName().data());
    }

    // Check that the elementary system is defined on a connectivity of same space dimension
    // Weak forms that can be defined on any dimension have dim = -1
    if(elementarySystem->getDim() != -1 && s->getDim() != elementarySystem->getDim()) {
      return feErrorMsg(FE_STATUS_ERROR,
                        "(Bi-)linear form %s should be defined on connectivity with dimension %d "
                        "(given dimension is %d)\n",
                        elementarySystem->getWeakFormName().data(), elementarySystem->getDim(),
                        s->getDim());
    }
  }

  // TEMPORARY: Check that all FE spaces have the same number of quadrature nodes
  // for(size_t i = 0; i < spaces.size(); ++i) {
  //   if(spaces[i]->getNumQuadPoints() != spaces[0]->getNumQuadPoints()){
  //     return feErrorMsg(FE_STATUS_ERROR, "FE spaces used to define a (bi-)linear form "
  //       "must have the same number of quadrature points\n");
  //   }
  // }

  form = new feBilinearForm(spaces, elementarySystem);

  feInfoCond(FE_VERBOSE > 0, "\t\tDefined (bi-)linear form %s on connectivity %s",
             elementarySystem->getWeakFormName().data(), spaces[0]->getCncGeoID().data());

  return FE_STATUS_OK;
}

static inline double **allocateMatrix(feInt m, feInt n)
{
  double *p = new double[m * n];
  double **A = new double *[m];
  double *q = p;
  for(feInt i = 0; i < m; i++, q += n) A[i] = q;
  return A;
}

static inline void setMatrixToZero(feInt m, feInt n, double **A)
{
  for(feInt i = 0; i < m * n; ++i) (*A)[i] = 0.0;
}

static inline void freeMatrix(feInt m, double **A)
{
  UNUSED(m);
  delete[] A[0];
  delete[] A;
}

static inline void printMatrix(feInt m, feInt n, double ***A)
{
#if defined(HAVE_PETSC)
  for(feInt i = 0; i < m; ++i)
    for(feInt j = 0; j < n; ++j) printf("A[%d][%d] = %+10.16e\n", i, j, (*A)[i][j]);
#else
  for(feInt i = 0; i < m; ++i)
    for(feInt j = 0; j < n; ++j) printf("A[%ld][%ld] = %+10.16e\n", i, j, (*A)[i][j]);
#endif
}

static inline void allocateResidual(feInt m, double **b)
{
  *b = (double *)calloc(m, sizeof *b);
  if(*b == nullptr) {
    printf("In feBilinearForm::allocateResidual : Error - NULL pointer.\n");
  }
}

static inline void setResidualToZero(feInt m, double **b)
{
  for(feInt i = 0; i < m; ++i) (*b)[i] = 0.0;
}

static inline void printResidual(feInt m, double **b)
{
#if defined(HAVE_PETSC)
  for(feInt i = 0; i < m; ++i) printf("b[%d] = %+10.16e\n", i, (*b)[i]);
#else
  for(feInt i = 0; i < m; ++i) printf("b[%ld] = %+10.16e\n", i, (*b)[i]);
#endif
}

static inline void freeResidual(double **b) { free(*b); }

feBilinearForm::feBilinearForm(std::vector<feSpace *> spaces, feSysElm *elementarySystem)
  : _sysElm(elementarySystem), _intSpaces(spaces), _cnc(spaces[0]->getCncGeo()),
    _geoSpace(_cnc->getFeSpace()), _J(_cnc->getJacobians()),
    _cncGeoTag(spaces[0]->getCncGeoTag()), _cncGeoID(spaces[0]->getCncGeoID())
{
  _geoCoord.resize(3 * _cnc->getNumVerticesPerElem());

  // Initialize the elementary system
  _sysElm->createElementarySystem(_intSpaces);
  _fieldsLayoutI = _sysElm->getFieldLayoutI();
  _fieldsLayoutJ = _sysElm->getFieldLayoutJ();

  // Copy the quadrature rule from the first FE space
  // to the elementary system
  // FIXME: they must have the same number of quad points...
  _sysElm->_nQuad = spaces[0]->getNumQuadPoints();
  _sysElm->_wQuad = spaces[0]->getQuadratureWeights();
  _sysElm->_rQuad = spaces[0]->getRQuadraturePoints();
  _sysElm->_sQuad = spaces[0]->getSQuadraturePoints();
  _sysElm->_tQuad = spaces[0]->getTQuadraturePoints();

  // Set the elementary system size (number of DOF in i and j)
  _M = _N = 0;
  for(size_t k = 0; k < _fieldsLayoutI.size(); ++k)
    _M += _intSpaces[_fieldsLayoutI[k]]->getNumFunctions();
  for(size_t k = 0; k < _fieldsLayoutJ.size(); ++k)
    _N += _intSpaces[_fieldsLayoutJ[k]]->getNumFunctions();

  _adrI.resize(_M);
  _adrJ.resize(_N);

  // Allocate memory for the elementary system
  _Ae = allocateMatrix(_M, _N);
  allocateResidual(_M, &_Be);

  // Set the computation method for the element matrix
  // Auxiliary residuals for finite differences computation
  // are allocated anyway, even if the exact Jacobian is computed.
  _h0 = pow(DBL_EPSILON, 1.0 / 2.0);
  allocateResidual(_M, &_R0);
  allocateResidual(_M, &_Rh);
  if(elementarySystem->computeMatrixWithFD()) {
    ptrComputeMatrix = &feBilinearForm::computeMatrixFiniteDifference;
  } else {
    ptrComputeMatrix = &feBilinearForm::computeMatrixAnalytical;
  }

  _adr.resize(_intSpaces.size());
  _sol.resize(_intSpaces.size());
  _solDot.resize(_intSpaces.size());
  _solPrev.resize(_intSpaces.size());
  _solNext.resize(_intSpaces.size());
  for(size_t k = 0; k < _intSpaces.size(); ++k) _adr[k].resize(_intSpaces[k]->getNumFunctions());
  for(size_t k = 0; k < _intSpaces.size(); ++k) _sol[k].resize(_intSpaces[k]->getNumFunctions());
  for(size_t k = 0; k < _intSpaces.size(); ++k) _solDot[k].resize(_intSpaces[k]->getNumFunctions());
  for(size_t k = 0; k < _intSpaces.size(); ++k)
    _solPrev[k].resize(_intSpaces[k]->getNumFunctions());
  for(size_t k = 0; k < _intSpaces.size(); ++k)
    _solNext[k].resize(_intSpaces[k]->getNumFunctions());
}

// Copy constructor
feBilinearForm::feBilinearForm(const feBilinearForm &f)
  : _intSpaces(f._intSpaces), _cnc(f._cnc), _geoSpace(f._geoSpace), _geoCoord(f._geoCoord),
  _J(f._J), _cncGeoTag(f._cncGeoTag), _cncGeoID(f._cncGeoID),
  _fieldsLayoutI(f._fieldsLayoutI), _fieldsLayoutJ(f._fieldsLayoutJ),
  _M(f._M), _N(f._N), _adrI(f._adrI), _adrJ(f._adrJ), _R0(f._R0), _Rh(f._Rh),
  _h0(f._h0), ptrComputeMatrix(f.ptrComputeMatrix), _adr(f._adr), _sol(f._sol),
  _solDot(f._solDot), _solPrev(f._solPrev), _solNext(f._solNext)
{
  _Ae = allocateMatrix(_M, _N);
  allocateResidual(_M, &_Be);
  allocateResidual(_M, &_Rh);
  allocateResidual(_M, &_R0);
  _h0 = pow(DBL_EPSILON, 1.0 / 2.0);

  // Elementary system has to be cloned to be thread-safe
  _sysElm = f._sysElm->clone();
  _sysElm->_nQuad = _intSpaces[0]->getNumQuadPoints();
  _sysElm->_wQuad = _intSpaces[0]->getQuadratureWeights();
  _sysElm->_rQuad = _intSpaces[0]->getRQuadraturePoints();
  _sysElm->_sQuad = _intSpaces[0]->getSQuadraturePoints();
  _sysElm->_tQuad = _intSpaces[0]->getTQuadraturePoints();
}

feBilinearForm::~feBilinearForm()
{
  freeMatrix(_M, _Ae);
  freeResidual(&_Be);
  freeResidual(&_Rh);
  freeResidual(&_R0);
  delete _sysElm;
}

void feBilinearForm::setComputeMatrixWithFD(bool flag)
{
  _sysElm->setComputeMatrixWithFD(flag);
  if(flag) {
    ptrComputeMatrix = &feBilinearForm::computeMatrixFiniteDifference;
  } else {
    ptrComputeMatrix = &feBilinearForm::computeMatrixAnalytical;
  }
}

void feBilinearForm::viewLocalMatrix() { printMatrix(_M, _N, &_Ae); }

void feBilinearForm::viewLocalResidual() { printResidual(_M, &_Be); }

//
// Initialize the connectivity vector for element numElem
// Used to determine the sparsity pattern (mask) of the global matrix
//
void feBilinearForm::initializeAddressingVectors(int numElem)
{
  // feInfo("adr.size = %d on thread %d", _adr.size(), omp_get_thread_num());
  // Initialize _adr for all FE spaces
  for(size_t i = 0; i < _intSpaces.size(); i++) {
    feSpace *fS = _intSpaces[i];
    fS->initializeAddressingVector(numElem, _adr[i]);
    // feInfo("adr[%d].size = %d on thread %d", i, _adr[i].size(), omp_get_thread_num());
  }

  // Initialize continuous addressing vectors in I and J
  // (one for all fields in I, and one for all fields in J)
  for(size_t i = 0, count = 0; i < _fieldsLayoutI.size(); ++i) {
    int M = _intSpaces[_fieldsLayoutI[i]]->getNumFunctions();
    for(int k = 0; k < M; ++k) {
      _adrI[count++] = _adr[_fieldsLayoutI[i]][k];
    };
  }
  for(size_t i = 0, count = 0; i < _fieldsLayoutJ.size(); ++i) {
    int N = _intSpaces[_fieldsLayoutJ[i]]->getNumFunctions();
    for(int k = 0; k < N; ++k) {
      _adrJ[count++] = _adr[_fieldsLayoutJ[i]][k];
    };
  }
}

//
// Initialize the form (feSysElm handler) on element numElem
// Called before computing a local residual or local matrix
//
void feBilinearForm::initialize(feSolution *sol, int numElem)
{
  _numElem = numElem;
  _c0 = sol->getC0();
  _tn = sol->getCurrentTime();
  _dt = sol->getTimeStep();
  _args.t = _tn;

  // Get the coordinates of this element's vertices
  _intSpaces[0]->_mesh->getCoord(_cncGeoTag, numElem, _geoCoord);

  // Get the element transformation once for this element
  // instead of on each quadrature node
  if(_cnc->hasConstantTransformation()) {
    _cnc->getElementTransformation(numElem, _transformation);
  }

  std::vector<double> &solArray = sol->getSolutionReference();

  for(size_t i = 0; i < _intSpaces.size(); i++)
  {
    feSpace *fS = _intSpaces[i];

    // //
    // // Initialize solution on neighbouring elements (for e.g. DG fluxes)
    // // FIXME: Only relevant for 1D meshes
    // //
    // if(numElem > 0) {
    //   fS->initializeAddressingVector(numElem - 1, _adr[i]);
    // } else {
    //   // Periodicity:
    //   fS->initializeAddressingVector(_cnc->getNumElements() - 1, _adr[i]);
    //   // 0 flux
    //   // fS->initializeAddressingVector(numElem, _adr[i]);
    // }
    // for(size_t k = 0; k < _adr[i].size(); ++k) {
    //   _solPrev[i][k] = solArray[_adr[i][k]];
    // }

    // if(numElem < _cnc->getNumElements() - 1) {
    //   fS->initializeAddressingVector(numElem + 1, _adr[i]);
    // } else {
    //   // Periodicity:
    //   fS->initializeAddressingVector(0, _adr[i]);
    //   // 0 flux
    //   // fS->initializeAddressingVector(numElem, _adr[i]);
    // }
    // for(size_t k = 0; k < _adr[i].size(); ++k) {
    //   _solNext[i][k] = solArray[_adr[i][k]];
    // }
    
    //
    // Initialize solution and its time derivative on given element
    //
    fS->initializeAddressingVector(numElem, _adr[i]);
    for(size_t k = 0; k < _adr[i].size(); ++k) {
      _sol[i][k] = solArray[_adr[i][k]];
      _solDot[i][k] = sol->getSolDotAtDOF(_adr[i][k]);
    }
  }

  //
  // Initialize the local connectivity ("address" vector)
  //
  for(size_t i = 0, count = 0; i < _fieldsLayoutI.size(); ++i) {
    int nielm = _intSpaces[_fieldsLayoutI[i]]->getNumFunctions();
    for(int k = 0; k < nielm; ++k) {
      _adrI[count++] = _adr[_fieldsLayoutI[i]][k];
    };
  }
  for(size_t i = 0, count = 0; i < _fieldsLayoutJ.size(); ++i) {
    int nielm = _intSpaces[_fieldsLayoutJ[i]]->getNumFunctions();
    for(int k = 0; k < nielm; ++k) {
      _adrJ[count++] = _adr[_fieldsLayoutJ[i]][k];
    };
  }
}

void feBilinearForm::computeMatrix(feSolution *sol, int numElem)
{
  (this->*feBilinearForm::ptrComputeMatrix)(sol, numElem);
}

void feBilinearForm::computeResidual(feSolution *sol, int numElem)
{
  this->initialize(sol, numElem);
  setResidualToZero(_M, &_Be);
  _sysElm->computeBe(this);
}

void feBilinearForm::computeMatrixAnalytical(feSolution *sol, int numElem)
{
  this->initialize(sol, numElem);
  setMatrixToZero(_M, _N, _Ae);
  _sysElm->computeAe(this);
}

void feBilinearForm::computeMatrixFiniteDifference(feSolution *sol, int numElem)
{
  this->initialize(sol, numElem);
  setMatrixToZero(_M, _N, _Ae);
  setResidualToZero(_M, &_Be);

  // Reference residual (non-perturbated)
  // Compute _Be and copy it in _R0
  _sysElm->computeBe(this);
  for(feInt i = 0; i < _M; i++) _R0[i] = _Be[i];

  // Apply perturbation to _sol and _solDot and recompute each residual
  feInt numColumn = 0;
  for(size_t k = 0; k < _fieldsLayoutJ.size(); k++) {
    feSpace *Unknowns = _intSpaces[_fieldsLayoutJ[k]];

    for(feInt j = 0; j < Unknowns->getNumFunctions(); j++, numColumn++) {
      double temp_sol = _sol[_fieldsLayoutJ[k]][j];
      double temp_soldot = _solDot[_fieldsLayoutJ[k]][j];
      double delta_h = _h0 * std::max(fabs(temp_sol), 1.0);
      double invdelta_h = 1.0 / delta_h;

      _sol[_fieldsLayoutJ[k]][j] = _sol[_fieldsLayoutJ[k]][j] + delta_h;
      _solDot[_fieldsLayoutJ[k]][j] = _solDot[_fieldsLayoutJ[k]][j] + delta_h * _c0;

      setResidualToZero(_M, &_Be);

      // Compute _Be and copy it in _Rh
      _sysElm->computeBe(this);
      for(feInt i = 0; i < _M; i++) _Rh[i] = _Be[i];

      for(feInt i = 0; i < _M; i++) {
        // Matrix and residual are of opposite sign
        _Ae[i][numColumn] = -(_Rh[i] - _R0[i]) * invdelta_h;
      }

      _sol[_fieldsLayoutJ[k]][j] = temp_sol;
      _solDot[_fieldsLayoutJ[k]][j] = temp_soldot;
    }
  }
}

double feBilinearForm::compareAnalyticalAndFDMatrices(feSolution *sol, int numElem)
{
  double **Aexact = allocateMatrix(_M, _N);
  double **Afd    = allocateMatrix(_M, _N);

  computeMatrixAnalytical(sol, numElem);

  for(feInt i = 0; i < _M; ++i) {
    for(feInt j = 0; j < _N; ++j) {
      Aexact[i][j] = _Ae[i][j];
    }
  }
  // feInfo("Exact Jacobian matrix:");
  // printMatrix(_M, _N, &_Ae);

  computeMatrixFiniteDifference(sol, numElem);

  double maxErrorRel = 0.;
  for(feInt i = 0; i < _M; ++i) {
    for(feInt j = 0; j < _N; ++j) {
      Afd[i][j] = _Ae[i][j];
      maxErrorRel = fmax(maxErrorRel, fabs(Aexact[i][j] - Afd[i][j]) / fmax(1e-16, fabs(Afd[i][j])));
    }
  }
  // feInfo("FD Jacobian matrix:");
  // printMatrix(_M, _N, &_Ae);

  for(feInt i = 0; i < _M; ++i) {
    for(feInt j = 0; j < _N; ++j) {
      double abs_error = fabs(Aexact[i][j] - Afd[i][j]);
      double rel_error = fabs(Aexact[i][j] - Afd[i][j]) / fmax(1e-16, fabs(Afd[i][j]));
#if defined(HAVE_PETSC)
      printf("Aex[%2d][%2d] = %+10.16e \t Afd[%2d][%2d] = %+10.16e \t abs error = %+-1.6e \t rel error = %+-1.6e\n",
#else
      printf("Aex[%2ld][%2ld] = %+10.16e \t Afd[%2ld][%2ld] = %+10.16e \t abs error = %+-1.6e \t rel error = %+-1.6e\n",
#endif
        i, j, Aexact[i][j], i, j, Afd[i][j], abs_error, rel_error);
    }
  }

  feInfo("max |Aexact - Afd|/|Afd| = %1.6e", maxErrorRel);

  freeMatrix(_M, Afd);
  freeMatrix(_M, Aexact);
  return maxErrorRel;
}

double feBilinearForm::getMatrixNorm()
{
  int res = 0.0;
  for(feInt i = 0; i < _M; ++i)
    for(feInt j = 0; j < _N; ++j) res += _Ae[i][j] * _Ae[i][j];
  return sqrt(res);
}

double feBilinearForm::getResidualNorm()
{
  int res = 0.0;
  for(feInt i = 0; i < _M; ++i) res += _Be[i] * _Be[i];
  return sqrt(res);
}