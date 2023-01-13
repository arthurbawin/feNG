#include "feBilinearForm.h"
#include "feQuadrature.h"

extern int FE_VERBOSE;

feStatus createBilinearForm(feBilinearForm *&form, const std::vector<feSpace *> &spaces, feSysElm *elementarySystem)
{
  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "(BI-)LINEAR FORM:");

  // Perform some checks and exit on error:

  // Check that all interpolation spaces are defined on the same connectivity
  int cncGeoTag = spaces[0]->getCncGeoTag();
  for(size_t i = 0; i < spaces.size(); ++i) {
    if(spaces[i]->getCncGeoTag() != cncGeoTag)
      return feErrorMsg(FE_STATUS_ERROR, "(Bi-)linear form %s is defined on more than one geometric connectivity.\n"
        "All finite element spaces in \"spaces\" should be defined on the same geometric connectivity.\n",
        elementarySystem->getWeakFormName().data());
  }

  // Check that the elementary system is defined on a connectivity of same space dimension
  for(size_t i = 0; i < spaces.size(); ++i) {
    if(spaces[i]->getDim() != elementarySystem->getDim()){
      return feErrorMsg(FE_STATUS_ERROR, "(Bi-)linear form %s should be defined on connectivity with dimension %d "
        "(given dimension is %d)\n",
        elementarySystem->getWeakFormName().data(),
        elementarySystem->getDim(),
        spaces[i]->getDim());
    }
  }

  // Check that the number of FE spaces match the required 
  // number of spaces to compute the weak form
  if(spaces.size() != elementarySystem->getNumFields()){
    return feErrorMsg(FE_STATUS_ERROR, "(Bi-)linear form %s expects %d FE space(s) (given: %d)\n",
        elementarySystem->getWeakFormName().data(),
        elementarySystem->getNumFields(),
        spaces.size());
  }

  // TEMPORARY: Check that all FE spaces have the same number of quadrature nodes
  for(size_t i = 0; i < spaces.size(); ++i) {
    if(spaces[i]->getNbQuadPoints() != spaces[0]->getNbQuadPoints()){
      return feErrorMsg(FE_STATUS_ERROR, "FE spaces used to define a (bi-)linear form "
        "must have the same number of quadrature points\n");
    }
  }


  form = new feBilinearForm(spaces, elementarySystem);

  feInfoCond(FE_VERBOSE > 0, "\t\tDefined (bi-)linear form %s on connectivity %s",
    elementarySystem->getWeakFormName().data(),
    spaces[0]->getCncGeoID().data());

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
  delete[] A[0];
  delete[] A;
}

static inline void printMatrix(feInt m, feInt n, double ***A)
{
  for(feInt i = 0; i < m; ++i)
    for(feInt j = 0; j < n; ++j) printf("A[%d][%d] = %+10.16e\n", i, j, (*A)[i][j]);
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
  for(feInt i = 0; i < m; ++i) printf("b[%d] = %+10.16e\n", i, (*b)[i]);
}

static inline void freeResidual(double **b) { free(*b); }

feBilinearForm::feBilinearForm(std::vector<feSpace *> spaces, feSysElm *elementarySystem)
  : _sysElm(elementarySystem)
  , _intSpaces(spaces)
  , _cnc(spaces[0]->getCncGeo())
  , _cncGeoID(spaces[0]->getCncGeoID())
  , _cncGeoTag(spaces[0]->getCncGeoTag())
  , _geoSpace(_cnc->getFeSpace())
{
  _geoCoord.resize(3 * _cnc->getNbNodePerElem());

  _J = _cnc->getJacobians();

  // Initialize the elementary system
  _sysElm->createElementarySystem(_intSpaces);
  _fieldsLayoutI = _sysElm->getFieldLayoutI();
  _fieldsLayoutJ = _sysElm->getFieldLayoutJ();

  // Copy the quadrature rule from the first FE space
  // to the elementary system
  // FIXME: they must have the same number of quad points...
  _sysElm->_nQuad = spaces[0]->getNbQuadPoints();
  _sysElm->_wQuad = spaces[0]->getQuadratureWeights(); 
  _sysElm->_rQuad = spaces[0]->getRQuadraturePoints();
  _sysElm->_sQuad = spaces[0]->getSQuadraturePoints();
  _sysElm->_tQuad = spaces[0]->getTQuadraturePoints();

  // Set the elementary system size (number of DOF in i and j)
  _M = _N = 0;
  for(size_t k = 0; k < _fieldsLayoutI.size(); ++k)
    _M += _intSpaces[_fieldsLayoutI[k]]->getNbFunctions();
  for(size_t k = 0; k < _fieldsLayoutJ.size(); ++k)
    _N += _intSpaces[_fieldsLayoutJ[k]]->getNbFunctions();

  _adrI.resize(_M);
  _adrJ.resize(_N);

  // Allocate memory for the elementary system
  _Ae = allocateMatrix(_M, _N);
  allocateResidual(_M, &_Be);

  // Set the computation method for the element matrix
  if(elementarySystem->computeMatrixWithFD()) {
    ptrComputeMatrix = &feBilinearForm::computeMatrixFiniteDifference;
      _R0 = nullptr;
      _Rh = nullptr;
      allocateResidual(_M, &_Rh);
      allocateResidual(_M, &_R0);
      _h0 = pow(DBL_EPSILON, 1.0 / 2.0);
  } else {
    ptrComputeMatrix = &feBilinearForm::computeMatrixAnalytical;
  }

  _adr.resize(_intSpaces.size());
  _sol.resize(_intSpaces.size());
  _solDot.resize(_intSpaces.size());
  _solPrev.resize(_intSpaces.size());
  _solNext.resize(_intSpaces.size());
  for(size_t k = 0; k < _intSpaces.size(); ++k)
    _adr[k].resize(_intSpaces[k]->getNbFunctions());
  for(size_t k = 0; k < _intSpaces.size(); ++k)
    _sol[k].resize(_intSpaces[k]->getNbFunctions());
  for(size_t k = 0; k < _intSpaces.size(); ++k)
    _solDot[k].resize(_intSpaces[k]->getNbFunctions());
  for(size_t k = 0; k < _intSpaces.size(); ++k)
    _solPrev[k].resize(_intSpaces[k]->getNbFunctions());
  for(size_t k = 0; k < _intSpaces.size(); ++k)
    _solNext[k].resize(_intSpaces[k]->getNbFunctions());
}

// Copy constructor
feBilinearForm::feBilinearForm(const feBilinearForm &f)
  : _intSpaces(f._intSpaces)
  , _cnc(f._cnc)
  , _cncGeoID(f._cncGeoID)
  , _cncGeoTag(f._cncGeoTag)
  , _geoSpace(f._geoSpace)
  , _fieldsLayoutI(f._fieldsLayoutI)
  , _fieldsLayoutJ(f._fieldsLayoutJ)
  , _M(f._M)
  , _N(f._N)
  , _adr(f._adr)
  , _adrI(f._adrI)
  , _adrJ(f._adrJ)
  , _sol(f._sol)
  , _solDot(f._solDot)
  , _solPrev(f._solPrev)
  , _solNext(f._solNext)
  , _R0(f._R0)
  , _Rh(f._Rh)
  , _h0(f._h0)
  , ptrComputeMatrix(f.ptrComputeMatrix)
  , _geoCoord(f._geoCoord)
{
  _Ae = allocateMatrix(_M, _N);
  allocateResidual(_M, &_Be);
  if(f._sysElm->computeMatrixWithFD()) {
    _R0 = nullptr;
    _Rh = nullptr;
    allocateResidual(_M, &_Rh);
    allocateResidual(_M, &_R0);
    _h0 = pow(DBL_EPSILON, 1.0 / 2.0);
  }

  switch(f._sysElm->getID()) {
    // case STIFFSPRING_0D :
    //   _sysElm=new feSysElm_0D_StiffSpring(f._sysElm);
    //   break;
    // case STIFF2_0D :
    //   _sysElm=new feSysElm_0D_Stiff2();
    //   break;
    // case STIFF3_0D :
    //   _sysElm=new feSysElm_0D_Stiff3();
    //   break;
    // case WEAKBC_0D :
    //   _sysElm=new feSysElm_0D_weakBC();
    //   break;
    // case WEAKBC_EDO1_0D :
    //   _sysElm=new feSysElm_0D_weakBC_edo1();
    //   break;
    // case WEAKBC_EDO1_V2_0D :
    //   _sysElm=new feSysElm_0D_weakBC_edo1_V2();
    //   break;
    // case WEAKBC_EDO2_0D :
    //   _sysElm=new feSysElm_0D_weakBC_edo2();
    //   break;
    // case MASSE_0D :
    //   _sysElm=new feSysElm_0D_Masse();
    //   break;
    // case SOURCE_0D :
    //   _sysElm=new feSysElm_0D_Source();
    //   break;
    // case SOURCE_CROSSED_0D :
    //   _sysElm=new feSysElm_0D_Source_crossed();
    //   break;
    // case WEAKBC_EDO1_1D :
    //   _sysElm=new feSysElm_1D_weakBC_edo1();
    //   break;
    // case SOURCE_1D :
    //   _sysElm=new feSysElm_1D_Source();
    //   break;
    // case DIFF_1D :
    //   _sysElm=new feSysElm_1D_Diffusion();
    //   break;
    // case MASSE_1D :
    //   _sysElm=new feSysElm_1D_Masse(();
    //     break;
    case MASS_2D :
      _sysElm = new feSysElm_2D_Masse(static_cast<feSysElm_2D_Masse &>(*f._sysElm));
      break;
    case SOURCE_2D:
      _sysElm = new feSysElm_2D_Source(static_cast<feSysElm_2D_Source &>(*f._sysElm));
      break;
    case DIFFUSION_2D:
      _sysElm = new feSysElm_2D_Diffusion(static_cast<feSysElm_2D_Diffusion &>(*f._sysElm));
      break;
    case NEUMANN_1D:
      _sysElm = new feSysElm_1D_NeumannBC(static_cast<feSysElm_1D_NeumannBC &>(*f._sysElm));
      break;
    case ADVECTION_2D:
      _sysElm = new feSysElm_2D_Advection(static_cast<feSysElm_2D_Advection &>(*f._sysElm));
      break;
    case STOKES_2D :
      _sysElm = new feSysElm_2D_Stokes(static_cast<feSysElm_2D_Stokes &>(*f._sysElm));
      break;
    case NAVIERSTOKES_2D:
      _sysElm = new feSysElm_2D_NavierStokes(static_cast<feSysElm_2D_NavierStokes &>(*f._sysElm));
      break;
    default:
      feWarning("Copy constructor for feBilinearForm - Elementary system \"%s\" is unknown. Exiting.",
        f._sysElm->getWeakFormName().data());
      exit(-1);
  }
}

feBilinearForm::~feBilinearForm()
{
  freeMatrix(_M, _Ae);
  freeResidual(&_Be);
  if(_sysElm->computeMatrixWithFD()) {
    freeResidual(&_Rh);
    freeResidual(&_R0);
  }
  delete _sysElm;
}

void feBilinearForm::initializeAddressingVectors(feMetaNumber *metaNumber, int numElem)
{
  // Initialize _adr for all FE spaces
  for(size_t i = 0; i < _intSpaces.size(); i++) {
    feSpace *fS = _intSpaces[i];
    fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), numElem, _adr[i]);
  }

  // Initialize continuous addressing vectors in I and J
  // (one for all fields in I, and one for all fields in J)
  for(size_t i = 0, count = 0; i < _fieldsLayoutI.size(); ++i) {
    int M = _intSpaces[_fieldsLayoutI[i]]->getNbFunctions();
    for(int k = 0; k < M; ++k) {
      _adrI[count++] = _adr[i][k];
    };
  }
  for(size_t i = 0, count = 0; i < _fieldsLayoutJ.size(); ++i) {
    int N = _intSpaces[_fieldsLayoutJ[i]]->getNbFunctions();
    for(int k = 0; k < N; ++k) {
      _adrJ[count++] = _adr[i][k];
    };
  }
}

void feBilinearForm::initialize(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol, int numElem)
{
  for(size_t i = 0; i<_intSpaces.size(); i++) {
    feSpace *fS = _intSpaces[i];
    std::vector<double> &solRef = sol->getSolutionReference();

    // Initialize solution on neighbouring elements (for e.g. DG fluxes)
    // if(numElem > 0){
    //   fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), numElem - 1, _adr[i]);
    // } else{
    //   // Periodicity:
    //   fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), _cnc->getNbElm() - 1, _adr[i]);
    // }

    // for(size_t k = 0; k < _adr[i].size(); ++k) {
    //   _solPrev[i][k] = solRef[_adr[i][k]];
    // }

    // if(numElem < _cnc->getNbElm() - 1){
    //   fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), numElem + 1, _adr[i]);
    // } else{
    //   // Periodicity:
    //   fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), 0, _adr[i]);
    // }

    // for(size_t k = 0; k < _adr[i].size(); ++k) {
    //   _solNext[i][k] = solRef[_adr[i][k]];
    // }

    // Initialize solution and its time derivative on current element
    fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), numElem, _adr[i]);
    for(size_t k = 0; k < _adr[i].size(); ++k) {
      _sol[i][k] = solRef[_adr[i][k]];
      _solDot[i][k] = sol->getSolDotAtDOF(_adr[i][k]);
    }
  }

  mesh->getCoord(_cncGeoTag, numElem, _geoCoord);

  // INITIALISERVADIJ
  for(size_t i = 0, count = 0; i < _fieldsLayoutI.size(); ++i) {
    int nielm = _intSpaces[_fieldsLayoutI[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k) {
      _adrI[count++] = _adr[_fieldsLayoutI[i]][k];
    };
  }
  for(size_t i = 0, count = 0; i < _fieldsLayoutJ.size(); ++i) {
    int nielm = _intSpaces[_fieldsLayoutJ[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k) {
      _adrJ[count++] = _adr[_fieldsLayoutJ[i]][k];
    };
  }
}

void feBilinearForm::computeMatrix(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol,
                                   int numElem)
{
  (this->*feBilinearForm::ptrComputeMatrix)(metaNumber, mesh, sol, numElem);
}

void feBilinearForm::computeResidual(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol,
                                     int numElem)
{
  this->initialize(metaNumber, mesh, sol, numElem);
  setResidualToZero(_M, &_Be);
  _numElem = numElem;
  _c0 = sol->getC0();
  _tn = sol->getCurrentTime();
  _dt = sol->getTimeStep();
  _sysElm->computeBe(this);
}

void feBilinearForm::computeMatrixAnalytical(feMetaNumber *metaNumber, feMesh *mesh,
                                             feSolution *sol, int numElem)
{
  this->initialize(metaNumber, mesh, sol, numElem);
  setMatrixToZero(_M, _N, _Ae);
  _numElem = numElem;
  _c0 = sol->getC0();
  _tn = sol->getCurrentTime();
  _dt = sol->getTimeStep();

  _sysElm->computeAe(this);
}

void feBilinearForm::computeMatrixFiniteDifference(feMetaNumber *metaNumber, feMesh *mesh,
                                                   feSolution *sol, int numElem)
{
  this->initialize(metaNumber, mesh, sol, numElem);
  setMatrixToZero(_M, _N, _Ae);

  _numElem = numElem;
  _c0 = sol->getC0();
  _tn = sol->getCurrentTime();
  _dt = sol->getTimeStep();

  // ==================================================================
  // Le résidu non perturbé
  // ==================================================================
  setResidualToZero(_M, &_Be);

  // Compute _Be and copy it in _R0
  _sysElm->computeBe(this);
  for(feInt i = 0; i < _M; i++)
    _R0[i] = _Be[i];
  
  // ==================================================================
  // Calcul des résidus perturbés
  // ==================================================================
  feInt numColumn = 0;
  for(size_t k = 0; k < _fieldsLayoutJ.size(); k++) {

    feSpace *Unknowns = _intSpaces[_fieldsLayoutJ[k]];

    for(feInt j = 0; j < Unknowns->getNbFunctions(); j++, numColumn++) {
      double temp_sol = _sol[_fieldsLayoutJ[k]][j]; 
      double temp_soldot = _solDot[_fieldsLayoutJ[k]][j];
      double delta_h = _h0 * std::max(fabs(temp_sol), 1.0);
      double invdelta_h = 1.0 / delta_h;

      _sol[_fieldsLayoutJ[k]][j]    = _sol[_fieldsLayoutJ[k]][j] + delta_h;
      _solDot[_fieldsLayoutJ[k]][j] = _solDot[_fieldsLayoutJ[k]][j] + delta_h * _c0;

      setResidualToZero(_M, &_Be);

      // Compute _Be and copy it in _Rh
      _sysElm->computeBe(this);
      for(feInt i = 0; i < _M; i++)
        _Rh[i] = _Be[i];

      for(feInt i = 0; i < _M; i++)
        _Ae[i][numColumn] = (-_Rh[i] + _R0[i]) * invdelta_h;

      _sol[_fieldsLayoutJ[k]][j]    = temp_sol;
      _solDot[_fieldsLayoutJ[k]][j] = temp_soldot;
    }
  }
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