#include "feBilinearForm.h"
#include "feQuadrature.h"

extern int FE_VERBOSE;

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

feStatus createBilinearForm(feBilinearForm *&form, std::vector<feSpace *> spaces,
  feMesh *mesh, int degQuad, feSysElm *elementarySystem)
{
  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "(BI-)LINEAR FORM:");

  // Perform some checks and exit on error
  // Check that all interpolation spaces are defined on the same connectivity
  int cncGeoTag = spaces[0]->getCncGeoTag();
  for(size_t i = 0; i < spaces.size(); ++i) {
    if(spaces[i]->getCncGeoTag() != cncGeoTag)
      return feErrorMsg(FE_STATUS_ERROR, "(Bi-)linear form %s is defined on more than one geometric connectivity.\n"
        "All finite element spaces in \"spaces\" should be defined on the same geometric connectivity.\n",
        elementarySystem->getIDName().data());
  }
  // Check that the elementary system is defined on a connectivity of same space dimension
  for(size_t i = 0; i < spaces.size(); ++i) {
    if(spaces[i]->getDim() != elementarySystem->getDim()){
      return feErrorMsg(FE_STATUS_ERROR, "(Bi-)linear form %s should be defined on connectivity with dimension %d "
        "(given dimension is %d)\n",
        elementarySystem->getIDName().data(),
        elementarySystem->getDim(),
        spaces[i]->getDim());
    }
  }

  form = new feBilinearForm(spaces, mesh, degQuad, elementarySystem);

  feInfoCond(FE_VERBOSE > 0, "\t\tDefined (bi-)linear form %s on connectivity %s",
    elementarySystem->getIDName().data(),
    spaces[0]->getCncGeoID().data());

  return FE_STATUS_OK;
}

/* Une feBilinearForm:
  Stocke la matrice et le residu elementaire
  Gère l'adressage dans la matrice globale
  Appelle le feSysElm pour calculer matrice et residu
  (elementary system is stored here but computed in feSysElm)
*/
feBilinearForm::feBilinearForm(std::vector<feSpace *> spaces, 
  feMesh *mesh, int degQuad, feSysElm *elementarySystem)
  : _sysElm(elementarySystem)
  , _intSpace(spaces)
  , _cnc(spaces[0]->getCncGeo())
  , _cncGeoID(spaces[0]->getCncGeoID())
  , _cncGeoTag(spaces[0]->getCncGeoTag())
  , _geoSpace(_cnc->getFeSpace())
  , _nGeoElm(_cnc->getNbElm())
  , _degQuad(degQuad)
{
  _nGeoNodes = _cnc->getNbNodePerElem();
  _geoCoord.resize(3*_nGeoNodes);

  // Initialize the elementary system computer
  _sysElm->createElementarySystem(_intSpace);
  _iVar = _sysElm->getIVar();
  _jVar = _sysElm->getJVar();

  // Set the elementary system size (number of DOF in i and j)
  _niElm = _njElm = 0;
  for(size_t k = 0; k < _iVar.size(); ++k)
    _niElm += _intSpace[_iVar[k]]->getNbFunctions();
  for(size_t k = 0; k < _jVar.size(); ++k)
    _njElm += _intSpace[_jVar[k]]->getNbFunctions();

  _adrI.resize(_niElm);
  _adrJ.resize(_njElm);

  // Allocate memory for the elementary system
  _Ae = allocateMatrix(_niElm, _njElm);
  allocateResidual(_niElm, &_Be);

  // Set the computation method for the element matrix
  if(elementarySystem->computeMatrixWithFD()) {
    ptrComputeMatrix = &feBilinearForm::computeMatrixFiniteDifference;
      _R0 = nullptr;
      _Rh = nullptr;
      allocateResidual(_niElm, &_Rh);
      allocateResidual(_niElm, &_R0);
      _h0 = pow(DBL_EPSILON, 1.0 / 2.0);
  } else {
    ptrComputeMatrix = &feBilinearForm::computeMatrixAnalytical;
  }

  _adr.resize(_intSpace.size());
  _sol.resize(_intSpace.size());
  _solDot.resize(_intSpace.size());
  _solPrev.resize(_intSpace.size());
  _solNext.resize(_intSpace.size());
  for(size_t k = 0; k < _intSpace.size(); ++k)
    _adr[k].resize(_intSpace[k]->getNbFunctions());
  for(size_t k = 0; k < _intSpace.size(); ++k)
    _sol[k].resize(_intSpace[k]->getNbFunctions());
  for(size_t k = 0; k < _intSpace.size(); ++k)
    _solDot[k].resize(_intSpace[k]->getNbFunctions());
  for(size_t k = 0; k < _intSpace.size(); ++k)
    _solPrev[k].resize(_intSpace[k]->getNbFunctions());
  for(size_t k = 0; k < _intSpace.size(); ++k)
    _solNext[k].resize(_intSpace[k]->getNbFunctions());
}

// Copy constructor
feBilinearForm::feBilinearForm(const feBilinearForm &f)
  : _intSpace(f._intSpace)
  , _cnc(f._cnc)
  , _cncGeoID(f._cncGeoID)
  , _cncGeoTag(f._cncGeoTag)
  , _geoSpace(f._geoSpace)
  , _nGeoElm(f._nGeoElm)
  , _degQuad(f._degQuad)
  , _nGeoNodes(f._nGeoNodes)
  , _iVar(f._iVar)
  , _jVar(f._jVar)
  , _niElm(f._niElm)
  , _njElm(f._njElm)
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
  _Ae = allocateMatrix(_niElm, _njElm);
  allocateResidual(_niElm, &_Be);
  if(f._sysElm->computeMatrixWithFD()) {
    _R0 = nullptr;
    _Rh = nullptr;
    allocateResidual(_niElm, &_Rh);
    allocateResidual(_niElm, &_R0);
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
        f._sysElm->getIDName().data());
      exit(-1);
  }
}

feBilinearForm::~feBilinearForm()
{
  freeMatrix(_niElm, _Ae);
  freeResidual(&_Be);
  if(_sysElm->computeMatrixWithFD()) {
    freeResidual(&_Rh);
    freeResidual(&_R0);
  }
  delete _sysElm;
}

void feBilinearForm::initialize_vadij_only(feMetaNumber *metaNumber, int numElem)
{
  for(size_t i = 0; i<_intSpace.size(); i++) {
    feSpace *fS = _intSpace[i];
    fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), numElem, _adr[i]);
  }

  // INITIALISERVADIJ
  for(size_t i = 0, count = 0; i < _iVar.size(); ++i) {
    int nielm = _intSpace[_iVar[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k) {
      _adrI[count++] = _adr[i][k];
    };
  }
  for(size_t i = 0, count = 0; i < _jVar.size(); ++i) {
    int nielm = _intSpace[_jVar[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k) {
      _adrJ[count++] = _adr[i][k];
    };
  }
}

void feBilinearForm::initialize(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol, int numElem)
{
  for(size_t i = 0; i<_intSpace.size(); i++) {
    feSpace *fS = _intSpace[i];
    std::vector<double> &solVec = sol->getSolutionReference();

    // Initialize solution on neighbouring elements (for e.g. DG fluxes)
    if(numElem > 0){
      fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), numElem - 1, _adr[i]);
    } else{
      // Periodicity:
      fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), _cnc->getNbElm() - 1, _adr[i]);
    }

    for(size_t k = 0; k < _adr[i].size(); ++k) {
      _solPrev[i][k] = solVec[_adr[i][k]];
    }

    if(numElem < _cnc->getNbElm() - 1){
      fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), numElem + 1, _adr[i]);
    } else{
      // Periodicity:
      fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), 0, _adr[i]);
    }

    for(size_t k = 0; k < _adr[i].size(); ++k) {
      _solNext[i][k] = solVec[_adr[i][k]];
    }

    // Initialize solution and its time derivative on current element
    fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), numElem, _adr[i]);
    for(size_t k = 0; k < _adr[i].size(); ++k) {
      _sol[i][k] = solVec[_adr[i][k]];
      _solDot[i][k] = sol->getSolDotAtDOF(_adr[i][k]);
    }
  }

  mesh->getCoord(_cncGeoTag, numElem, _geoCoord);

  // INITIALISERVADIJ
  for(size_t i = 0, count = 0; i < _iVar.size(); ++i) {
    int nielm = _intSpace[_iVar[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k) {
      _adrI[count++] = _adr[_iVar[i]][k];
    };
  }
  for(size_t i = 0, count = 0; i < _jVar.size(); ++i) {
    int nielm = _intSpace[_jVar[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k) {
      _adrJ[count++] = _adr[_jVar[i]][k];
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
  setResidualToZero(_niElm, &_Be);
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
  setMatrixToZero(_niElm, _njElm, _Ae);
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
  setMatrixToZero(_niElm, _njElm, _Ae);

  _numElem = numElem;
  _c0 = sol->getC0();
  _tn = sol->getCurrentTime();
  _dt = sol->getTimeStep();

  // ==================================================================
  // Le résidu non perturbé
  // ==================================================================
  setResidualToZero(_niElm, &_Be);

  // Compute _Be and copy it in _R0
  _sysElm->computeBe(this);
  for(feInt i = 0; i < _niElm; i++)
    _R0[i] = _Be[i];
  
  // ==================================================================
  // Calcul des résidus perturbés
  // ==================================================================
  feInt numColumn = 0;
  for(size_t k = 0; k < _jVar.size(); k++) {

    feSpace *Unknowns = _intSpace[_jVar[k]];

    for(feInt j = 0; j < Unknowns->getNbFunctions(); j++, numColumn++) {
      double temp_sol = _sol[_jVar[k]][j]; 
      double temp_soldot = _solDot[_jVar[k]][j];
      double delta_h = _h0 * std::max(fabs(temp_sol), 1.0);
      double invdelta_h = 1.0 / delta_h;

      _sol[_jVar[k]][j]    = _sol[_jVar[k]][j] + delta_h;
      _solDot[_jVar[k]][j] = _solDot[_jVar[k]][j] + delta_h * _c0;

      setResidualToZero(_niElm, &_Be);

      // Compute _Be and copy it in _Rh
      _sysElm->computeBe(this);
      for(feInt i = 0; i < _niElm; i++)
        _Rh[i] = _Be[i];

      for(feInt i = 0; i < _niElm; i++)
        _Ae[i][numColumn] = (-_Rh[i] + _R0[i]) * invdelta_h;

      _sol[_jVar[k]][j]    = temp_sol;
      _solDot[_jVar[k]][j] = temp_soldot;
    }
  }
}

double feBilinearForm::getMatrixNorm()
{
  int res = 0.0;
  for(feInt i = 0; i < _niElm; ++i)
    for(feInt j = 0; j < _njElm; ++j) res += _Ae[i][j] * _Ae[i][j];
  return sqrt(res);
}

double feBilinearForm::getResidualNorm()
{
  int res = 0.0;
  for(feInt i = 0; i < _niElm; ++i) res += _Be[i] * _Be[i];
  return sqrt(res);
}