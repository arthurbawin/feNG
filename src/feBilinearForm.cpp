#include "feBilinearForm.h"
#include "feQuadrature.h"

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
    for(feInt j = 0; j < n; ++j) printf("A[%ld][%ld] = %+10.16e\n", i, j, (*A)[i][j]);
  // std::cout<<"A["<<i<<"]["<<j<<"] = "<<(*A)[i][j]<<std::endl;
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
  for(feInt i = 0; i < m; ++i) printf("b[%ld] = %+10.16e\n", i, (*b)[i]);
}

static inline void freeResidual(double **b) { free(*b); }

feBilinearForm::feBilinearForm(std::vector<feSpace *> space, feMesh *mesh, int degQuad,
                               feSysElm *sysElm)
  : _sysElm(sysElm), _intSpace(space), _cnc(space[0]->getCncGeo()),
    _cncGeoID(space[0]->getCncGeoID()), _cncGeoTag(space[0]->getCncGeoTag()),
    _geoSpace(_cnc->getFeSpace()), _nGeoElm(_cnc->getNbElm()), _degQuad(degQuad)

{
  _nCoord = mesh->getDim();
  _nGeoNodes = _cnc->getNbNodePerElem();
  _geoCoord.resize(3*_nGeoNodes);

  if(_geoSpace->getCncGeoTag() != _cncGeoTag)
    printf("In feBilinearForm::feBilinearForm : Erreur - Connectivité géométrique ne correspond "
           "pas à la connectivité de l'espace d'interpolation.\n");

  for(size_t i = 0; i < _intSpace.size(); ++i) {
    if(_intSpace[i]->getCncGeoTag() != _cncGeoTag)
      printf("In feBilinearForm::feBilinearForm : Erreur - Forme définie sur plus d'une "
             "connectivité.\n");
  }
  _sysElm->createElementarySystem(_intSpace);
  
  _iVar = _sysElm->getIVar();
  _jVar = _sysElm->getJVar();

  // CREERVADIJ
  _niElm = _njElm = 0;
  for(size_t k = 0; k < _iVar.size(); ++k) _niElm += _intSpace[_iVar[k]]->getNbFunctions();
  for(size_t k = 0; k < _jVar.size(); ++k) _njElm += _intSpace[_jVar[k]]->getNbFunctions();

  // _adr.resize(_intSpace[0]->getNbFunctions());
  // _sol.resize(_intSpace[0]->getNbFunctions());
  // _solDot.resize(_intSpace[0]->getNbFunctions());

  _adr.resize(_intSpace.size());
  _sol.resize(_intSpace.size());
  _solDot.resize(_intSpace.size());
  for(size_t k = 0; k < _intSpace.size(); ++k) _adr[k].resize(_intSpace[k]->getNbFunctions());
  for(size_t k = 0; k < _intSpace.size(); ++k) _sol[k].resize(_intSpace[k]->getNbFunctions());
  for(size_t k = 0; k < _intSpace.size(); ++k) _solDot[k].resize(_intSpace[k]->getNbFunctions());
  
  _adrI.resize(_niElm);
  _adrJ.resize(_njElm);

  _Ae = allocateMatrix(_niElm, _njElm);
  allocateResidual(_niElm, &_Be);

  _R0 = NULL;
  _Rh = NULL;
  allocateResidual(_niElm, &_Rh);
  allocateResidual(_niElm, &_R0);
  _h0 = pow(DBL_EPSILON, 1.0 / 2.0);

  // Set the computation method for the element matrix
  if(sysElm->getMatrixAnalyticalStatus()) {
    ptrComputeMatrix = &feBilinearForm::computeMatrixAnalytical;
  } else {
    ptrComputeMatrix = &feBilinearForm::computeMatrixFiniteDifference;
  }
}

feBilinearForm::feBilinearForm(const feBilinearForm &f)
  : _intSpace(f._intSpace), _cnc(f._cnc), _cncGeoID(f._cncGeoID), _cncGeoTag(f._cncGeoTag),
    _geoSpace(f._geoSpace), _nGeoElm(f._nGeoElm), _degQuad(f._degQuad), _nCoord(f._nCoord),
    _nGeoNodes(f._nGeoNodes), _iVar(f._iVar), _jVar(f._jVar), _niElm(f._niElm), _njElm(f._njElm),
    _adr(f._adr), _adrI(f._adrI), _adrJ(f._adrJ), _sol(f._sol), _solDot(f._solDot), _R0(f._R0),
    _Rh(f._Rh), _h0(f._h0), ptrComputeMatrix(f.ptrComputeMatrix), _geoCoord(f._geoCoord)
{
  _Ae = allocateMatrix(_niElm, _njElm);
  allocateResidual(_niElm, &_Be);

  _R0 = NULL;
  _Rh = NULL;
  allocateResidual(_niElm, &_Rh);
  allocateResidual(_niElm, &_R0);
  _h0 = pow(DBL_EPSILON, 1.0 / 2.0);

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
    case WEAKBC_0D :
      _sysElm=new feSysElm_0D_weakBC(static_cast<feSysElm_0D_weakBC &>(*f._sysElm));
      break;
    case WEAKBC_EDO1_0D :
      _sysElm=new feSysElm_0D_weakBC_edo1(static_cast<feSysElm_0D_weakBC_edo1 &>(*f._sysElm));
      break;
    // case WEAKBC_EDO1_V2_0D :
    //   _sysElm=new feSysElm_0D_weakBC_edo1_V2();
    //   break;
    case WEAKBC_EDO2_0D :
      _sysElm=new feSysElm_0D_weakBC_edo2(static_cast<feSysElm_0D_weakBC_edo2 &>(*f._sysElm));
      break;
    case MASSE_0D :
      _sysElm=new feSysElm_0D_Masse(static_cast<feSysElm_0D_Masse &>(*f._sysElm));
      break;
    case SOURCE_0D :
      _sysElm=new feSysElm_0D_Source(static_cast<feSysElm_0D_Source &>(*f._sysElm));
      break;
    // case SOURCE_CROSSED_0D :
    //   _sysElm=new feSysElm_0D_Source_crossed();
    //   break;
    case WEAKBC_EDO1_1D :
      _sysElm=new feSysElm_1D_weakBC_edo1(static_cast<feSysElm_1D_weakBC_edo1 &>(*f._sysElm));
      break;
    case SOURCE_1D :
      _sysElm=new feSysElm_1D_Source(static_cast<feSysElm_1D_Source &>(*f._sysElm));
      break;
    case DIFF_1D :
      _sysElm=new feSysElm_1D_Diffusion(static_cast<feSysElm_1D_Diffusion &>(*f._sysElm));
      break;
    case MASSE_1D :
      _sysElm=new feSysElm_1D_Masse(static_cast<feSysElm_1D_Masse &>(*f._sysElm));
        break;
    case MASSE_2D :
      _sysElm=new feSysElm_2D_Masse(static_cast<feSysElm_2D_Masse &>(*f._sysElm));
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
      _sysElm=new feSysElm_2D_Stokes(static_cast<feSysElm_2D_Stokes &>(*f._sysElm));
      break;
    case NAVIERSTOKES_2D:
      _sysElm = new feSysElm_2D_NavierStokes(static_cast<feSysElm_2D_NavierStokes &>(*f._sysElm));
      break;
    // case DIRECTIONALDONOTHING_1D :
    //   _sysElm=new feSysElm_1D_DirectionalDoNothing();
    //   break;
    default:
      feInfo("No elementary system => default nullptr");
      // _sysElm=new feSysElm_2D_Diffusion(static_cast<feSysElm_2D_Diffusion&>(*f._sysElm));
      _sysElm = nullptr;
      exit(-1);
      break;
  }
}

feBilinearForm::~feBilinearForm()
{
  // freeMatrix(_niElm, &_Ae);
  freeMatrix(_niElm, _Ae);
  freeResidual(&_Be);
  freeResidual(&_Rh);
  freeResidual(&_R0);
  delete _sysElm;
}

void feBilinearForm::initialize_vadij_only(feMetaNumber *metaNumber, int numElem)
{
  for(int i=0; i<_intSpace.size(); i++) {
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

void feBilinearForm::initialize(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol,
                                int numElem)
{
  for(int i=0; i<_intSpace.size(); i++) {
    feSpace *fS = _intSpace[i];
    fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), numElem, _adr[i]);

    std::vector<double> &solVec = sol->getSolutionReference();
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
  // this->initialize(metaNumber, mesh, sol, numElem);
  // setMatrixToZero(_niElm, _njElm, &_Ae);
  // _sysElm->computeAe(_intSpace, _geoSpace, _geoCoord, sol->getC0(), sol->getCurrentTime(), _Ae);
  // printMatrix(_niElm, _njElm, &_Ae);
  (this->*feBilinearForm::ptrComputeMatrix)(metaNumber, mesh, sol, numElem);
}

void feBilinearForm::computeMatrixAnalytical(feMetaNumber *metaNumber, feMesh *mesh,
                                             feSolution *sol, int numElem)
{
  this->initialize(metaNumber, mesh, sol, numElem);
  setMatrixToZero(_niElm, _njElm, _Ae);
  std::vector<double> &J = _cnc->getJacobians();
  _sysElm->computeAe(J, numElem, _intSpace, _geoSpace, _geoCoord, sol->getC0(),
                     sol->getCurrentTime(), _Ae, _sol, _solDot);
  // printMatrix(_niElm, _njElm, &_Ae);
}
// void feBilinearForm::computeMatrixFiniteDifference(feMetaNumber *metaNumber, feMesh *mesh,
// feSolution *sol, int numElem) {
//   this->initialize(metaNumber, mesh, sol, numElem);
//   setMatrixToZero(_niElm, _njElm, &_Ae);
//   _sysElm->computeAe(_intSpace, _geoSpace, _geoCoord, sol->getC0(), sol->getCurrentTime(), _Ae);
// }

void feBilinearForm::computeResidual(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol,
                                     int numElem)
{
  this->initialize(metaNumber, mesh, sol, numElem);
  setResidualToZero(_niElm, &_Be);

  std::vector<double> &J = _cnc->getJacobians();
  _sysElm->computeBe(J, numElem, _intSpace, _geoSpace, _geoCoord, sol->getC0(),
                     sol->getCurrentTime(), sol->getTimeStep(), _Be, _sol, _solDot);
  // printResidual(_niElm, &_Be);
}

void feBilinearForm::computeMatrixFiniteDifference(feMetaNumber *metaNumber, feMesh *mesh,
                                                   feSolution *sol, int numElem)
{
  // feInfo("feBilinearForm::computeMatrixFiniteDifference");
  
  this->initialize(metaNumber, mesh, sol, numElem);
  setMatrixToZero(_niElm, _njElm, _Ae);

  double C0 = sol->getC0();
  // ==================================================================
  // Le résidu non perturbé
  // ==================================================================
  setResidualToZero(_niElm, &_R0);
  std::vector<double> &J = _cnc->getJacobians();
  _sysElm->computeBe(J, numElem, _intSpace, _geoSpace, _geoCoord, C0, sol->getCurrentTime(),
                     sol->getTimeStep(), _R0, _sol, _solDot);
  //_sysElm->computeAe(_intSpace, _geoSpace, _geoCoord, sol->getC0(), sol->getCurrentTime(), _Ae);
  
  // On boucle sur les fespace des inconnues
  // ==================================================================
  // Caclul des résidus perturbés
  // ==================================================================
  feInt numColumn = 0;
  for(size_t k = 0; k < _jVar.size(); k++) {
    feSpace *Unknowns = _intSpace[_jVar[k]];
    // On boucle sur les coefficients de la fonction d'interpolation
    for(feInt j = 0; j < Unknowns->getNbFunctions(); j++, numColumn++) {
      double temp_sol = _sol[_jVar[k]][j]; 
      double temp_soldot = _solDot[_jVar[k]][j];
      double delta_h = _h0 * std::max(fabs(temp_sol), 1.0);
      double invdelta_h = 1.0 / delta_h;

      _sol[_jVar[k]][j]    = _sol[_jVar[k]][j] + delta_h;
      _solDot[_jVar[k]][j] = _solDot[_jVar[k]][j] + delta_h * C0;

      setResidualToZero(_niElm, &_Rh);
      _sysElm->computeBe(J, numElem, _intSpace, _geoSpace, _geoCoord, C0, sol->getCurrentTime(),
                         sol->getTimeStep(), _Rh, _sol, _solDot);

      for(feInt i = 0; i < _niElm; i++) _Ae[i][numColumn] = (-_Rh[i] + _R0[i]) * invdelta_h;

      _sol[_jVar[k]][j]    = temp_sol;
      _solDot[_jVar[k]][j] = temp_soldot;
      // printMatrix(_niElm, _njElm, &_Ae);
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