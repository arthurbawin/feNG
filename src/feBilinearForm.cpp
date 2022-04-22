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

  if(_geoSpace->getCncGeoTag() != _cncGeoTag)
    printf("In feBilinearForm::feBilinearForm : Erreur - Connectivité géométrique ne correspond "
           "pas à la connectivité de l'espace d'interpolation.\n");

  for(size_t i = 1; i < _intSpace.size(); ++i) {
    if(_intSpace[i]->getCncGeoTag() != _cncGeoTag)
      printf("In feBilinearForm::feBilinearForm : Erreur - Forme définie sur plus d'une "
             "connectivité.\n");
  }

  // CHANGED : The quadrature rules are set when each fespace is created
  // (Re-)initialize the interpolation functions at quadrature nodes
  // feQuadrature *rule =
  //   new feQuadrature(_degQuad, space[0]->getDim(), space[0]->getCncGeo()->getForme());

  // for(feSpace *fS : _intSpace) fS->setQuadratureRule(rule);
  // _geoSpace->setQuadratureRule(rule);
  // delete rule;

  _sysElm->createElementarySystem(_intSpace);
  _iVar = _sysElm->getIVar();
  _jVar = _sysElm->getJVar();

  // CREERVADIJ
  _niElm = _njElm = 0;
  for(size_t k = 0; k < _iVar.size(); ++k) _niElm += _intSpace[_iVar[k]]->getNbFunctions();
  for(size_t k = 0; k < _jVar.size(); ++k) _njElm += _intSpace[_jVar[k]]->getNbFunctions();

  _adrI.resize(_niElm);
  _adrJ.resize(_njElm);

  _Ae = allocateMatrix(_niElm, _njElm);
  allocateResidual(_niElm, &_Be);

  // ==================================================================
  // Pour le calcul par différences finies de la matrice
  // ==================================================================
  R0 = NULL;
  Rh = NULL;
  allocateResidual(_niElm, &Rh);
  allocateResidual(_niElm, &R0);
  h0 = pow(DBL_EPSILON, 1.0 / 2.0);
  // ==================================================================
  // Pour choisir la méthode de calcul de la matrice
  // ==================================================================
  ptrComputeMatrix = &feBilinearForm::computeMatrixFiniteDifference;
  if(sysElm->getMatrixAnalyticalStatus())
    ptrComputeMatrix = &feBilinearForm::computeMatrixAnalytical;
}

feBilinearForm::feBilinearForm(const feBilinearForm &f)
  : _intSpace(f._intSpace), _cnc(f._cnc), _cncGeoID(f._cncGeoID), _cncGeoTag(f._cncGeoTag),
    _geoSpace(f._geoSpace), _nGeoElm(f._nGeoElm), _degQuad(f._degQuad), _nCoord(f._nCoord),
    _nGeoNodes(f._nGeoNodes), _iVar(f._iVar), _jVar(f._jVar), _niElm(f._niElm), _njElm(f._njElm),
    _adrI(f._adrI), _adrJ(f._adrJ), R0(f.R0), Rh(f.Rh), h0(f.h0),
    ptrComputeMatrix(f.ptrComputeMatrix)
{
  _Ae = allocateMatrix(_niElm, _njElm);
  allocateResidual(_niElm, &_Be);

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
    // case MASSE_2D :
    //   _sysElm=new feSysElm_2D_Masse();
    //   break;
    case SOURCE_2D:
      _sysElm = new feSysElm_2D_Source(static_cast<feSysElm_2D_Source &>(*f._sysElm));
      break;
    case DIFFUSION_2D:
      _sysElm = new feSysElm_2D_Diffusion(static_cast<feSysElm_2D_Diffusion &>(*f._sysElm));
      break;
    case NEUMANN_1D:
      _sysElm = new feSysElm_1D_NeumannBC(static_cast<feSysElm_1D_NeumannBC &>(*f._sysElm));
      break;
    // case ADVECTION_2D :
    //   _sysElm=new feSysElm_2D_Advection();
    //   break;
    // case STOKES_2D :
    //   _sysElm=new feSysElm_2D_Stokes();
    //   break;
    // case NAVIERSTOKES_2D :
    //   _sysElm=new feSysElm_2D_NavierStokes());
    //   break;
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
  freeResidual(&Rh);
  freeResidual(&R0);
  delete _sysElm;
}

void feBilinearForm::initialize_vadij_only(feMetaNumber *metaNumber, int numElem)
{
  _adr.resize(_intSpace[0]->getNbFunctions());
  for(feSpace *fS : _intSpace) {
    fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), numElem, _adr);
    fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), numElem);
  }

  // INITIALISERVADIJ
  for(size_t i = 0, count = 0; i < _iVar.size(); ++i) {
    int nielm = _intSpace[_iVar[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k) _adrI[count++] = _adr[k];
    // for(int k = 0; k < nielm; ++k) _adrI[count++] =
    // _intSpace[_iVar[i]]->getAddressingVectorAt(k);
  }

  for(size_t i = 0, count = 0; i < _jVar.size(); ++i) {
    int nielm = _intSpace[_jVar[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k) _adrJ[count++] = _adr[k];
    // for(int k = 0; k < nielm; ++k) _adrJ[count++] =
    // _intSpace[_jVar[i]]->getAddressingVectorAt(k);
  }
}

void feBilinearForm::initialize(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol,
                                int numElem)
{
  _adr.resize(_intSpace[0]->getNbFunctions());
  for(feSpace *fS : _intSpace) {
    fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), numElem, _adr);
    fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), numElem);
    fS->initializeSolution(sol);
    fS->initializeSolutionDot(sol);
  }

  _sol.resize(_adr.size());
  std::vector<double> &solVec = sol->getSolutionReference();
  for(size_t i = 0; i < _adr.size(); ++i) {
    _sol[i] = solVec[_adr[i]];
  }

  _geoCoord = mesh->getCoord(_cncGeoTag, numElem);

  // INITIALISERVADIJ
  for(size_t i = 0, count = 0; i < _iVar.size(); ++i) {
    int nielm = _intSpace[_iVar[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k) {
      //_adrI[count++] = _intSpace[_iVar[i]]->getAddressingVectorAt(k);
      _adrI[count++] = _adr[k];
    };
  }
  for(size_t i = 0, count = 0; i < _jVar.size(); ++i) {
    int nielm = _intSpace[_jVar[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k) {
      //_adrJ[count++] = _intSpace[_jVar[i]]->getAddressingVectorAt(k);
      _adrJ[count++] = _adr[k];
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
                     sol->getCurrentTime(), _Ae);
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
                     sol->getCurrentTime(), sol->getTimeStep(), _Be, _sol);
  // #pragma omp critical
  // {
  // printf("%d;",numElem);for(int jj = 0; jj < 6; ++jj)std::cout<<_Be[jj]<<" ";
  // std::cout<<std::endl;
  // }
  // printResidual(_niElm, &_Be);
}

void feBilinearForm::computeMatrixFiniteDifference(feMetaNumber *metaNumber, feMesh *mesh,
                                                   feSolution *sol, int numElem)
{
  // printf("feBilinearForm::computeMatrixFiniteDifference\n");
  this->initialize(metaNumber, mesh, sol, numElem);
  // setMatrixToZero(_niElm, _njElm, &_Ae);
  setMatrixToZero(_niElm, _njElm, _Ae);

  double C0 = sol->getC0();
  // ==================================================================
  // Le résidu non perturbé
  // ==================================================================
  setResidualToZero(_niElm, &R0);
  std::vector<double> &J = _cnc->getJacobians();
  _sysElm->computeBe(J, numElem, _intSpace, _geoSpace, _geoCoord, C0, sol->getCurrentTime(),
                     sol->getTimeStep(), R0, _sol);
  //_sysElm->computeAe(_intSpace, _geoSpace, _geoCoord, sol->getC0(), sol->getCurrentTime(), _Ae);
  // On boucle sur les fespace des inconnues
  // ==================================================================
  // Caclul des résidus perturbés
  // ==================================================================
  feInt numColumn = 0;
  for(size_t k = 0; k < _jVar.size(); k++) {
    feSpace *Unknowns = _intSpace[_jVar[k]];
    std::vector<double> &solVec = Unknowns->getSolutionReference();
    std::vector<double> &soldotVec = Unknowns->getSolutionReferenceDot();
    // On boucle sur les coefficients de la fonction d'interpolation
    for(feInt j = 0; j < Unknowns->getNbFunctions(); j++, numColumn++) {
      double temp_sol = solVec[j];
      double temp_soldot = soldotVec[j];
      double delta_h = h0 * std::max(fabs(temp_sol), 1.0);
      double invdelta_h = 1.0 / delta_h;

      solVec[j] = solVec[j] + delta_h;
      soldotVec[j] = soldotVec[j] + delta_h * C0;

      setResidualToZero(_niElm, &Rh);
      _sysElm->computeBe(J, numElem, _intSpace, _geoSpace, _geoCoord, C0, sol->getCurrentTime(),
                         sol->getTimeStep(), Rh, _sol);

      for(feInt i = 0; i < _niElm; i++) _Ae[i][numColumn] = (-Rh[i] + R0[i]) * invdelta_h;

      solVec[j] = temp_sol;
      soldotVec[j] = temp_soldot;
    }
  }
  // for(auto val : _adrI)
  //   std::cout<<val<<" ";
  // std::cout<<std::endl;
  // for(auto val : _adrJ)
  //   std::cout<<val<<" ";
  // std::cout<<std::endl;
  // printMatrix(_niElm, _njElm, &_Ae);
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