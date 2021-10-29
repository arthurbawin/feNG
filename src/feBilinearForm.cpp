#include "feBilinearForm.h"
#include "feQuadrature.h"

// static inline void allocateMatrix(feInt m, feInt n, double*** A){
//   *A = (double**) calloc(m, sizeof **A);
//   if(*A == nullptr){
//     printf("In feBilinearForm::allocateMatrix : Error - NULL pointer.\n");
//     return;
//   }
//   for(feInt i = 0; i < m; ++i)
//     (*A)[i] = (double*) calloc(n, sizeof ***A);
// }

// static inline void setMatrixToZero(feInt m, feInt n, double*** A){
//   for(feInt i = 0; i < m; ++i)
//     for(feInt j = 0; j < n; ++j)
//       (*A)[i][j] = 0.0;
// };

// static inline void freeMatrix(feInt m, double*** A){
//   for(feInt i = 0; i < m; i++)
//     free((*A)[i]);
//   free(*A);
// }

static inline double **allocateMatrix(feInt m, feInt n) {
  double *p = new double[m * n];
  double **A = new double *[m];
  double *q = p;
  for(feInt i = 0; i < m; i++, q += n) A[i] = q;
  return A;
}

static inline void setMatrixToZero(feInt m, feInt n, double **A) {
  for(feInt i = 0; i < m * n; ++i) (*A)[i] = 0.0;
}

static inline void freeMatrix(feInt m, double **A) {
  delete[] A[0];
  delete[] A;
}

static inline void printMatrix(feInt m, feInt n, double ***A) {
  for(feInt i = 0; i < m; ++i)
    for(feInt j = 0; j < n; ++j) printf("A[%ld][%ld] = %+10.16e\n", i, j, (*A)[i][j]);
  // std::cout<<"A["<<i<<"]["<<j<<"] = "<<(*A)[i][j]<<std::endl;
}

static inline void allocateResidual(feInt m, double **b) {
  *b = (double *)calloc(m, sizeof *b);
  if(*b == nullptr) { printf("In feBilinearForm::allocateResidual : Error - NULL pointer.\n"); }
}

static inline void setResidualToZero(feInt m, double **b) {
  for(feInt i = 0; i < m; ++i) (*b)[i] = 0.0;
}

static inline void printResidual(feInt m, double **b) {
  for(feInt i = 0; i < m; ++i) printf("b[%ld] = %+10.16e\n", i, (*b)[i]);
}

static inline void freeResidual(double **b) { free(*b); }

feBilinearForm::feBilinearForm(std::vector<feSpace *> &space, feMesh *mesh, int degQuad,
                               feSysElm *sysElm)
  : _sysElm(sysElm), _intSpace(space), _cnc(space[0]->getCncGeo()), _cncGeoID(space[0]->getCncGeoID()),
    _cncGeoTag(space[0]->getCncGeoTag()), _geoSpace(_cnc->getFeSpace()), _nGeoElm(_cnc->getNbElm()),
    _degQuad(degQuad)

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

  // (Re-)initialize the interpolation functions at quadrature nodes
  feQuadrature *rule =
    new feQuadrature(_degQuad, space[0]->getDim(), space[0]->getCncGeo()->getForme());

  for(feSpace *fS : _intSpace) fS->setQuadratureRule(rule);
  _geoSpace->setQuadratureRule(rule);
  delete rule;

  _sysElm->createElementarySystem(_intSpace);
  _iVar = _sysElm->getIVar();
  _jVar = _sysElm->getJVar();

  // CREERVADIJ
  _niElm = 0;
  for(size_t k = 0; k < _iVar.size(); ++k) _niElm += _intSpace[_iVar[k]]->getNbFunctions();
  _njElm = 0;
  for(size_t k = 0; k < _jVar.size(); ++k) _njElm += _intSpace[_jVar[k]]->getNbFunctions();

  _adrI.resize(_niElm);
  _adrJ.resize(_njElm);

  _Ae = allocateMatrix(_niElm, _njElm);
  // allocateMatrix(_niElm, _njElm, &_Ae);
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

feBilinearForm::~feBilinearForm() {
  // freeMatrix(_niElm, &_Ae);
  freeMatrix(_niElm, _Ae);
  freeResidual(&_Be);
  freeResidual(&Rh);
  freeResidual(&R0);
  delete _sysElm;
}

void feBilinearForm::initialize_vadij_only(feMetaNumber *metaNumber, int numElem) {
  for(feSpace *fS : _intSpace) {
    fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), numElem);
  }

  // INITIALISERVADIJ
  for(size_t i = 0, count = 0; i < _iVar.size(); ++i) {
    int nielm = _intSpace[_iVar[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k) _adrI[count++] = _intSpace[_iVar[i]]->getAddressingVectorAt(k);
  }

  for(size_t i = 0, count = 0; i < _jVar.size(); ++i) {
    int nielm = _intSpace[_jVar[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k) _adrJ[count++] = _intSpace[_jVar[i]]->getAddressingVectorAt(k);
  }
}

void feBilinearForm::initialize(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol,
                                int numElem) {
  for(feSpace *fS : _intSpace) {
    fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), numElem);
    fS->initializeSolution(sol);
    fS->initializeSolutionDot(sol);
  }

  _geoCoord = mesh->getCoord(_cncGeoTag, numElem);

  // INITIALISERVADIJ
  for(size_t i = 0, count = 0; i < _iVar.size(); ++i) {
    int nielm = _intSpace[_iVar[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k) _adrI[count++] = _intSpace[_iVar[i]]->getAddressingVectorAt(k);
  }
  for(size_t i = 0, count = 0; i < _jVar.size(); ++i) {
    int nielm = _intSpace[_jVar[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k) _adrJ[count++] = _intSpace[_jVar[i]]->getAddressingVectorAt(k);
  }
}

void feBilinearForm::computeMatrix(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol,
                                   int numElem) {
  // this->initialize(metaNumber, mesh, sol, numElem);
  // setMatrixToZero(_niElm, _njElm, &_Ae);
  // _sysElm->computeAe(_intSpace, _geoSpace, _geoCoord, sol->getC0(), sol->getCurrentTime(), _Ae);
  // printMatrix(_niElm, _njElm, &_Ae);
  (this->*feBilinearForm::ptrComputeMatrix)(metaNumber, mesh, sol, numElem);
}

void feBilinearForm::computeMatrixAnalytical(feMetaNumber *metaNumber, feMesh *mesh,
                                             feSolution *sol, int numElem) {
  this->initialize(metaNumber, mesh, sol, numElem);
  // setMatrixToZero(_niElm, _njElm, &_Ae);
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
                                     int numElem) {
  this->initialize(metaNumber, mesh, sol, numElem);
  setResidualToZero(_niElm, &_Be);
  std::vector<double> &J = _cnc->getJacobians();
  _sysElm->computeBe(J, numElem, _intSpace, _geoSpace, _geoCoord, sol->getC0(),
                     sol->getCurrentTime(), sol->getTimeStep(), _Be);
  // printResidual(_niElm, &_Be);
}

void feBilinearForm::computeMatrixFiniteDifference(feMetaNumber *metaNumber, feMesh *mesh,
                                                   feSolution *sol, int numElem) {
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
  _sysElm->computeBe(J, numElem, _intSpace, _geoSpace, _geoCoord, C0, sol->getCurrentTime(), sol->getTimeStep(), R0);
  //_sysElm->computeAe(_intSpace, _geoSpace, _geoCoord, sol->getC0(), sol->getCurrentTime(), _Ae);
  // On boucle sur les fespace des inconnues
  // ==================================================================
  // Caclul des résidus perturbés
  // ==================================================================
  feInt numColumn = 0;
  for(size_t k = 0; k < _jVar.size(); k++) {
    feSpace *Unknowns = _intSpace[_jVar[k]];
    std::vector<double> &_sol = Unknowns->getSolutionReference();
    std::vector<double> &_soldot = Unknowns->getSolutionReferenceDot();
    // On boucle sur les coefficients de la fonction d'interpolation
    for(feInt j = 0; j < Unknowns->getNbFunctions(); j++, numColumn++) {
      double temp_sol = _sol[j];
      double temp_soldot = _soldot[j];
      double delta_h = h0 * std::max(fabs(temp_sol), 1.0);
      double invdelta_h = 1.0 / delta_h;

      _sol[j] = _sol[j] + delta_h;
      _soldot[j] = _soldot[j] + delta_h * C0;

      setResidualToZero(_niElm, &Rh);
      _sysElm->computeBe(J, numElem, _intSpace, _geoSpace, _geoCoord, C0, sol->getCurrentTime(), sol->getTimeStep(), Rh);

      for(feInt i = 0; i < _niElm; i++) _Ae[i][numColumn] = (-Rh[i] + R0[i]) * invdelta_h;

      _sol[j] = temp_sol;
      _soldot[j] = temp_soldot;
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

double feBilinearForm::getMatrixNorm() {
  int res = 0.0;
  for(feInt i = 0; i < _niElm; ++i)
    for(feInt j = 0; j < _njElm; ++j) res += _Ae[i][j] * _Ae[i][j];

  return sqrt(res);
}

double feBilinearForm::getResidualNorm() {
  int res = 0.0;
  for(feInt i = 0; i < _niElm; ++i) res += _Be[i] * _Be[i];

  return sqrt(res);
}