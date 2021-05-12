#include "feBilinearForm.h"
#include "feQuadrature.h"

static inline void allocateMatrix(feInt m, feInt n, double*** A){
  *A = (double**) calloc(m, sizeof **A);
  if(*A == nullptr){
    printf("In feBilinearForm::allocateMatrix : Error - NULL pointer.\n");
    return;
  }
  for(feInt i = 0; i < m; ++i)
    (*A)[i] = (double*) calloc(n, sizeof ***A);
} 

static inline void setMatrixToZero(feInt m, feInt n, double*** A){
  for(feInt i = 0; i < m; ++i)
    for(feInt j = 0; j < n; ++j)
      (*A)[i][j] = 0.0;
};

static inline void printMatrix(feInt m, feInt n, double*** A){
  for(feInt i = 0; i < m; ++i)
    for(feInt j = 0; j < n; ++j)
      std::cout<<"A["<<i<<"]["<<j<<"] = "<<(*A)[i][j]<<std::endl;
};

static inline void freeMatrix(feInt m, double*** A){
  for(feInt i = 0; i < m; i++)
    free((*A)[i]);
  free(*A); 
}

static inline void allocateResidual(feInt m, double** b){
  *b = (double*) calloc(m, sizeof *b);
  if(*b == nullptr){
    printf("In feBilinearForm::allocateResidual : Error - NULL pointer.\n");
  }
} 

static inline void setResidualToZero(feInt m, double** b){
  for(feInt i = 0; i < m; ++i)
    (*b)[i] = 0.0;
};

static inline void printResidual(feInt m, double** b){
  for(feInt i = 0; i < m; ++i)
    std::cout<<"res["<<i<<"] = "<<(*b)[i]<<std::endl;
};

static inline void freeResidual(double** b){
  free(*b); 
}

feBilinearForm::feBilinearForm(std::vector<feSpace*> &space, feMesh *mesh, int nQuadraturePoints, feSysElm *sysElm)
  : _sysElm(sysElm), _intSpace(space), _cncGeoTag(space[0]->getCncGeoTag()),
  _geoSpace(mesh->getCncGeoByTag(_cncGeoTag)->getFeSpace()),
   _nGeoElm(mesh->getCncGeoByTag(_cncGeoTag)->getNbElm()), _nQuad(nQuadraturePoints)
{
  _nCoord = mesh->getDim();
  _nGeoNodes = mesh->getCncGeoByTag(_cncGeoTag)->getNbNodePerElem();
  _geoCoord.resize(_nCoord * _nGeoNodes);

  if(_geoSpace->getCncGeoTag() != _cncGeoTag)
    printf("In feBilinearForm::feBilinearForm : Erreur - Connectivité géométrique ne correspond pas à la connectivité de l'espace d'interpolation.\n");

  for(size_t i = 1; i < _intSpace.size(); ++i){
    if(_intSpace[i]->getCncGeoTag() != _cncGeoTag)
      printf("In feBilinearForm::feBilinearForm : Erreur - Forme définie sur plus d'une connectivité.\n");
  }
  
  // (Re-)initialize the interpolation functions at quadrature nodes
  feQuadrature2 *rule = new feQuadrature2(_nQuad);  // TODO : change this, choose the degree
  for(feSpace *fS : _intSpace)
    fS->setQuadratureRule(rule);
  _geoSpace->setQuadratureRule(rule);
  delete rule;

  _sysElm->createElementarySystem(_intSpace);
  _iVar = _sysElm->getIVar();
  _jVar = _sysElm->getJVar();

  // CREERVADIJ
  _niElm = 0;
  for(size_t k = 0; k <  _iVar.size(); ++k) _niElm += _intSpace[_iVar[k]]->getNbFunctions();
  _njElm = 0;
  for(size_t k = 0; k < _jVar.size(); ++k) _njElm += _intSpace[_jVar[k]]->getNbFunctions();
  _adrI.resize(_niElm);
  _adrJ.resize(_njElm);

  allocateMatrix(_niElm, _njElm, &_Ae);
  allocateResidual(_niElm, &_Be);
};

feBilinearForm::~feBilinearForm(){
    // if(_sysElm) delete _sysElm;
    freeMatrix(_niElm, &_Ae);
    freeResidual(&_Be);
  }

void feBilinearForm::initialize_vadij_only(feMetaNumber *metaNumber, int numElem){
  for(feSpace *fS : _intSpace){
    fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), numElem);
  }

  // INITIALISERVADIJ
  for(size_t i = 0, count = 0; i < _iVar.size(); ++i){
    int nielm = _intSpace[_iVar[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k)
      _adrI[count++] = _intSpace[_iVar[i]]->getAddressingVectorAt(k);
  }
  for(size_t i = 0, count = 0; i < _jVar.size(); ++i){
    int nielm = _intSpace[_jVar[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k)
      _adrJ[count++] = _intSpace[_jVar[i]]->getAddressingVectorAt(k);
  }
}

void feBilinearForm::initialize(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol, int numElem){
  for(feSpace *fS : _intSpace){
    fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), numElem);
    fS->initializeSolution(sol);
    fS->initializeSolutionDot(sol);
  }

  _geoCoord = mesh->getCoord(_cncGeoTag, numElem);

  // INITIALISERVADIJ
  for(size_t i = 0, count = 0; i < _iVar.size(); ++i){
    int nielm = _intSpace[_iVar[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k)
      _adrI[count++] = _intSpace[_iVar[i]]->getAddressingVectorAt(k);
  }
  for(size_t i = 0, count = 0; i < _jVar.size(); ++i){
    int nielm = _intSpace[_jVar[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k)
      _adrJ[count++] = _intSpace[_jVar[i]]->getAddressingVectorAt(k);
  }
}

void feBilinearForm::computeMatrix(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol, int numElem){
  this->initialize(metaNumber, mesh, sol, numElem);
  setMatrixToZero(_niElm, _njElm, &_Ae);
  _sysElm->computeAe(_intSpace, _geoSpace, _geoCoord, sol->getC0(), sol->getCurrentTime(), _Ae);
}

void feBilinearForm::computeResidual(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol, int numElem){
  this->initialize(metaNumber, mesh, sol, numElem);
  setResidualToZero(_niElm, &_Be);
  _sysElm->computeBe(_intSpace, _geoSpace, _geoCoord, sol->getC0(), sol->getCurrentTime(), _Be);
}