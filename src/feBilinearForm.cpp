#include "feBilinearForm.h"

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

  for(int i = 1; i < _intSpace.size(); ++i){
    if(_intSpace[i]->getCncGeoTag() != _cncGeoTag)
      printf("In feBilinearForm::feBilinearForm : Erreur - Forme définie sur plus d'une connectivité.\n");
  }
  
  // (Re-)initialize the interpolation functions at quadrature nodes
  feQuadrature *rule = new feQuadrature(_nQuad);  // TODO : change this, choose the degree
  for(feSpace *fS : _intSpace)
    fS->setQuadratureRule(rule);
  _geoSpace->setQuadratureRule(rule);
  delete rule;

  _sysElm->createElementarySystem(_intSpace);
  _iVar = _sysElm->getIVar();
  _jVar = _sysElm->getJVar();

  // CREERVADIJ
  _niElm = 0;
  for(int k = 0; k <  _iVar.size(); ++k) _niElm += _intSpace[_iVar[k]]->getNbFunctions();
  _njElm = 0;
  for(int k = 0; k < _jVar.size(); ++k) _njElm += _intSpace[_jVar[k]]->getNbFunctions();
  _adrI.resize(_niElm);
  _adrJ.resize(_njElm);

  _Ae.resize(_niElm * _njElm);
  _Be.resize(_niElm);
};

void feBilinearForm::initialize(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol, int numElem){
  for(feSpace *fS : _intSpace){
    fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), numElem);
    fS->initializeSolution(sol);
    fS->initializeSolutionDot(sol);
  }

  _geoCoord = mesh->getCoord(_cncGeoTag, numElem);

  // INITIALISERVADIJ
  for(int i = 0, count = 0; i < _iVar.size(); ++i){
    int nielm = _intSpace[_iVar[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k)
      _adrI[count++] = _intSpace[_iVar[i]]->getAddressingVectorAt(k);
  }
  for(int i = 0, count = 0; i < _jVar.size(); ++i){
    int nielm = _intSpace[_jVar[i]]->getNbFunctions();
    for(int k = 0; k < nielm; ++k)
      _adrJ[count++] = _intSpace[_jVar[i]]->getAddressingVectorAt(k);
  }
}

void feBilinearForm::computeMatrix(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol, int numElem){
  this->initialize(metaNumber, mesh, sol, numElem);
  for(int i = 0; i < _niElm * _njElm; ++i)
    _Ae[i] = 0.0;
  _sysElm->computeAe(_intSpace, _geoSpace, _geoCoord, sol->getC0(), sol->getCurrentTime(), _Ae);
}

void feBilinearForm::computeRHS(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol, int numElem){
  this->initialize(metaNumber, mesh, sol, numElem);
  for(int i = 0; i < _niElm; ++i){
    _Be[i] = 0.0;
  }
  _sysElm->computeBe(_intSpace, _geoSpace, _geoCoord, sol->getC0(), sol->getCurrentTime(), _Be);
}