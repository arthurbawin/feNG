#include "feSpace.h"
#include "feMesh.h"
#include "feNumber.h"
#include "feSolution.h"

feSpace::feSpace(class feMesh *mesh, std::string fieldID, std::string cncGeoID,
  std::function<double(const double, const std::vector<double> &)> fct)
  : _mesh(mesh), _fieldID(fieldID), _fieldTag(-1), _cncGeoID(cncGeoID), _cncGeoTag(-1), 
  _fct(fct), _nQuad(-1), _nFunctions(0)
{
  if(mesh != nullptr) // Maillage existe et le fespace est associe a la bonne cncGeo
    _cncGeoTag = mesh->getCncGeoTag(cncGeoID);
};

int feSpace::getDim(){
  return _mesh->getDim();
}

int feSpace::getNbElm(){
  return _mesh->getNbElm(_cncGeoTag);
}

int feSpace::getNbNodePerElem(){
  return _mesh->getNbNodePerElem(_cncGeoTag);
}

void feSpace::setQuadratureRule(feQuadrature *quad){
  _nQuad = (_nQuad != -1) ? _nQuad : quad->getNQuad();
  _wQuad = quad->getWeights();
  _xQuad = quad->getPoints();

  _L.resize(_nFunctions*_nQuad);
  _dLdr.resize(_nFunctions*_nQuad); // TODO : multiplier par la dimension
  for(int i = 0; i < _xQuad.size(); ++i){
    double r[3] = {_xQuad[i], 0., 0.};
    std::vector<double>    l =    L(r);
    std::vector<double> dldr = dLdr(r);
    for(int j = 0; j < _nFunctions; ++j)
      _L[_nFunctions*i+j] = l[j];
    for(int j = 0; j < _nFunctions; ++j)
      _dLdr[_nFunctions*i+j] = dldr[j];
  }
}

void feSpace::initializeSolution(feSolution *sol){
  _sol.resize(_adr.size());
  for(int i = 0; i < _adr.size(); ++i)
    _sol[i] = sol->getSolAtDOF(_adr[i]);
}

void feSpace::initializeSolutionDot(feSolution *sol){
  _soldot.resize(_adr.size());
  for(int i = 0; i < _adr.size(); ++i)
    _soldot[i] = sol->getSolDotAtDOF(_adr[i]);
}

double feSpace::interpolateSolutionAtQuadNode(int iNode){
  double res = 0.0;
  for(int i = 0; i < _nFunctions; ++i)
    res += _sol[i]*_L[_nFunctions*iNode+i];
  return res;
}

double feSpace::interpolateSolutionAtQuadNode_rDerivative(int iNode){
  double res = 0.0;
  for(int i = 0; i < _nFunctions; ++i)
    res += _sol[i]*_dLdr[_nFunctions*iNode+i];
  return res;
}

double feSpace::interpolateField(std::vector<double> field, double r[3]){
  double res = 0.0;
  if(field.size() != _nFunctions){
    printf(" In feSpace::interpolateField : Erreur - Nombre de valeurs nodales non compatible avec le nombre d'interpolants de l'espace.\n");
    return res;
  }
  for(int i = 0; i < _nFunctions; ++i)
    res += field[i]*L(r)[i];
  return res;
}

double feSpace::interpolateFieldAtQuadNode(std::vector<double> field, int iNode){
  double res = 0.0;
  if(field.size() != _nFunctions){
    printf(" In feSpace::interpolateFieldAtQuadNode : Erreur - Nombre de valeurs nodales non compatible avec le nombre d'interpolants de l'espace.\n");
    return res;
  }
  for(int i = 0; i < _nFunctions; ++i)
    res += field[i]*_L[_nFunctions*iNode+i];
  return res;
}

double feSpace::interpolateField_rDerivative(std::vector<double> field, double r[3]){
  double res = 0.0;
  if(field.size() != _nFunctions){
    printf(" In feSpace::interpolateField : Erreur - Nombre de valeurs nodales non compatible avec le nombre d'interpolants de l'espace.\n");
    return res;
  }
  for(int i = 0; i < _nFunctions; ++i)
    res += field[i]*dLdr(r)[i];
  return res;
}

double feSpace::interpolateFieldAtQuadNode_rDerivative(std::vector<double> field, int iNode){
  double res = 0.0;
  if(field.size() != _nFunctions){
    printf(" In feSpace::interpolateFieldAtQuadNode : Erreur - Nombre de valeurs nodales non compatible avec le nombre d'interpolants de l'espace.\n");
    return res;
  }
  for(int i = 0; i < _nFunctions; ++i)
    res += field[i]*_dLdr[_nFunctions*iNode+i];
  return res;
}

void feSpace::printL(){
  std::cout<<"taille = "<<_L.size()<<std::endl;
  for(double l : _L)
    std::cout<<l<<" ";
  std::cout<<std::endl;
}

void feSpace::printdLdr(){
  std::cout<<"taille = "<<_dLdr.size()<<std::endl;
  for(double l : _dLdr)
    std::cout<<l<<" ";
  std::cout<<std::endl;
}

void feSpace1DP0::initializeNumberingUnknowns(feNumber *number){
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i)
    number->defDDLSommet(_mesh, _cncGeoID, i, 0);
}

void feSpace1DP0::initializeNumberingEssential(feNumber *number){
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i)
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 0);
}

void feSpace1DP0::initializeAddressingVector(feNumber *number, int numElem){
  // for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i)
  _adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
}

void feSpace1DP1::initializeNumberingUnknowns(feNumber *number){
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i){
    number->defDDLSommet(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet(_mesh, _cncGeoID, i, 1);
  }
}

void feSpace1DP1::initializeNumberingEssential(feNumber *number){
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i){
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 1);
  }
}

void feSpace1DP1::initializeAddressingVector(feNumber *number, int numElem){
  // _adr.resize(_nFunctions);
  // for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i){
    _adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
    _adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
  // }
}

void feSpace1DP2::initializeNumberingUnknowns(feNumber *number){
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i){
    number->defDDLSommet(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet(_mesh, _cncGeoID, i, 1);
    number->defDDLElement(_mesh, _cncGeoID, i, 1);
  }
}

void feSpace1DP2::initializeNumberingEssential(feNumber *number){
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i){
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 1);
    number->defDDLElement_essentialBC(_mesh, _cncGeoID, i);
  }
}

void feSpace1DP2::initializeAddressingVector(feNumber *number, int numElem){
  // _adr.resize(_nFunctions);
  // for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i){
    _adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
    _adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
    _adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
  // }
}

void feSpace1DP3::initializeNumberingUnknowns(feNumber *number){
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i){
    number->defDDLSommet(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet(_mesh, _cncGeoID, i, 1);
    number->defDDLElement(_mesh, _cncGeoID, i, 2);
  }
}

void feSpace1DP3::initializeNumberingEssential(feNumber *number){
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i){
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 1);
    number->defDDLElement_essentialBC(_mesh, _cncGeoID, i);
  }
}

void feSpace1DP3::initializeAddressingVector(feNumber *number, int numElem){
  // _adr.resize(_nFunctions);
  // for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i){
    _adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
    _adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
    _adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
    _adr[3] = number->getDDLElement(_mesh, _cncGeoID, numElem, 1);
  // }
}