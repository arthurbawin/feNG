#include "feSpace.h"
#include "feMesh.h"
#include "feNumber.h"
#include "feSolution.h"

feSpace::feSpace(feMesh *mesh, std::string fieldID, std::string cncGeoID, feFunction *fct)
  : _mesh(mesh), _fieldID(fieldID), _fieldTag(-1), _cncGeoID(cncGeoID), _cncGeoTag(-1), 
  _nQuad(-1), _nFunctions(0), _fct(fct)
{
  if(mesh != nullptr){ // Maillage existe et le fespace est associe a la bonne cncGeo
    _cncGeoTag = mesh->getCncGeoTag(cncGeoID);
    if(_cncGeoTag == -1)
      printf("In feSpace::feSpace : Error - The geometric connectivity \"%s\" does not exist in the mesh.\n", cncGeoID.c_str());
  }
};

int feSpace::getDim(){
  return _mesh->getDim(); // This is constant for all spaces on the same mesh...
  // return this->getCncGeo()->getDim();
}

int feSpace::getNbElm(){
  return _mesh->getNbElm(_cncGeoTag);
}

int feSpace::getNbNodePerElem(){
  return _mesh->getNbNodePerElem(_cncGeoTag);
}

feCncGeo* const feSpace::getCncGeo(){
  return _mesh->getCncGeoByTag(_cncGeoTag);
}

void feSpace::setQuadratureRule(feQuadrature *quad){
  int dim = quad->getDim();
  if(_nQuad == 1){
    printf("Warning : nQuad = 1 was already set.\n");
  }
  _nQuad = quad->getNQuad();
  _wQuad = quad->getWeights();
  _xQuad = quad->getPoints();
  _yQuad = quad->getYPoints();
  _zQuad = quad->getZPoints();

  _L.resize(_nFunctions*_nQuad);
  _dLdr.resize(_nFunctions*_nQuad);
  _dLds.resize(_nFunctions*_nQuad);
  _dLdt.resize(_nFunctions*_nQuad);

  for(size_t i = 0; i < _xQuad.size(); ++i){
    double r[3] = {_xQuad[i], _yQuad[i], _zQuad[i]};
    std::vector<double>    l =    L(r);
    std::vector<double> dldr = dLdr(r);
    std::vector<double> dlds = dLds(r);
    std::vector<double> dldt = dLdt(r);
    for(int j = 0; j < _nFunctions; ++j)
      _L[_nFunctions*i+j] = l[j];
    for(int j = 0; j < _nFunctions; ++j)
      _dLdr[_nFunctions*i+j] = dldr[j];
    for(int j = 0; j < _nFunctions; ++j)
      _dLds[_nFunctions*i+j] = dlds[j];
    for(int j = 0; j < _nFunctions; ++j)
      _dLdt[_nFunctions*i+j] = dldt[j];
  }

  // If the space is a geometric interpolant, precompute jacobians of the elements
  if(_fieldID == "GEO"){
    this->getCncGeo()->computeJacobians();
  }
}

void feSpace::initializeSolution(feSolution *sol){
  _sol.resize(_adr.size());
  for(size_t i = 0; i < _adr.size(); ++i){
    _sol[i] = sol->getSolAtDOF(_adr[i]);
  }
}

void feSpace::initializeSolutionDot(feSolution *sol){
  _soldot.resize(_adr.size());
  for(size_t i = 0; i < _adr.size(); ++i)
    _soldot[i] = sol->getSolDotAtDOF(_adr[i]);
}

double feSpace::interpolateSolutionAtQuadNode(int iNode){
  double res = 0.0;
  for(int i = 0; i < _nFunctions; ++i)
    res += _sol[i]*_L[_nFunctions*iNode+i];
  return res;
}

double feSpace::interpolateSolutionDotAtQuadNode(int iNode){
  double res = 0.0;
  for(int i = 0; i < _nFunctions; ++i)
    res += _soldot[i]*_L[_nFunctions*iNode+i];
  return res;
}

double feSpace::interpolateSolutionAtQuadNode_rDerivative(int iNode){
  double res = 0.0;
  for(int i = 0; i < _nFunctions; ++i)
    res += _sol[i]*_dLdr[_nFunctions*iNode+i];
  return res;
}

double feSpace::interpolateSolutionAtQuadNode_sDerivative(int iNode){
  double res = 0.0;
  for(int i = 0; i < _nFunctions; ++i)
    res += _sol[i]*_dLds[_nFunctions*iNode+i];
  return res;
}

double feSpace::interpolateField(std::vector<double> field, double r[3]){
  double res = 0.0;
  if(field.size() != (unsigned) _nFunctions){
    printf(" In feSpace::interpolateField : Erreur - Nombre de valeurs nodales non compatible avec le nombre d'interpolants de l'espace.\n");
    return res;
  }
  for(int i = 0; i < _nFunctions; ++i)
    res += field[i]*L(r)[i];
  return res;
}

double feSpace::interpolateFieldAtQuadNode(std::vector<double> field, int iNode){
  double res = 0.0;
  if(field.size() != (unsigned) _nFunctions){
    printf(" In feSpace::interpolateFieldAtQuadNode : Erreur - Nombre de valeurs nodales non compatible avec le nombre d'interpolants de l'espace.\n");
    return res;
  }
  for(int i = 0; i < _nFunctions; ++i)
    res += field[i]*_L[_nFunctions*iNode+i];
  return res;
}

double feSpace::interpolateField_rDerivative(std::vector<double> field, double r[3]){
  double res = 0.0;
  if(field.size() != (unsigned) _nFunctions){
    printf(" In feSpace::interpolateField : Erreur - Nombre de valeurs nodales non compatible avec le nombre d'interpolants de l'espace.\n");
    return res;
  }
  for(int i = 0; i < _nFunctions; ++i)
    res += field[i]*dLdr(r)[i];
  return res;
}

double feSpace::interpolateFieldAtQuadNode_rDerivative(std::vector<double> field, int iNode){
  double res = 0.0;
  if(field.size() != (unsigned) _nFunctions){
    printf(" In feSpace::interpolateFieldAtQuadNode : Erreur - Nombre de valeurs nodales non compatible avec le nombre d'interpolants de l'espace.\n");
    return res;
  }
  for(int i = 0; i < _nFunctions; ++i)
    res += field[i]*_dLdr[_nFunctions*iNode+i];
  return res;
}

void feSpace::interpolateVectorField(std::vector<double> field, double r[3], std::vector<double>& res){
  // Field structure :
  // [fx0 fy0 fz0 fx1 fy1 fz1 ... fxn fyn fzn]
  res[0] = res[1] = res[2] = 0.0;
  if(field.size() != 3 * (unsigned) _nFunctions){
    printf(" In feSpace::interpolateVectorField : Erreur - Nombre de valeurs nodales non compatible avec le nombre d'interpolants de l'espace.\n");
    return;
  }
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < _nFunctions; ++j){
      // res[i] += field[3*i+j]*L(r)[j];
      res[i] += field[3*j+i]*L(r)[j];
    }
  }
}

void feSpace::interpolateVectorFieldAtQuadNode(std::vector<double> field, int iNode, std::vector<double>& res){
  // Field structure :
  // [fx0 fy0 fz0 fx1 fy1 fz1 ... fxn fyn fzn]
  res[0] = res[1] = res[2] = 0.0;
  if(field.size() != 3 * (unsigned) _nFunctions){
    printf(" In feSpace::interpolateVectorFieldAtQuadNode : Erreur - Nombre de valeurs nodales non compatible avec le nombre d'interpolants de l'espace.\n");
    return;
  }
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < _nFunctions; ++j){
      // res[i] += field[3*i+j]*_L[_nFunctions*iNode+j];
      res[i] += field[3*j+i]*_L[_nFunctions*iNode+j];
    }
  }
}

void feSpace::interpolateVectorFieldAtQuadNode_rDerivative(std::vector<double> field, int iNode, std::vector<double>& res){
  // Field structure :
  // [fx0 fy0 fz0 fx1 fy1 fz1 ... fxn fyn fzn]
  res[0] = res[1] = res[2] = 0.0;
  if(field.size() != 3 * (unsigned) _nFunctions){
    printf(" In feSpace::interpolateVectorFieldAtQuadNode_rDerivative : Erreur - Nombre de valeurs nodales non compatible avec le nombre d'interpolants de l'espace.\n");
    return;
  }
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < _nFunctions; ++j){
      // res[i] += field[3*i+j]*_dLdr[_nFunctions*iNode+j];
      res[i] += field[3*j+i]*_dLdr[_nFunctions*iNode+j];
    }
  }
}

void feSpace::interpolateVectorFieldAtQuadNode_sDerivative(std::vector<double> field, int iNode, std::vector<double>& res){
  // Field structure :
  // [fx0 fy0 fz0 fx1 fy1 fz1 ... fxn fyn fzn]
  res[0] = res[1] = res[2] = 0.0;
  if(field.size() != 3 * (unsigned) _nFunctions){
    printf(" In feSpace::interpolateVectorFieldAtQuadNode_sDerivative : Erreur - Nombre de valeurs nodales non compatible avec le nombre d'interpolants de l'espace.\n");
    return;
  }
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < _nFunctions; ++j){
      // res[i] += field[3*i+j]*_dLdr[_nFunctions*iNode+j];
      res[i] += field[3*j+i]*_dLds[_nFunctions*iNode+j];
    }
  }
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
  _adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  _adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
}

void feSpace1DP2::initializeNumberingUnknowns(feNumber *number){
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i){
    number->defDDLSommet(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet(_mesh, _cncGeoID, i, 1);
    number->defDDLElement(_mesh, _cncGeoID, i, 1);
    // If the line is a boundary element, the edge should be set by the interior element
    // number->defDDLEdge(_mesh, _cncGeoID, i, 0, 0);
  }
}

void feSpace1DP2::initializeNumberingEssential(feNumber *number){
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i){
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 1);
    number->defDDLElement_essentialBC(_mesh, _cncGeoID, i);
    // Set essential BC on the edge if the line is a boundary element
    // If the line is an interior element, there is nothing in cncGeo->connecEdges
    number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 0);
  }
}

void feSpace1DP2::initializeAddressingVector(feNumber *number, int numElem){
  _adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  _adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
  _adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
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
    number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 0);
  }
}

void feSpace1DP3::initializeAddressingVector(feNumber *number, int numElem){
  _adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  _adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
  _adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
  _adr[3] = number->getDDLElement(_mesh, _cncGeoID, numElem, 1);
}

void feSpace1DP4::initializeNumberingUnknowns(feNumber *number){
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i){
    number->defDDLSommet(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet(_mesh, _cncGeoID, i, 1);
    number->defDDLElement(_mesh, _cncGeoID, i, 3);
  }
}

void feSpace1DP4::initializeNumberingEssential(feNumber *number){
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i){
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 1);
    number->defDDLElement_essentialBC(_mesh, _cncGeoID, i);
    number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 0);
  }
}

void feSpace1DP4::initializeAddressingVector(feNumber *number, int numElem){
  _adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  _adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
  _adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
  _adr[3] = number->getDDLElement(_mesh, _cncGeoID, numElem, 1);
  _adr[4] = number->getDDLElement(_mesh, _cncGeoID, numElem, 2);
}