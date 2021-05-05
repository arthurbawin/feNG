#include <feSysElm.h>

void feSysElm::createElementarySystem(std::vector<feSpace*> &space){
  _iVar[0] = 0;
  _jTest[0] = 0;
}

void feSysElm_1D_Source::createElementarySystem(std::vector<feSpace*> &space){
  _idU = 0;
  _iVar[0] = _idU;
  _jTest[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
}

void feSysElm_1D_Source::computeBe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, 
  std::vector<double> &geoCoord, double c0, double tn, std::vector<double> &Be)
{
  int                nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double            rho = _par;
  int        nFunctions = intSpace[_idU]->getNbFunctions();
  int               dim = intSpace[_idU]->getDim();
  int      nNodePerElem = intSpace[_idU]->getNbNodePerElem();

  double J, x;
  for(int k = 0; k < nG; ++k){
    std::vector<double> xGeoCoord(nNodePerElem,0.0); 
    for(int i = 0; i < nNodePerElem; ++i)
      xGeoCoord[i] = geoCoord[dim*i];

    J = geoSpace->interpolateFieldAtQuadNode_rDerivative(xGeoCoord, k);
    x = geoSpace->interpolateFieldAtQuadNode(xGeoCoord, k);

    double S = _fct(tn, {x, 0., 0.});

    for(int i = 0; i < nFunctions; ++i){
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i,k);
      Be[i] -= _feU[i] * S * J * w[k];
    }
  }
}

void feSysElm_1D_Diffusion::createElementarySystem(std::vector<feSpace*> &space){
  _idU = 0;
  _iVar[0] = _idU;
  _jTest[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
}

void feSysElm_1D_Diffusion::computeAe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, 
  std::vector<double> &geoCoord, double c0, double tn, std::vector<double> &Ae){
  int                nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double             kD = _par;
  int        nFunctions = intSpace[_idU]->getNbFunctions();
  int               dim = intSpace[_idU]->getDim();
  int      nNodePerElem = intSpace[_idU]->getNbNodePerElem();

  double J, dudx;
  for(int k = 0; k < nG; ++k){
    std::vector<double> xGeoCoord(nNodePerElem,0.0); 
    for(int i = 0; i < nNodePerElem; ++i)
      xGeoCoord[i] = geoCoord[dim*i];

    J = geoSpace->interpolateFieldAtQuadNode_rDerivative(xGeoCoord, k);

    for(int i = 0; i < nFunctions; ++i)
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i,k) / J;

    for(int i = 0; i < nFunctions; ++i)
      for(int j = 0; j < nFunctions; ++j)
        Ae[nFunctions*i+j] += _feUdx[i] * kD * _feUdx[j] * J * w[k];
  }
}

void feSysElm_1D_Diffusion::computeBe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, 
  std::vector<double> &geoCoord, double c0, double tn, std::vector<double> &Be)
{
  int                nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double             kD = _par;
  int        nFunctions = intSpace[_idU]->getNbFunctions();
  int               dim = intSpace[_idU]->getDim();
  int      nNodePerElem = intSpace[_idU]->getNbNodePerElem();

  double J, dudx;
  for(int k = 0; k < nG; ++k){
    std::vector<double> xGeoCoord(nNodePerElem,0.0); 
    for(int i = 0; i < nNodePerElem; ++i)
      xGeoCoord[i] = geoCoord[dim*i];

    J = geoSpace->interpolateFieldAtQuadNode_rDerivative(xGeoCoord, k);
    dudx = intSpace[_idU]->interpolateSolutionAtQuadNode_rDerivative(k);
    dudx /= J;

    for(int i = 0; i < nFunctions; ++i){
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i,k);
      _feUdx[i] /= J;
      Be[i] -= _feUdx[i] * kD * dudx * J * w[k];
    }
  }
}