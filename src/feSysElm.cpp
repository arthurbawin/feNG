#include <feSysElm.h>

void feSysElm::createElementarySystem(std::vector<feSpace*> &space){
  _iVar[0] = 0;
  _jVar[0] = 0;
}

void feSysElm_1D_Source::createElementarySystem(std::vector<feSpace*> &space){
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
}

void feSysElm_1D_Source::computeBe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, 
  std::vector<double> &geoCoord, double c0, double tn, double* Be){
  int                nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  int        nFunctions = intSpace[_idU]->getNbFunctions();
  int               dim = intSpace[_idU]->getDim();
  int      nNodePerElem = intSpace[_idU]->getNbNodePerElem();

  double J;
  for(int k = 0; k < nG; ++k){
    // std::vector<double> xGeoCoord(nNodePerElem,0.0); 
    // for(int i = 0; i < nNodePerElem; ++i)
    //   xGeoCoord[i] = geoCoord[dim*i];

    // J = geoSpace->interpolateFieldAtQuadNode_rDerivative(xGeoCoord, k);
    std::vector<double> j(3, 0.0);
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, j); // TODO : complete for dim > 1
    J = j[0];
    // x = geoSpace->interpolateFieldAtQuadNode(xGeoCoord, k);
    std::vector<double> x(3,0.0);
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

    // double S = _fct->eval(tn, {x, 0., 0.});
    double S = _fct->eval(tn, x);

    for(int i = 0; i < nFunctions; ++i){
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i,k);
      Be[i] -= _feU[i] * S * J * w[k];
    }
  }
}

void feSysElm_1D_Diffusion::createElementarySystem(std::vector<feSpace*> &space){
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
}

void feSysElm_1D_Diffusion::computeAe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, 
  std::vector<double> &geoCoord, double c0, double tn, double** Ae){
  int                nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double             kD = _par;
  int        nFunctions = intSpace[_idU]->getNbFunctions();
  int               dim = intSpace[_idU]->getDim();
  int      nNodePerElem = intSpace[_idU]->getNbNodePerElem();

  double J;
  for(int k = 0; k < nG; ++k){
    // std::vector<double> xGeoCoord(nNodePerElem,0.0); 
    // for(int i = 0; i < nNodePerElem; ++i)
    //   xGeoCoord[i] = geoCoord[dim*i];

    // J = geoSpace->interpolateFieldAtQuadNode_rDerivative(xGeoCoord, k);
    std::vector<double> j(3, 0.0);
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, j); // TODO : complete for dim > 1
    J = j[0];

    for(int i = 0; i < nFunctions; ++i)
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i,k) / J;

    for(int i = 0; i < nFunctions; ++i)
      for(int j = 0; j < nFunctions; ++j)
        Ae[i][j] += _feUdx[i] * kD * _feUdx[j] * J * w[k];
        // Ae[nFunctions*i+j] += _feUdx[i] * kD * _feUdx[j] * J * w[k];
  }
}

void feSysElm_1D_Diffusion::computeBe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, 
  std::vector<double> &geoCoord, double c0, double tn, double* Be){
  int                nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double             kD = _par;
  int        nFunctions = intSpace[_idU]->getNbFunctions();
  int               dim = intSpace[_idU]->getDim();
  int      nNodePerElem = intSpace[_idU]->getNbNodePerElem();

  double J, dudx;
  for(int k = 0; k < nG; ++k){
    // std::vector<double> xGeoCoord(nNodePerElem,0.0); 
    // for(int i = 0; i < nNodePerElem; ++i)
    //   xGeoCoord[i] = geoCoord[dim*i];

    // J = geoSpace->interpolateFieldAtQuadNode_rDerivative(xGeoCoord, k);
    std::vector<double> j(3, 0.0);
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, j); // TODO : complete for dim > 1
    J = j[0];
    dudx = intSpace[_idU]->interpolateSolutionAtQuadNode_rDerivative(k);
    dudx /= J;

    for(int i = 0; i < nFunctions; ++i){
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i,k);
      _feUdx[i] /= J;
      Be[i] -= _feUdx[i] * kD * dudx * J * w[k];
    }
  }
}

void feSysElm_1D_Masse::createElementarySystem(std::vector<feSpace*> &space){
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
}

void feSysElm_1D_Masse::computeAe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, 
  std::vector<double> &geoCoord, double c0, double tn, double** Ae){
  int                nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double            rho = _par;
  int        nFunctions = intSpace[_idU]->getNbFunctions();
  int               dim = intSpace[_idU]->getDim();
  int      nNodePerElem = intSpace[_idU]->getNbNodePerElem();

  double J;
  for(int k = 0; k < nG; ++k){
    // std::vector<double> xGeoCoord(nNodePerElem,0.0); 
    // for(int i = 0; i < nNodePerElem; ++i)
    //   xGeoCoord[i] = geoCoord[dim*i];

    // J = geoSpace->interpolateFieldAtQuadNode_rDerivative(xGeoCoord, k);
    std::vector<double> j(3, 0.0);
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, j); // TODO : complete for dim > 1
    J = j[0];

    for(int i = 0; i < nFunctions; ++i)
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i,k);

    for(int i = 0; i < nFunctions; ++i)
      for(int j = 0; j < nFunctions; ++j)
        Ae[i][j] += _feU[i] * rho * c0 * _feU[j] * J * w[k];
        // Ae[nFunctions*i+j] += _feU[i] * rho * c0 * _feU[j] * J * w[k];
  }
}

void feSysElm_1D_Masse::computeBe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, 
  std::vector<double> &geoCoord, double c0, double tn, double* Be){
  int                nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double            rho = _par;
  int        nFunctions = intSpace[_idU]->getNbFunctions();
  int               dim = intSpace[_idU]->getDim();
  int      nNodePerElem = intSpace[_idU]->getNbNodePerElem();

  double J, uDot;
  for(int k = 0; k < nG; ++k){
    // std::vector<double> xGeoCoord(nNodePerElem,0.0); 
    // for(int i = 0; i < nNodePerElem; ++i)
    //   xGeoCoord[i] = geoCoord[dim*i];

    // J = geoSpace->interpolateFieldAtQuadNode_rDerivative(xGeoCoord, k);
    std::vector<double> j(3, 0.0);
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, j); // TODO : complete for dim > 1
    J = j[0];
    uDot = intSpace[_idU]->interpolateSolutionDotAtQuadNode(k);

    for(int i = 0; i < nFunctions; ++i){
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i,k);
      Be[i] -= _feU[i] * rho * uDot * J * w[k];
    }
  }
}

void feSysElm_2D_Source::createElementarySystem(std::vector<feSpace*> &space){
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
  _feUdy.resize(space[_idU]->getNbFunctions());
}

void feSysElm_2D_Source::computeBe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, 
  std::vector<double> &geoCoord, double c0, double tn, double* Be){
  int                nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  int        nFunctions = intSpace[_idU]->getNbFunctions();
  // int               dim = intSpace[_idU]->getDim();
  int      nNodePerElem = intSpace[_idU]->getNbNodePerElem();

  double J;
  for(int k = 0; k < nG; ++k){
    std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
    std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr); // TODO : À discuter
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
    J = dxdr[0]*dxds[1] - dxdr[1]*dxds[0];
    
    std::vector<double> x(3,0.0);
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

    double S = _fct->eval(tn, x);

    for(int i = 0; i < nFunctions; ++i){
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i,k);
      Be[i] -= _feU[i] * S * J * w[k];
    }
  }
}

void feSysElm_2D_Diffusion::createElementarySystem(std::vector<feSpace*> &space){
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
  _feUdy.resize(space[_idU]->getNbFunctions());
}

void feSysElm_2D_Diffusion::computeAe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, 
  std::vector<double> &geoCoord, double c0, double tn, double** Ae){
  int                nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double             kD = _par;
  int        nFunctions = intSpace[_idU]->getNbFunctions();
  // int               dim = intSpace[_idU]->getDim();
  int      nNodePerElem = intSpace[_idU]->getNbNodePerElem();

  double J;
  for(int k = 0; k < nG; ++k){
    std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
    std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
    J = dxdr[0]*dxds[1] - dxdr[1]*dxds[0];

    double drdx =  dxds[1]/J;
    double drdy = -dxds[0]/J;
    double dsdx = -dxdr[1]/J;
    double dsdy =  dxdr[0]/J;

    for(int i = 0; i < nFunctions; ++i){
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i,k) * drdx + 
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i,k) * dsdx;
      _feUdy[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i,k) * drdy +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i,k) * dsdy;
    }

    for(int i = 0; i < nFunctions; ++i){
      for(int j = 0; j < nFunctions; ++j){
        Ae[i][j] += (_feUdx[i] * _feUdx[j] + _feUdy[i] * _feUdy[j]) * kD * J * w[k];
      }
    }
  }
}

void feSysElm_2D_Diffusion::computeBe(std::vector<feSpace*> &intSpace, feSpace *geoSpace, 
  std::vector<double> &geoCoord, double c0, double tn, double* Be)
{
  int                nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double             kD = _par;
  int        nFunctions = intSpace[_idU]->getNbFunctions();
  int               dim = intSpace[_idU]->getDim();
  int      nNodePerElem = intSpace[_idU]->getNbNodePerElem();

  double J, dudx, dudy;
  for(int k = 0; k < nG; ++k){
    std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
    std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
    J = dxdr[0]*dxds[1] - dxdr[1]*dxds[0];

    double drdx =  dxds[1]/J;
    double drdy = -dxds[0]/J;
    double dsdx = -dxdr[1]/J;
    double dsdy =  dxdr[0]/J;

    dudx = intSpace[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdx +
           intSpace[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdx;
    dudy = intSpace[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdy +
           intSpace[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdy;

    for(int i = 0; i < nFunctions; ++i){
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i,k) * drdx + 
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i,k) * dsdx;
      _feUdy[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i,k) * drdy +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i,k) * dsdy;
      Be[i] -= (_feUdx[i] * dudx + _feUdy[i] * dudy) * kD * J * w[k];
    }
  }
}