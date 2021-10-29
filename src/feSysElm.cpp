#include <feSysElm.h>

void feSysElm_0D_Masse::createElementarySystem(std::vector<feSpace *> &space) {
  // _idU = 0;
  _iVar.resize(space.size());
  _jVar.resize(space.size());
  for(int i = 0; i < space.size(); i++) {
    _iVar[i] = i;
    _jVar[i] = i;
  }
}

void feSysElm_0D_Masse::computeAe(std::vector<double> &J, int numElem,
                                  std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                  std::vector<double> &geoCoord, double c0, double tn,
                                  double **Ae) {
  // int                nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double rho = _par;
  // int        nFunctions = intSpace[_idU]->getNbFunctions();

  for(int i = 0; i < intSpace.size(); ++i) { Ae[i][i] = rho * c0; }
}

void feSysElm_0D_Masse::computeBe(std::vector<double> &J, int numElem,
                                  std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                  std::vector<double> &geoCoord, double c0, double tn, double dt, double *Be) {
  // int                nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double rho = _par;
  // int        nFunctions = intSpace[_idU]->getNbFunctions();

  double uDot;

  for(int i = 0; i < intSpace.size(); ++i) {
    uDot = intSpace[i]->interpolateSolutionDotAtQuadNode(0);
    Be[i] -= rho * uDot;
  }
}

void feSysElm_0D_Source::createElementarySystem(std::vector<feSpace *> &space) {
  _idU = 0;
  _feU.resize(space[_idU]->getNbFunctions());
}

void feSysElm_0D_Source::computeBe(std::vector<double> &J, int numElem,
                                   std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                   std::vector<double> &geoCoord, double c0, double tn,
                                   double dt, double *Be) {
  // int                nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  int nFunctions = intSpace[_idU]->getNbFunctions();
  std::vector<double> x(3, 0.0);
  double S = _fct->eval(tn, x);
  Be[0] -= S;
}

void feSysElm_0D_Source_crossed::createElementarySystem(std::vector<feSpace *> &space) {
  _iVar.resize(space.size());
  _jVar.resize(space.size());
  for(int i = 0; i < space.size(); i++) {
    _iVar[i] = i;
    _jVar[i] = i;
  }
}

void feSysElm_0D_Source_crossed::computeBe(std::vector<double> &J, int numElem,
                                           std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                           std::vector<double> &geoCoord, double c0, double tn, double dt,
                                           double *Be) {
  // int                nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  // int        nFunctions = intSpace[_idU]->getNbFunctions();
  std::vector<double> x(3, 0.0);
  std::vector<double> f(intSpace.size(), 0);
  if(_fct != nullptr) _fct->eval(tn, x, f);
  for(int i = 0; i < intSpace.size(); i++) { Be[i] = -f[i]; }
}

void feSysElm_1D_Source::createElementarySystem(std::vector<feSpace *> &space) {
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
}

void feSysElm_1D_Source::computeBe(std::vector<double> &J, int numElem,
                                   std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                   std::vector<double> &geoCoord, double c0, double tn,
                                   double dt, double *Be) {
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  int nFunctions = intSpace[_idU]->getNbFunctions();

  std::vector<double> x(3, 0.0);

  for(int k = 0; k < nG; ++k) {
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
    double S = _fct->eval(tn, x);
    // std::cout<< "-----S = "<<S<<"---------"<<std::endl;
    for(int i = 0; i < nFunctions; ++i) {
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
      Be[i] -= _feU[i] * S * J[nG * numElem + k] * w[k];
      // std::cout<<"feU["<<i<<"] = "<<_feU[i]<<std::endl;
      // std::cout<<"J = "<<J[nG * numElem + k]<<std::endl;
      // std::cout<<"w["<<k<<"] = "<<w[k]<<std::endl;
      // std::cout<<"Be["<<i<<"] = "<<Be[i]<<std::endl;
    }
  }
}

void feSysElm_1D_Diffusion::createElementarySystem(std::vector<feSpace *> &space) {
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
}

void feSysElm_1D_Diffusion::computeAe(std::vector<double> &J, int numElem,
                                      std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                      std::vector<double> &geoCoord, double c0, double tn,
                                      double **Ae) {
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double kD = _par;
  int nFunctions = intSpace[_idU]->getNbFunctions();

  double jac;
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * numElem + k];

    for(int i = 0; i < nFunctions; ++i)
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) / jac;

    for(int i = 0; i < nFunctions; ++i) {
      for(int j = 0; j < nFunctions; ++j) { Ae[i][j] += _feUdx[i] * kD * _feUdx[j] * jac * w[k]; }
    }
  }
}

void feSysElm_1D_Diffusion::computeBe(std::vector<double> &J, int numElem,
                                      std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                      std::vector<double> &geoCoord, double c0, double tn, double dt,
                                      double *Be) {
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double kD = _par;
  int nFunctions = intSpace[_idU]->getNbFunctions();

  double jac, dudx;
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * numElem + k];

    dudx = intSpace[_idU]->interpolateSolutionAtQuadNode_rDerivative(k);
    dudx /= jac;

    for(int i = 0; i < nFunctions; ++i) {
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k);
      _feUdx[i] /= jac;
      Be[i] -= _feUdx[i] * kD * dudx * jac * w[k];
      // std::cout<<"Be Diff["<<i<<"] = "<<Be[i]<<std::endl;
    }
  }
}

void feSysElm_1D_Masse::createElementarySystem(std::vector<feSpace *> &space) {
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
}

void feSysElm_1D_Masse::computeAe(std::vector<double> &J, int numElem,
                                  std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                  std::vector<double> &geoCoord, double c0, double tn,
                                  double **Ae) {
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double rho = _par;
  int nFunctions = intSpace[_idU]->getNbFunctions();

  double jac;
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * numElem + k];

    for(int i = 0; i < nFunctions; ++i){
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
    }

    for(int i = 0; i < nFunctions; ++i){
      for(int j = 0; j < nFunctions; ++j){
        Ae[i][j] += _feU[i] * rho * c0 * _feU[j] * jac * w[k];
      }
    }
  }
}

void feSysElm_1D_Masse::computeBe(std::vector<double> &J, int numElem,
                                  std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                  std::vector<double> &geoCoord, double c0, double tn, double dt, double *Be) {
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double rho = _par;
  int nFunctions = intSpace[_idU]->getNbFunctions();

  double jac, uDot;
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * numElem + k];

    uDot = intSpace[_idU]->interpolateSolutionDotAtQuadNode(k);

    for(int i = 0; i < nFunctions; ++i) {
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
      Be[i] -= _feU[i] * rho * uDot * jac * w[k];
      // std::cout<<"Be Mass["<<i<<"] = "<<Be[i]<<std::endl;
    }
  }
}

void feSysElm_2D_Source::createElementarySystem(std::vector<feSpace *> &space) {
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
  _feUdy.resize(space[_idU]->getNbFunctions());
}

void feSysElm_2D_Source::computeBe(std::vector<double> &J, int numElem,
                                   std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                   std::vector<double> &geoCoord, double c0, double tn,
                                   double dt, double *Be) {
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  int nFunctions = intSpace[_idU]->getNbFunctions();

  double jac;
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * numElem + k];

    std::vector<double> x(3, 0.0);
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

    double S = _fct->eval(tn, x);

    for(int i = 0; i < nFunctions; ++i) {
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
      Be[i] -= _feU[i] * S * jac * w[k];
    }
  }
}

void feSysElm_2D_Diffusion::createElementarySystem(std::vector<feSpace *> &space) {
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
  _feUdy.resize(space[_idU]->getNbFunctions());
}

void feSysElm_2D_Diffusion::computeAe(std::vector<double> &Ja, int numElem,
                                      std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                      std::vector<double> &geoCoord, double c0, double tn,
                                      double **Ae) {
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double kD = _par;
  int nFunctions = intSpace[_idU]->getNbFunctions();

  double J;
  for(int k = 0; k < nG; ++k) {
    std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
    std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
    J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];

    double drdx = dxds[1] / J;
    double drdy = -dxds[0] / J;
    double dsdx = -dxdr[1] / J;
    double dsdy = dxdr[0] / J;

    for(int i = 0; i < nFunctions; ++i) {
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
    }

    for(int i = 0; i < nFunctions; ++i) {
      for(int j = 0; j < nFunctions; ++j) {
        Ae[i][j] += (_feUdx[i] * _feUdx[j] + _feUdy[i] * _feUdy[j]) * kD * J * w[k];
      }
    }
  }
}

void feSysElm_2D_Diffusion::computeBe(std::vector<double> &Ja, int numElem,
                                      std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                      std::vector<double> &geoCoord, double c0, double tn,
                                      double dt, double *Be) {
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double kD = _par;
  int nFunctions = intSpace[_idU]->getNbFunctions();

  double J, dudx, dudy;
  for(int k = 0; k < nG; ++k) {
    std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
    std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
    J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];

    double drdx = dxds[1] / J;
    double drdy = -dxds[0] / J;
    double dsdx = -dxdr[1] / J;
    double dsdy = dxdr[0] / J;

    dudx = intSpace[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdx +
           intSpace[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdx;
    dudy = intSpace[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdy +
           intSpace[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdy;

    for(int i = 0; i < nFunctions; ++i) {
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
      Be[i] -= (_feUdx[i] * dudx + _feUdy[i] * dudy) * kD * J * w[k];
    }
  }
}

void feSysElm_2D_Advection::createElementarySystem(std::vector<feSpace *> &space) {
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
  _feUdy.resize(space[_idU]->getNbFunctions());
}

void feSysElm_2D_Advection::computeAe(std::vector<double> &Ja, int numElem,
                                      std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                      std::vector<double> &geoCoord, double c0, double tn,
                                      double **Ae) {
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double kD = _par;
  int nFunctions = intSpace[_idU]->getNbFunctions();

  double J;
  for(int k = 0; k < nG; ++k) {
    std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
    std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
    J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];

    double drdx = dxds[1] / J;
    double drdy = -dxds[0] / J;
    double dsdx = -dxdr[1] / J;
    double dsdy = dxdr[0] / J;

    for(int i = 0; i < nFunctions; ++i) {
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
    }

    for(int i = 0; i < nFunctions; ++i) {
      for(int j = 0; j < nFunctions; ++j) {
        Ae[i][j] += (_feUdx[i] * _feUdx[j] + _feUdy[i] * _feUdy[j]) * kD * J * w[k];
      }
    }
  }
}

void feSysElm_2D_Advection::computeBe(std::vector<double> &Ja, int numElem,
                                      std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                      std::vector<double> &geoCoord, double c0, double tn,
                                      double dt, double *Be) {
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double kD = _par;
  int nFunctions = intSpace[_idU]->getNbFunctions();

  double J, u, dudt, dudx, dudy;
  for(int k = 0; k < nG; ++k) {

    // Evaluate the exterior velocity field
    std::vector<double> v(2, 0.0);
    std::vector<double> x(3, 0.0);
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
    _fct->eval(tn, x, v);

    std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
    std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
    J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];

    double drdx = dxds[1] / J;
    double drdy = -dxds[0] / J;
    double dsdx = -dxdr[1] / J;
    double dsdy = dxdr[0] / J;

    // Compute SUPG parameter
    // double tau1 = 0.0, tau2 = 0.0, tau3 = 0.0, c = 0.0, cT = 0.0, kT = 0.0;
    // for(int i = 0; i < nFunctions; ++i) {
    //   dudt = intSpace[_idU]->interpolateSolutionDotAtQuadNode(k);
    //   dudx = intSpace[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdx +
    //          intSpace[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdx;
    //   dudy = intSpace[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdy +
    //          intSpace[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdy;

    //   _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
    //   _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
    //               intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
    //   _feUdy[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
    //               intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
    //   c  += J * w[k] * _feU[i] * (v[0] * dudx + v[1] * dudy);
    //   cT += J * w[k] * (v[0] * _feUdx[i] + v[1] * _feUdy[i]) * dudt;
    //   kT += J * w[k] * (v[0] * _feUdx[i] + v[1] * _feUdy[i]) * (v[0] * dudx + v[1] * dudy);
    // }
    double dx01 = fabs(geoCoord[3]-geoCoord[0]);
    double dx02 = fabs(geoCoord[6]-geoCoord[0]);
    double dx12 = fabs(geoCoord[6]-geoCoord[3]);
    double dy01 = fabs(geoCoord[4]-geoCoord[1]);
    double dy02 = fabs(geoCoord[7]-geoCoord[1]);
    double dy12 = fabs(geoCoord[7]-geoCoord[4]);
    double h01 = sqrt(dx01*dx01 + dy01*dy01);
    double h02 = sqrt(dx02*dx02 + dy02*dy02);
    double h12 = sqrt(dx12*dx12 + dy12*dy12);
    double h = fmax(h12,fmax(h01,h02));
    double normeU = sqrt(v[0]*v[0]+v[1]*v[1]);
    // double Peh = normeU * h / 2.0;

    // double tau = 1.0/( sqrt( 4.0/(dt*dt) + 4.0*normeU*normeU/h/h + 9.0*16.0/(h*h*h*h)) );
    double tau = 1.0/( sqrt( 4.0/(dt*dt) + 4.0*normeU*normeU/h/h ) );
    // tau = 0.01;
    // double tau = h/2.0/sqrt(v[0]*v[0]+v[1]*v[1]) * (cosh(Peh)/sinh(Peh) - 1.0/Peh);
    // double tau = (cosh(Peh)/sinh(Peh) - 1.0/Peh);
    // double tau = 0.01;
    // std::cout<<tau<<std::endl;
    // double deltaSUPG = c/kT;
    // double Re = (v[0]*v[0]+v[1]*v[1]) * c/kT;
    // printf("Re = %+-4.4e - c = %+-4.4e - cT = %+-4.4e - kT = %+-4.4e - delta = %+-4.4e\n", Re, c, cT, kT, deltaSUPG);

    for(int i = 0; i < nFunctions; ++i) {
      u = intSpace[_idU]->interpolateSolutionAtQuadNode(k);
      dudt = intSpace[_idU]->interpolateSolutionDotAtQuadNode(k);
      dudx = intSpace[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdx +
             intSpace[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdx;
      dudy = intSpace[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdy +
             intSpace[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdy;

      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
      // Be[i] -= (v[0] * _feUdx[i] + v[1] * _feUdy[i]) * u * J * w[k];
      // SUPG ?
      // double delta = geoCoord[3]-geoCoord[0];
      // double deltaSUPG = 0.01;
      // Be[i] -= ((v[0] * dudx + v[1] * dudy) * _feU[i] + delta * (v[0] * dudx + v[1] * dudy) * (v[0] * _feUdx[i] + v[1] * _feUdy[i])) * J * w[k];
      // Be[i] -= J * w[k] * ((v[0] * dudx + v[1] * dudy) * _feU[i] + deltaSUPG * (v[0] * _feUdx[i] + v[1] * _feUdy[i]) * (dudt + v[0] * dudx + v[1] * dudy));
      Be[i] -= J * w[k] * ((v[0] * dudx + v[1] * dudy) * _feU[i] + tau * (v[0] * _feUdx[i] + v[1] * _feUdy[i]) * (dudt + v[0] * dudx + v[1] * dudy));
    }
  }
}

void feSysElm_2D_Stokes::createElementarySystem(std::vector<feSpace *> &space) {
  _idU = 0;
  _idV = 1;
  _idP = 2;
  _iVar = {_idU, _idV, _idP};
  _jVar = {_idU, _idV, _idP};
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
  _feUdy.resize(space[_idU]->getNbFunctions());
  _feV.resize(space[_idV]->getNbFunctions());
  _feVdx.resize(space[_idV]->getNbFunctions());
  _feVdy.resize(space[_idV]->getNbFunctions());
  _feP.resize(space[_idP]->getNbFunctions());
}

void feSysElm_2D_Stokes::computeBe(std::vector<double> &Ja, int numElem,
                                   std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                   std::vector<double> &geoCoord, double c0, double tn, double dt,
                                   double *Be) {
  // TODO : Verifier que les règles d'intégration soient identiques ou adapter
  // On prend les poids de geoSpace, ok pour le jacobien seulement ?
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  // double             rho = _par[0];
  double mu = _par[1];
  int nFunctionsU = intSpace[_idU]->getNbFunctions();
  int nFunctionsV = intSpace[_idV]->getNbFunctions();
  int nFunctionsP = intSpace[_idP]->getNbFunctions();

  double J, u, v, p, dudx, dudy, dvdx, dvdy, Sxx, Sxy, Syx, Syy;
  std::vector<double> x(3, 0.);
  std::vector<double> f(2, 0.); // Forces volumiques

  for(int k = 0; k < nG; ++k) {
    // Forces volumiques
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
    if(_fct != nullptr) _fct->eval(tn, x, f);

    // Jacobien : à remplacer
    std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
    std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
    J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];

    double drdx = dxds[1] / J;
    double drdy = -dxds[0] / J;
    double dsdx = -dxdr[1] / J;
    double dsdy = dxdr[0] / J;

    u = intSpace[_idU]->interpolateSolutionAtQuadNode(k);
    v = intSpace[_idV]->interpolateSolutionAtQuadNode(k);
    p = intSpace[_idP]->interpolateSolutionAtQuadNode(k);

    dudx = intSpace[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdx +
           intSpace[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdx;
    dudy = intSpace[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdy +
           intSpace[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdy;
    dvdx = intSpace[_idV]->interpolateSolutionAtQuadNode_rDerivative(k) * drdx +
           intSpace[_idV]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdx;
    dvdy = intSpace[_idV]->interpolateSolutionAtQuadNode_rDerivative(k) * drdy +
           intSpace[_idV]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdy;

    Sxx = -p + 2. * mu * dudx;
    Sxy = mu * (dudy + dvdx);
    Syx = mu * (dudy + dvdx);
    Syy = -p + 2. * mu * dvdy;

    int cnt = 0;
    // Residu pour u
    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(
        i, k); // Fonction test de u : uniquement pour les forces volumiques
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
      Be[cnt++] -= (Sxx * _feUdx[i] + Sxy * _feUdy[i] - f[0] * _feU[i]) * J * w[k];
    }
    // Residu pour v
    for(int i = 0; i < nFunctionsV; ++i) {
      _feV[i] = intSpace[_idV]->getFunctionAtQuadNode(
        i, k); // Fonction test de v : uniquement pour les forces volumiques
      _feVdx[i] = intSpace[_idV]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  intSpace[_idV]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feVdy[i] = intSpace[_idV]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  intSpace[_idV]->getdFunctiondsAtQuadNode(i, k) * dsdy;
      Be[cnt++] -= (Syx * _feVdx[i] + Syy * _feVdy[i] - f[1] * _feV[i]) * J * w[k];
    }
    // Residu pour p
    for(int i = 0; i < nFunctionsP; ++i) {
      _feP[i] = intSpace[_idP]->getFunctionAtQuadNode(i, k);
      Be[cnt++] += _feP[i] * (dudx + dvdy) * J * w[k];
    }
  }
}

void feSysElm_2D_Stokes::computeAe(std::vector<double> &Ja, int numElem,
                                   std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                   std::vector<double> &geoCoord, double c0, double tn,
                                   double **Ae) {
  // TODO : Verifier que les règles d'intégration soient identiques ou adapter
  // On prend les poids de geoSpace, ok pour le jacobien seulement ?
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  // double             rho = _par[0];
  double mu = _par[1];
  int nFunctionsU = intSpace[_idU]->getNbFunctions();
  int nFunctionsV = intSpace[_idV]->getNbFunctions();
  int nFunctionsP = intSpace[_idP]->getNbFunctions();

  double jac, dudx, dudy, dvdx, dvdy;
  std::vector<double> x(3, 0.);
  std::vector<double> f(2, 0.); // Body forces

  for(int k = 0; k < nG; ++k) {
    // Body forces
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
    if(_fct != nullptr) _fct->eval(tn, x, f);

    // Jacobien : à remplacer
    std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
    std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
    jac = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];

    double drdx = dxds[1] / jac;
    double drdy = -dxds[0] / jac;
    double dsdx = -dxdr[1] / jac;
    double dsdy = dxdr[0] / jac;

    dudx = intSpace[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdx +
           intSpace[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdx;
    dudy = intSpace[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdy +
           intSpace[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdy;
    dvdx = intSpace[_idV]->interpolateSolutionAtQuadNode_rDerivative(k) * drdx +
           intSpace[_idV]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdx;
    dvdy = intSpace[_idV]->interpolateSolutionAtQuadNode_rDerivative(k) * drdy +
           intSpace[_idV]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdy;

    int I = 0, J;

    // Set up interpolation functions
    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
    }
    for(int i = 0; i < nFunctionsV; ++i) {
      _feV[i] = intSpace[_idV]->getFunctionAtQuadNode(i, k);
      _feVdx[i] = intSpace[_idV]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feVdy[i] = intSpace[_idV]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
    }
    for(int i = 0; i < nFunctionsP; ++i) { _feP[i] = intSpace[_idP]->getFunctionAtQuadNode(i, k); }

    // Équations pour u
    for(int i = 0; i < nFunctionsU; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) {
        Ae[I][J++] += jac * w[k] * (_feUdx[i] * mu * 2. * _feUdx[j] + _feUdy[i] * mu * _feUdy[j]);
      }
      for(int j = 0; j < nFunctionsV; ++j) {
        Ae[I][J++] += jac * w[k] * (_feUdy[i] * mu * _feVdx[j]);
      }
      for(int j = 0; j < nFunctionsP; ++j) { Ae[I][J++] -= jac * w[k] * (_feUdx[i] * _feP[j]); }
      ++I;
    }

    // Équations pour v
    for(int i = 0; i < nFunctionsV; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) {
        Ae[I][J++] += jac * w[k] * (_feVdx[i] * mu * _feUdy[j]);
      }
      for(int j = 0; j < nFunctionsV; ++j) {
        Ae[I][J++] += jac * w[k] * (_feVdy[i] * mu * 2. * _feVdy[j] + _feVdx[i] * mu * _feVdx[j]);
      }
      for(int j = 0; j < nFunctionsP; ++j) { Ae[I][J++] -= jac * w[k] * (_feVdy[i] * _feP[j]); }
      ++I;
    }

    // Équations pour p
    for(int i = 0; i < nFunctionsP; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) { Ae[I][J++] -= jac * w[k] * (_feP[i] * _feUdx[j]); }
      for(int j = 0; j < nFunctionsV; ++j) { Ae[I][J++] -= jac * w[k] * (_feP[i] * _feVdy[j]); }
      for(int j = 0; j < nFunctionsP; ++j) { Ae[I][J++] += 0.0; }
      ++I;
    }
  }
}

void feSysElm_2D_NavierStokes::createElementarySystem(std::vector<feSpace *> &space) {
  _idU = 0;
  _idV = 1;
  _idP = 2;
  _iVar = {_idU, _idV, _idP};
  _jVar = {_idU, _idV, _idP};
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
  _feUdy.resize(space[_idU]->getNbFunctions());
  _feV.resize(space[_idV]->getNbFunctions());
  _feVdx.resize(space[_idV]->getNbFunctions());
  _feVdy.resize(space[_idV]->getNbFunctions());
  _feP.resize(space[_idP]->getNbFunctions());
}

void feSysElm_2D_NavierStokes::computeBe(std::vector<double> &Ja, int numElem,
                                         std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                         std::vector<double> &geoCoord, double c0, double tn, double dt,
                                         double *Be) {
  // TODO : Verifier que les règles d'intégration soient identiques ou adapter
  // On prend les poids de geoSpace, ok pour le jacobien seulement ?
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double rho = _par[0];
  double mu = _par[1];
  int nFunctionsU = intSpace[_idU]->getNbFunctions();
  int nFunctionsV = intSpace[_idV]->getNbFunctions();
  int nFunctionsP = intSpace[_idP]->getNbFunctions();

  double J, u, v, p, dudt, dvdt, dudx, dudy, dvdx, dvdy, Sxx, Sxy, Syx, Syy;
  std::vector<double> x(3, 0.);
  std::vector<double> f(2, 0.); // Forces volumiques

  for(int k = 0; k < nG; ++k) {
    // Viscosité dynamique
    // double mu = (_viscosityFct != nullptr) ? _viscosityFct->eval(tn, x) : _par[1];
    // Forces volumiques
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
    if(_fct != nullptr) _fct->eval(tn, x, f);

    // Jacobien : à remplacer
    std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
    std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
    J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];

    double drdx = dxds[1] / J;
    double drdy = -dxds[0] / J;
    double dsdx = -dxdr[1] / J;
    double dsdy = dxdr[0] / J;

    u = intSpace[_idU]->interpolateSolutionAtQuadNode(k);
    v = intSpace[_idV]->interpolateSolutionAtQuadNode(k);
    p = intSpace[_idP]->interpolateSolutionAtQuadNode(k);

    dudt = intSpace[_idU]->interpolateSolutionDotAtQuadNode(k);
    dvdt = intSpace[_idV]->interpolateSolutionDotAtQuadNode(k);

    dudx = intSpace[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdx +
           intSpace[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdx;
    dudy = intSpace[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdy +
           intSpace[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdy;
    dvdx = intSpace[_idV]->interpolateSolutionAtQuadNode_rDerivative(k) * drdx +
           intSpace[_idV]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdx;
    dvdy = intSpace[_idV]->interpolateSolutionAtQuadNode_rDerivative(k) * drdy +
           intSpace[_idV]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdy;

    Sxx = -p + 2. * mu * dudx;
    Sxy = mu * (dudy + dvdx);
    Syx = mu * (dudy + dvdx);
    Syy = -p + 2. * mu * dvdy;

    // Test : Directional Do-Nothing condition


    int cnt = 0;
    // Residu pour u
    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
      Be[cnt++] -= (_feU[i] * rho * (dudt + u * dudx + v * dudy) + Sxx * _feUdx[i] +
                    Sxy * _feUdy[i] - f[0] * _feU[i]) *
                   J * w[k];
    }
    // Residu pour v
    for(int i = 0; i < nFunctionsV; ++i) {
      _feV[i] = intSpace[_idV]->getFunctionAtQuadNode(i, k);
      _feVdx[i] = intSpace[_idV]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  intSpace[_idV]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feVdy[i] = intSpace[_idV]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  intSpace[_idV]->getdFunctiondsAtQuadNode(i, k) * dsdy;
      Be[cnt++] -= (_feV[i] * rho * (dvdt + u * dvdx + v * dvdy) + Syx * _feVdx[i] +
                    Syy * _feVdy[i] - f[1] * _feV[i]) *
                   J * w[k];
    }
    // Residu pour p
    for(int i = 0; i < nFunctionsP; ++i) {
      _feP[i] = intSpace[_idP]->getFunctionAtQuadNode(i, k);
      Be[cnt++] += _feP[i] * (dudx + dvdy) * J * w[k];
    }
  }
}

void feSysElm_2D_NavierStokes::computeAe(std::vector<double> &Ja, int numElem,
                                         std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                         std::vector<double> &geoCoord, double c0, double tn,
                                         double **Ae) {
  // TODO : Verifier que les règles d'intégration soient identiques ou adapter
  // On prend les poids de geoSpace, ok pour le jacobien seulement ?
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double rho = _par[0];
  double mu = _par[1];
  int nFunctionsU = intSpace[_idU]->getNbFunctions();
  int nFunctionsV = intSpace[_idV]->getNbFunctions();
  int nFunctionsP = intSpace[_idP]->getNbFunctions();

  double jac, u, v, dudt, dvdt, dudx, dudy, dvdx, dvdy;
  std::vector<double> x(3, 0.);
  std::vector<double> f(2, 0.); // Forces volumiques

  for(int k = 0; k < nG; ++k) {
    // Viscosité dynamique
    // double mu = (_viscosityFct != nullptr) ? _viscosityFct->eval(tn, x) : _par[1];
    // Forces volumiques
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
    if(_fct != nullptr) _fct->eval(tn, x, f);

    // Jacobien : à remplacer
    std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
    std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
    jac = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];

    double drdx = dxds[1] / jac;
    double drdy = -dxds[0] / jac;
    double dsdx = -dxdr[1] / jac;
    double dsdy = dxdr[0] / jac;

    u = intSpace[_idU]->interpolateSolutionAtQuadNode(k);
    v = intSpace[_idV]->interpolateSolutionAtQuadNode(k);

    dudt = intSpace[_idU]->interpolateSolutionDotAtQuadNode(k);
    dvdt = intSpace[_idV]->interpolateSolutionDotAtQuadNode(k);

    dudx = intSpace[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdx +
           intSpace[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdx;
    dudy = intSpace[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdy +
           intSpace[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdy;
    dvdx = intSpace[_idV]->interpolateSolutionAtQuadNode_rDerivative(k) * drdx +
           intSpace[_idV]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdx;
    dvdy = intSpace[_idV]->interpolateSolutionAtQuadNode_rDerivative(k) * drdy +
           intSpace[_idV]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdy;

    int I = 0, J;

    // Set up interpolation functions
    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
    }
    for(int i = 0; i < nFunctionsV; ++i) {
      _feV[i] = intSpace[_idV]->getFunctionAtQuadNode(i, k);
      _feVdx[i] = intSpace[_idV]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feVdy[i] = intSpace[_idV]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
    }
    for(int i = 0; i < nFunctionsP; ++i) { _feP[i] = intSpace[_idP]->getFunctionAtQuadNode(i, k); }

    // Équations pour u
    for(int i = 0; i < nFunctionsU; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) {
        Ae[I][J++] +=
          jac * w[k] *
          (rho * _feU[i] * (c0 * _feU[j] + u * _feUdx[j] + v * _feUdy[j] + _feU[j] * dudx) +
           _feUdx[i] * mu * 2. * _feUdx[j] + _feUdy[i] * mu * _feUdy[j]);
      }
      for(int j = 0; j < nFunctionsV; ++j) {
        Ae[I][J++] += jac * w[k] * (rho * _feU[i] * (_feV[j] * dudy) + _feUdy[i] * mu * _feVdx[j]);
      }
      for(int j = 0; j < nFunctionsP; ++j) { Ae[I][J++] -= jac * w[k] * (_feUdx[i] * _feP[j]); }
      ++I;
    }

    // Équations pour v
    for(int i = 0; i < nFunctionsV; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) {
        Ae[I][J++] += jac * w[k] * (rho * _feV[i] * (_feU[j] * dvdx) + _feVdx[i] * mu * _feUdy[j]);
      }
      for(int j = 0; j < nFunctionsV; ++j) {
        Ae[I][J++] +=
          jac * w[k] *
          (rho * _feV[i] * (c0 * _feV[j] + u * _feVdx[j] + v * _feVdy[j] + _feV[j] * dvdy) +
           _feVdy[i] * mu * 2. * _feVdy[j] + _feVdx[i] * mu * _feVdx[j]);
      }
      for(int j = 0; j < nFunctionsP; ++j) { Ae[I][J++] -= jac * w[k] * (_feVdy[i] * _feP[j]); }
      ++I;
    }

    // Équations pour p
    for(int i = 0; i < nFunctionsP; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) { Ae[I][J++] -= jac * w[k] * (_feP[i] * _feUdx[j]); }
      for(int j = 0; j < nFunctionsV; ++j) { Ae[I][J++] -= jac * w[k] * (_feP[i] * _feVdy[j]); }
      for(int j = 0; j < nFunctionsP; ++j) { Ae[I][J++] += 0.0; }
      ++I;
    }
  }
}

void feSysElm_1D_DirectionalDoNothing::createElementarySystem(std::vector<feSpace *> &space) {
  _idU = 0;
  _idV = 1;
  _iVar = {_idU, _idV};
  _jVar = {_idU, _idV};
  _feU.resize(space[_idU]->getNbFunctions());
  _feV.resize(space[_idV]->getNbFunctions());
}

// Using finite differences instead
void feSysElm_1D_DirectionalDoNothing::computeAe(std::vector<double> &Ja, int numElem,
                                         std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                         std::vector<double> &geoCoord, double c0, double tn,
                                         double **Ae) {
}

void feSysElm_1D_DirectionalDoNothing::computeBe(std::vector<double> &Ja, int numElem,
                                         std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                         std::vector<double> &geoCoord, double c0, double tn, double dt,
                                         double *Be) {
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  int nFunctionsU = intSpace[_idU]->getNbFunctions();
  int nFunctionsV = intSpace[_idV]->getNbFunctions();

  double J, u, v, nx, ny, dx, dy, norm, udotn, s;

  for(int k = 0; k < nG; ++k) {
    J = Ja[nG * numElem + k];

    u = intSpace[_idU]->interpolateSolutionAtQuadNode(k);
    v = intSpace[_idV]->interpolateSolutionAtQuadNode(k);

    // Outward pointing normal vector to the boundary element
    // geoCoord = [x0 y0 z0 x1 y1 z1 ...], 0 and 1 are the extremities even for higher order lines
    dx = geoCoord[3] - geoCoord[0];
    dy = geoCoord[4] - geoCoord[1];
    nx = dy;
    ny = -dx;
    norm = sqrt(nx*nx + ny*ny);
    nx /= norm;
    ny /= norm;

    udotn = u*nx + v*ny;

    s = 0.5 * (1. - tanh(udotn/0.05));

    double alpha = 1.;
    double beta = 0.;
    double K = (u*u+v*v)/2.;

    int cnt = 0;
    // Residu pour u
    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
      // Be[cnt++] += 0.5 * s * _feU[i] * u * J * w[k];
      // Be[cnt++] += 0.5 * s * nx * _feU[i] * K * J * w[k];
      Be[cnt++] += 0.5 * (alpha * K * nx + beta * udotn * u) * _feU[i] * s * J * w[k];
    }
    // Residu pour v
    for(int i = 0; i < nFunctionsV; ++i) {
      _feV[i] = intSpace[_idV]->getFunctionAtQuadNode(i, k);
      // Be[cnt++] += 0.5 * s * _feV[i] * v * J * w[k];
      // Be[cnt++] += 0.5 * s * ny * _feV[i] * (u*u+v*v) * J * w[k];
      Be[cnt++] += 0.5 * (alpha * K * ny + beta * udotn * v) * _feV[i] * s * J * w[k];
    }
  }
}