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
                                  std::vector<double> &geoCoord, double c0, double tn, double *Be) {
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
                                   double *Be) {
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
                                           std::vector<double> &geoCoord, double c0, double tn,
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
                                   double *Be) {
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  int nFunctions = intSpace[_idU]->getNbFunctions();

  std::vector<double> x(3, 0.0);

  for(int k = 0; k < nG; ++k) {
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
    double S = _fct->eval(tn, x);

    for(int i = 0; i < nFunctions; ++i) {
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
      Be[i] -= _feU[i] * S * J[nG * numElem + k] * w[k];
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
                                      std::vector<double> &geoCoord, double c0, double tn,
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

    for(int i = 0; i < nFunctions; ++i) _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);

    for(int i = 0; i < nFunctions; ++i)
      for(int j = 0; j < nFunctions; ++j) Ae[i][j] += _feU[i] * rho * c0 * _feU[j] * jac * w[k];
  }
}

void feSysElm_1D_Masse::computeBe(std::vector<double> &J, int numElem,
                                  std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                  std::vector<double> &geoCoord, double c0, double tn, double *Be) {
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
                                   double *Be) {
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
                                      double *Be) {
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
                                   std::vector<double> &geoCoord, double c0, double tn,
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
                                         std::vector<double> &geoCoord, double c0, double tn,
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