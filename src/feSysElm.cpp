#include "feSysElm.h"
#include "feBilinearForm.h"

// Migrated to feSysElm_0D, 1D and 2D on on Jan 6, 2023

void feSysElm_0D_StiffSpring::createElementarySystem(std::vector<feSpace *> &space)
{
  _idX = 0;
  _idV = 1;
  _iVar.resize(space.size());
  _jVar.resize(space.size());
  for(int i = 0; i < space.size(); i++) {
    _iVar[i] = i;
    _jVar[i] = i;
  }
}

void feSysElm_0D_StiffSpring::computeAe(std::vector<double> &J, int numElem,
                                        std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                        std::vector<double> &geoCoord, double c0, double tn,
                                        double **Ae, std::vector<std::vector<double>> &sol,
                                        std::vector<std::vector<double>> &solDot)
{
  double m = _par[0];
  double c = _par[1];
  double k = _par[2];
  Ae[0][0] = c0;
  Ae[0][1] = -1.;
  Ae[1][0] = k;
  Ae[1][1] = c0 * m + c;
}

void feSysElm_0D_StiffSpring::computeBe(std::vector<double> &J, int numElem,
                                        std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                        std::vector<double> &geoCoord, double c0, double tn,
                                        double dt, double *Be, std::vector<std::vector<double>> &sol,
                                        std::vector<std::vector<double>> &solDot,
                                        std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  double m = _par[0];
  double c = _par[1];
  double k = _par[2];
  double xDot = intSpace[_idX]->interpolateFieldAtQuadNode(solDot[_idX], 0);
  double vDot = intSpace[_idV]->interpolateFieldAtQuadNode(solDot[_idV], 0);
  double x = intSpace[_idX]->interpolateFieldAtQuadNode(sol[_idX], 0);
  double v = intSpace[_idV]->interpolateFieldAtQuadNode(sol[_idV], 0);
  Be[0] -= xDot - v;
  Be[1] -= m * vDot + x * k + c * v;
}

void feSysElm_0D_Stiff2::createElementarySystem(std::vector<feSpace *> &space)
{
  _idX = 0;
  _idY = 1;
  _idZ = 2;
  _iVar.resize(space.size());
  _jVar.resize(space.size());
  for(int i = 0; i < space.size(); i++) {
    _iVar[i] = i;
    _jVar[i] = i;
  }
}

void feSysElm_0D_Stiff2::computeAe(std::vector<double> &J, int numElem,
                                   std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                   std::vector<double> &geoCoord, double c0, double tn, double **Ae,
                                   std::vector<std::vector<double>> &sol, std::vector<std::vector<double>> &solDot)
{
  double omega = _par;
  Ae[0][0] = c0 + 1.;
  Ae[0][1] = Ae[0][2] = Ae[0][2] = Ae[1][0] = Ae[2][0] = 0.;
  Ae[1][1] = c0;
  Ae[2][2] = c0;
  Ae[2][1] = +omega;
  Ae[1][2] = -omega;
}

void feSysElm_0D_Stiff2::computeBe(std::vector<double> &J, int numElem,
                                   std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                   std::vector<double> &geoCoord, double c0, double tn, double dt,
                                   double *Be, std::vector<std::vector<double>> &sol,
                                   std::vector<std::vector<double>> &solDot,
                                   std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  double omega = _par;
  
  double xDot = intSpace[_idX]->interpolateFieldAtQuadNode(solDot[_idX], 0);
  double yDot = intSpace[_idY]->interpolateFieldAtQuadNode(solDot[_idY], 0);
  double zDot = intSpace[_idZ]->interpolateFieldAtQuadNode(solDot[_idZ], 0);
  double x = intSpace[_idX]->interpolateFieldAtQuadNode(sol[_idX], 0);
  double y = intSpace[_idY]->interpolateFieldAtQuadNode(sol[_idY], 0);
  double z = intSpace[_idZ]->interpolateFieldAtQuadNode(sol[_idZ], 0);
  Be[0] -= xDot + x;
  Be[1] -= yDot - omega * z;
  Be[2] -= zDot + omega * y;
}

void feSysElm_0D_Stiff3::createElementarySystem(std::vector<feSpace *> &space)
{
  _idX = 0;
  _idY = 1;
  _idZ = 2;
  _iVar.resize(space.size());
  _jVar.resize(space.size());
  for(int i = 0; i < space.size(); i++) {
    _iVar[i] = i;
    _jVar[i] = i;
  }
}

void feSysElm_0D_Stiff3::computeAe(std::vector<double> &J, int numElem,
                                   std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                   std::vector<double> &geoCoord, double c0, double tn, double **Ae,
                                   std::vector<std::vector<double>> &sol, std::vector<std::vector<double>> &solDot)
{
  double omega = _par;
  Ae[0][0] = c0 + 1.;
  Ae[1][0] = Ae[2][0] = 0.;
  Ae[0][1] = -1.;
  Ae[0][2] = -omega;
  Ae[1][1] = c0;
  Ae[2][2] = c0;
  Ae[2][1] = +omega;
  Ae[1][2] = -omega;
}

void feSysElm_0D_Stiff3::computeBe(std::vector<double> &J, int numElem,
                                   std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                   std::vector<double> &geoCoord, double c0, double tn, double dt,
                                   double *Be, std::vector<std::vector<double>> &sol,
                                   std::vector<std::vector<double>> &solDot,
                                   std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  double omega = _par;
  
  double xDot = intSpace[_idX]->interpolateFieldAtQuadNode(solDot[_idX], 0);
  double yDot = intSpace[_idY]->interpolateFieldAtQuadNode(solDot[_idY], 0);
  double zDot = intSpace[_idZ]->interpolateFieldAtQuadNode(solDot[_idZ], 0);
  double x = intSpace[_idX]->interpolateFieldAtQuadNode(sol[_idX], 0);
  double y = intSpace[_idY]->interpolateFieldAtQuadNode(sol[_idY], 0);
  double z = intSpace[_idZ]->interpolateFieldAtQuadNode(sol[_idZ], 0);
  Be[0] -= xDot + x - y - omega * z;
  Be[1] -= yDot - omega * z;
  Be[2] -= zDot + omega * y;
}

void feSysElm_0D_weakBC::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idL = 1;
  _iVar.resize(space.size());
  _jVar.resize(space.size());
  for(int i = 0; i < space.size(); i++) {
    _iVar[i] = i;
    _jVar[i] = i;
  }
}

void feSysElm_0D_weakBC::computeAe(std::vector<double> &J, int numElem,
                                   std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                   std::vector<double> &geoCoord, double c0, double tn, double **Ae,
                                   std::vector<std::vector<double>> &sol, std::vector<std::vector<double>> &solDot)
{
  Ae[0][1] = Ae[1][0] = 1.;
}

void feSysElm_0D_weakBC::computeBe(std::vector<double> &J, int numElem,
                                   std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                   std::vector<double> &geoCoord, double c0, double tn, double dt,
                                   double *Be, std::vector<std::vector<double>> &sol,
                                   std::vector<std::vector<double>> &solDot,
                                   std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  std::vector<double> x(3);
  geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, 0, x);
  double gamma = _fct->eval(tn, x);
  double lambda = intSpace[_idL]->interpolateFieldAtQuadNode(sol[_idL], 0);
  double u = intSpace[_idU]->interpolateFieldAtQuadNode(sol[_idU], 0);
  Be[0] -= lambda;
  Be[1] -= (u - gamma);
}

void feSysElm_0D_weakBC_edo1::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idL = 1;
  _idV = 2;
  _iVar.resize(space.size());
  _jVar.resize(space.size());
  for(int i = 0; i < space.size(); i++) {
    _iVar[i] = i;
    _jVar[i] = i;
  }
}

void feSysElm_0D_weakBC_edo1::computeAe(std::vector<double> &J, int numElem,
                                        std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                        std::vector<double> &geoCoord, double c0, double tn,
                                        double **Ae, std::vector<std::vector<double>> &sol,
                                        std::vector<std::vector<double>> &solDot)
{
  Ae[0][1] = Ae[1][0] = 1.;
  Ae[1][2] = -1.;
  Ae[2][2] = c0;
}

void feSysElm_0D_weakBC_edo1::computeBe(std::vector<double> &J, int numElem,
                                        std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                        std::vector<double> &geoCoord, double c0, double tn,
                                        double dt, double *Be, std::vector<std::vector<double>> &sol,
                                        std::vector<std::vector<double>> &solDot,
                                        std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  std::vector<double> x(3);
  geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, 0, x);
  double gammaDot = _fct->eval(tn, x);
  double lambda = intSpace[_idL]->interpolateFieldAtQuadNode(sol[_idL], 0);
  double u = intSpace[_idU]->interpolateFieldAtQuadNode(sol[_idU], 0);
  double v = intSpace[_idV]->interpolateFieldAtQuadNode(sol[_idV], 0);
  double vDot = intSpace[_idV]->interpolateFieldAtQuadNode(solDot[_idV], 0);
  // std::cout<< "v vaut "<< v << std::endl;
  // std::cout<<"gammaDot= "<<gammaDot<< " x= "<<x[0]<< "tn= " <<tn<<std::endl;
  // printf("%16.6e \n" , v);
  Be[0] -= lambda;
  Be[1] -= (u - v);
  Be[2] -= (vDot - gammaDot);
}

void feSysElm_0D_weakBC_edo1_V2::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idL = 1;
  _idV = 2;
  _iVar.resize(space.size());
  _jVar.resize(space.size());
  for(int i = 0; i < space.size(); i++) {
    _iVar[i] = i;
    _jVar[i] = i;
  }
}

void feSysElm_0D_weakBC_edo1_V2::computeAe(std::vector<double> &J, int numElem,
                                           std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                           std::vector<double> &geoCoord, double c0, double tn,
                                           double **Ae, std::vector<std::vector<double>> &sol,
                                           std::vector<std::vector<double>> &solDot)
{
  Ae[0][1] = Ae[1][0] = 1.;
  Ae[1][2] = -1.;
  Ae[2][2] = c0;
}

void feSysElm_0D_weakBC_edo1_V2::computeBe(std::vector<double> &J, int numElem,
                                           std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                           std::vector<double> &geoCoord, double c0, double tn,
                                           double dt, double *Be, std::vector<std::vector<double>> &sol,
                                           std::vector<std::vector<double>> &solDot,
                                           std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  std::vector<double> x(3);
  geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, 0, x);
  double gamma = _fct->eval(tn, x);
  double lambda = intSpace[_idL]->interpolateFieldAtQuadNode(sol[_idL], 0);
  double u = intSpace[_idU]->interpolateFieldAtQuadNode(sol[_idU], 0);
  double v = intSpace[_idV]->interpolateFieldAtQuadNode(sol[_idV], 0);
  double vDot = intSpace[_idV]->interpolateFieldAtQuadNode(solDot[_idV], 0);
  // std::cout<<"gamma= "<<gamma<< " x= "<<x[0]<< "tn= " <<tn<<std::endl;
  printf("%16.6e \n", v);
  Be[0] -= lambda;
  Be[1] -= (u - v - gamma);
  Be[2] -= (vDot);
}

void feSysElm_0D_weakBC_edo2::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idL = 1;
  _idV = 2;
  _idW = 3;
  _iVar.resize(space.size());
  _jVar.resize(space.size());
  for(int i = 0; i < space.size(); i++) {
    _iVar[i] = i;
    _jVar[i] = i;
  }
}

void feSysElm_0D_weakBC_edo2::computeAe(std::vector<double> &J, int numElem,
                                        std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                        std::vector<double> &geoCoord, double c0, double tn,
                                        double **Ae, std::vector<std::vector<double>> &sol,
                                        std::vector<std::vector<double>> &solDot)
{
  Ae[0][1] = Ae[1][0] = 1.;
  Ae[1][2] = Ae[2][3] = -1.;
  Ae[2][2] = Ae[3][3] = c0;
}

void feSysElm_0D_weakBC_edo2::computeBe(std::vector<double> &J, int numElem,
                                        std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                        std::vector<double> &geoCoord, double c0, double tn,
                                        double dt, double *Be, std::vector<std::vector<double>> &sol,
                                        std::vector<std::vector<double>> &solDot,
                                        std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  std::vector<double> x(3);
  geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, 0, x);
  double gammaDotDot = _fct->eval(tn, x);
  double lambda = intSpace[_idL]->interpolateFieldAtQuadNode(sol[_idL], 0);
  double u = intSpace[_idU]->interpolateFieldAtQuadNode(sol[_idU], 0);
  double v = intSpace[_idV]->interpolateFieldAtQuadNode(sol[_idV], 0);
  double vDot = intSpace[_idV]->interpolateFieldAtQuadNode(solDot[_idV], 0);
  double w = intSpace[_idW]->interpolateFieldAtQuadNode(sol[_idW], 0);
  double wDot = intSpace[_idW]->interpolateFieldAtQuadNode(solDot[_idW], 0);
  Be[0] -= lambda;
  Be[1] -= (u - v);
  Be[2] -= (vDot - w);
  Be[3] -= (wDot - gammaDotDot);
}

void feSysElm_0D_Masse::createElementarySystem(std::vector<feSpace *> &space)
{
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
                                  std::vector<double> &geoCoord, double c0, double tn, double **Ae,
                                  std::vector<std::vector<double>> &sol, std::vector<std::vector<double>> &solDot)
{
  // int                nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double rho = _par;
  // int        nFunctions = intSpace[_idU]->getNbFunctions();

  for(int i = 0; i < intSpace.size(); ++i) {
    Ae[i][i] = rho * c0;
  }
}

void feSysElm_0D_Masse::computeBe(std::vector<double> &J, int numElem,
                                  std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                  std::vector<double> &geoCoord, double c0, double tn, double dt,
                                  double *Be, std::vector<std::vector<double>> &sol, std::vector<std::vector<double>> &solDot,
                                  std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  // int                nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double rho = _par;
  // int        nFunctions = intSpace[_idU]->getNbFunctions();

  double uDot;

  for(int i = 0; i < intSpace.size(); ++i) {
    uDot = intSpace[i]->interpolateFieldAtQuadNode(solDot[i], 0);
    Be[i] -= rho * uDot;
  }
}

void feSysElm_0D_Source::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _feU.resize(space[_idU]->getNbFunctions());
}

void feSysElm_0D_Source::computeBe(std::vector<double> &J, int numElem,
                                   std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                   std::vector<double> &geoCoord, double c0, double tn, double dt,
                                   double *Be, std::vector<std::vector<double>> &sol,
                                   std::vector<std::vector<double>> &solDot,
                                   std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  std::vector<double> x(3, 0.0);
  double S = _fct->eval(tn, x);
  Be[0] -= S;
}

void feSysElm_0D_Source_crossed::createElementarySystem(std::vector<feSpace *> &space)
{
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
                                           double dt, double *Be, std::vector<std::vector<double>> &sol,
                                           std::vector<std::vector<double>> &solDot,
                                           std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  // int                nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  // int        nFunctions = intSpace[_idU]->getNbFunctions();
  std::vector<double> x(3, 0.0);
  std::vector<double> f(intSpace.size(), 0);
  if(_fct != nullptr) _fct->eval(tn, x, f);
  for(int i = 0; i < intSpace.size(); i++) {
    Be[i] = -f[i];
  }
}

void feSysElm_1D_weakBC_edo1::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idL = 1;
  _idV = 2;
  _iVar.resize(space.size());
  _jVar.resize(space.size());
  _feU.resize(space[_idU]->getNbFunctions());
  _feL.resize(space[_idL]->getNbFunctions());
  _feV.resize(space[_idV]->getNbFunctions());
  for(int i = 0; i < space.size(); i++) {
    _iVar[i] = i;
    _jVar[i] = i;
  }
}

// void feSysElm_1D_weakBC_edo1::computeAe(std::vector<double> &J, int numElem,
//                                         std::vector<feSpace *> &intSpace, feSpace *geoSpace,
//                                         std::vector<double> &geoCoord, double c0, double tn,
//                                         double **Ae, std::vector<std::vector<double>> &sol,
//                                         std::vector<std::vector<double>> &solDot) {
//   int nG = geoSpace->getNbQuadPoints();
//   std::vector<double> w = geoSpace->getQuadratureWeights();
//   int nFunctions = intSpace[_idU]->getNbFunctions();
//   double jac;

//   for(int k = 0; k < nG; ++k) {
//     std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
//     geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
//     jac = sqrt(dxdr[0] * dxdr[0]+dxdr[1] * dxdr[1]);
//     // jac = J[nG * numElem + k];
//     for(int i = 0; i < nFunctions; ++i){
//           _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
//           _feL[i] = intSpace[_idL]->getFunctionAtQuadNode(i, k);
//           _feV[i] = intSpace[_idV]->getFunctionAtQuadNode(i, k);
//     }
//     for(int i = 0; i < nFunctions; ++i){
//       for(int j = 0; j < nFunctions; ++j) {
//         Ae[i][nFunctions +j] += _feU[i]*_feL[j]*jac *w[k];
//       }
//     }
//     for(int i = 0; i < nFunctions; ++i){
//       for(int j = 0; j < nFunctions; ++j) {
//         Ae[nFunctions +i][j] += _feU[j]*_feL[i]*jac *w[k];
//       }
//     }
//     for(int i = 0; i < nFunctions; ++i){
//       for(int j = 0; j < nFunctions; ++j) {
//         Ae[nFunctions +i][2*nFunctions +j] += -_feL[i]*_feV[j]*jac *w[k];
//       }
//     }

//     for(int i = 0; i < nFunctions; ++i){
//       for(int j = 0; j < nFunctions; ++j) {
//         Ae[2*nFunctions +i][2*nFunctions +j]  += c0*_feV[i]*_feV[j] *jac *w[k];
//       }
//     }
//   }
// }

void feSysElm_1D_weakBC_edo1::computeAe(std::vector<double> &J, int numElem,
                                        std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                        std::vector<double> &geoCoord, double c0, double tn,
                                        double **Ae, std::vector<std::vector<double>> &sol,
                                        std::vector<std::vector<double>> &solDot)
{
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  int nFunctionsU = intSpace[_idU]->getNbFunctions();
  int nFunctionsL = intSpace[_idL]->getNbFunctions();
  double jac;

  for(int k = 0; k < nG; ++k) {
    std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
    jac = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);
    // jac = J[nG * numElem + k];
    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
      _feV[i] = intSpace[_idV]->getFunctionAtQuadNode(i, k);
    }

    for(int i = 0; i < nFunctionsL; ++i) {
      _feL[i] = intSpace[_idL]->getFunctionAtQuadNode(i, k);
    }
    for(int i = 0; i < nFunctionsU; ++i) {
      for(int j = 0; j < nFunctionsL; ++j) {
        Ae[i][nFunctionsU + j] += _feU[i] * _feL[j] * jac * w[k];
      }
    }
    for(int i = 0; i < nFunctionsL; ++i) {
      for(int j = 0; j < nFunctionsU; ++j) {
        Ae[nFunctionsU + i][j] += _feU[j] * _feL[i] * jac * w[k];
      }
    }
    for(int i = 0; i < nFunctionsL; ++i) {
      for(int j = 0; j < nFunctionsU; ++j) {
        Ae[nFunctionsU + i][(nFunctionsU + nFunctionsL) + j] += -_feL[i] * _feV[j] * jac * w[k];
      }
    }

    for(int i = 0; i < nFunctionsU; ++i) {
      for(int j = 0; j < nFunctionsU; ++j) {
        Ae[(nFunctionsU + nFunctionsL) + i][(nFunctionsU + nFunctionsL) + j] +=
          c0 * _feV[i] * _feV[j] * jac * w[k];
      }
    }
  }
}

// void feSysElm_1D_weakBC_edo1::computeBe(std::vector<double> &J, int numElem,
//                                         std::vector<feSpace *> &intSpace, feSpace *geoSpace,
//                                         std::vector<double> &geoCoord, double c0, double tn,
//                                         double *Be, std::vector<std::vector<double>> &sol, std::vector<double>
//                                         &solDot) {
//   int nG = geoSpace->getNbQuadPoints();
//   std::vector<double> w = geoSpace->getQuadratureWeights();
//   int nFunctions = intSpace[_idU]->getNbFunctions();

//   std::vector<double> x(3, 0.0);
//   double u, l, v, vDot, jac;
//   // std::cout<<"gammaDot= "<<gammaDot<< " x= "<<x[0]<< "tn= " <<tn<<std::endl;
//   // printf("%16.6e \n" , v);

//   for(int k = 0; k < nG; ++k) {
//     std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
//     geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
//     jac = sqrt(dxdr[0] * dxdr[0]+dxdr[1] * dxdr[1]);
//     geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
//     double gammaDot = _fct->eval(tn, x);
//     for(int i = 0; i < nFunctions; ++i) {
//       u = intSpace[_idU]->interpolateFieldAtQuadNode(sol, k);
//       v = intSpace[_idV]->interpolateFieldAtQuadNode(sol, k);
//       l = intSpace[_idL]->interpolateFieldAtQuadNode(sol, k);
//       vDot = intSpace[_idV]->interpolateFieldAtQuadNode(solDot, k);
//       _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
//       _feL[i] = intSpace[_idL]->getFunctionAtQuadNode(i, k);
//       _feV[i] = intSpace[_idV]->getFunctionAtQuadNode(i, k);

//       Be[i] -=_feU[i] * l * jac * w[k];
//       Be[nFunctions +i] -= _feL[i] *(u - v) * jac* w[k] ;
//       Be[2*nFunctions +i] -= _feV[i] *(vDot - gammaDot) * jac * w[k];
//     }
//   }

// }

void feSysElm_1D_weakBC_edo1::computeBe(std::vector<double> &J, int numElem,
                                        std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                        std::vector<double> &geoCoord, double c0, double tn,
                                        double dt, double *Be, std::vector<std::vector<double>> &sol,
                                        std::vector<std::vector<double>> &solDot,
                                        std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  // Make the difference between different type of space
  int nFunctionsU = intSpace[_idU]->getNbFunctions();
  int nFunctionsL = intSpace[_idL]->getNbFunctions();

  std::vector<double> x(3, 0.0);
  double u, l, v, vDot, jac;
  // std::cout<<"gammaDot= "<<gammaDot<< " x= "<<x[0]<< "tn= " <<tn<<std::endl;
  // printf("%16.6e \n" , v);
  for(int k = 0; k < nG; ++k) {
    std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
    jac = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
    double gammaDot = _fct->eval(tn, x);
    u = intSpace[_idU]->interpolateFieldAtQuadNode(sol[_idU], k);
    v = intSpace[_idV]->interpolateFieldAtQuadNode(sol[_idV], k);
    l = intSpace[_idL]->interpolateFieldAtQuadNode(sol[_idL], k);
    vDot = intSpace[_idV]->interpolateFieldAtQuadNode(solDot[_idV], k);
    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
      Be[i] -= _feU[i] * l * jac * w[k];
    }
    for(int i = 0; i < nFunctionsL; ++i) {
      _feL[i] = intSpace[_idL]->getFunctionAtQuadNode(i, k);
      Be[nFunctionsU + i] -= _feL[i] * (u - v) * jac * w[k];
    }
    for(int i = 0; i < nFunctionsU; ++i) {
      _feV[i] = intSpace[_idV]->getFunctionAtQuadNode(i, k);
      Be[(nFunctionsU + nFunctionsL) + i] -= _feV[i] * (vDot - gammaDot) * jac * w[k];
    }
  }
}
// void feSysElm_1D_weakBC_Vec::createElementarySystem(std::vector<feSpace *> &space) {
//   _idU = 0;
//   _idU = 1;
//   _idL = 2;
//   _idBu = 3;
//   _idBv = 4 ;
//   _iVar = {_idU, _idV, _idL , _idBu, _idBv};
//   _jVar = {_idU, _idV, _idL , _idBu, _idBv};
//   _feU.resize(space[_idU]->getNbFunctions());
//   _feV.resize(space[_idV]->getNbFunctions());
//   _feL.resize(space[_idL]->getNbFunctions());
//   _feBu.resize(space[_idBv]->getNbFunctions());
//   _feBv.resize(space[_idBu]->getNbFunctions());
// }

// void feSysElm_1D_weakBC_Vec::computeAe(std::vector<double> &J, int numElem,
//                                         std::vector<feSpace *> &intSpace, feSpace *geoSpace,
//                                         std::vector<double> &geoCoord, double c0, double tn,
//                                         double **Ae, std::vector<std::vector<double>> &sol,
//                                         std::vector<std::vector<double>> &solDot) {
//   int nG = geoSpace->getNbQuadPoints();
//   std::vector<double> w = geoSpace->getQuadratureWeights();
//   int nFunctionsU = intSpace[_idU]->getNbFunctions();
//   int nFunctionsV = intSpace[_idV]->getNbFunctions();
//   int nFunctionsL = intSpace[_idL]->getNbFunctions();
//   int nFunctionsBv = intSpace[_idBu]->getNbFunctions();
//   int nFunctionsBu = intSpace[_idBv]->getNbFunctions();
//   double jac;

//   for(int k = 0; k < nG; ++k) {
//     std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
//     geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
//     jac = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);
//     // jac = J[nG * numElem + k];
//     int I = 0, J;
//     // // Set up interpolation functions
//     for(int i = 0; i < nFunctionsU; ++i) {
//       _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
//       _feV[i] = intSpace[_idV]->getFunctionAtQuadNode(i, k);
//       _feBv[i] = intSpace[_idBu]->getFunctionAtQuadNode(i, k);
//       _feBu[i] = intSpace[_idBv]->getFunctionAtQuadNode(i, k);
//     }
//     for(int i = 0; i < nFunctionsL; ++i) { _feL[i] = intSpace[_idL]->getFunctionAtQuadNode(i, k);
//     }

//     for(int i = 0; i < nFunctionsU; ++i) {
//       J=nFunctionsU + nFunctionsV;
//       for(int j = 0; j < nFunctionsL; ++j) {
//         Ae[I][J++] += _feU[i] * _feL[j] * jac * w[k];
//       }
//       I++;
//     }
//     for(int i = 0; i < nFunctionsU; ++i) {
//       J=0;
//       I=nFunctionsU + nFunctionsV ;
//       for(int j = 0; j < nFunctionsL; ++j) {
//         Ae[I][J++] += _feU[j] * _feL[i] * jac * w[k];
//       }
//       I++;
//     }

//     for(int i = 0; i < nFunctionsV; ++i) {
//       J=nFunctionsU + nFunctionsV;
//       I=nFunctionsU;
//       for(int j = 0; j < nFunctionsL; ++j) {
//         Ae[I][J++] += _feV[i] * _feL[j] * jac * w[k];
//       }
//       I++;
//     }
//     for(int i = 0; i < nFunctionsV; ++i) {
//       J=nFunctionsU;
//       I=nFunctionsU + nFunctionsV;
//       for(int j = 0; j < nFunctionsL; ++j) {
//         Ae[I][J++] += _feV[j] * _feL[i] * jac * w[k];
//       }
//       I++;
//     }

//     for(int i = 0; i < nFunctionsL; ++i) {
//       J=nFunctionsU + nFunctionsV +nFunctionsL;
//       I=nFunctionsU + nFunctionsV;
//       for(int j = 0; j < nFunctionsU; ++j) {
//         Ae[I][J++] += -_feL[i] * _feBu[j] * jac * w[k];
//       }
//       I++;
//     }
//     for(int i = 0; i < nFunctionsL; ++i) {
//       J=nFunctionsU + nFunctionsV +nFunctionsL + nFunctionsBu;
//       I=nFunctionsU + nFunctionsV;
//       for(int j = 0; j < nFunctionsU; ++j) {
//         Ae[I][J++] += -_feL[i] * _feBv[j] * jac * w[k];
//       }
//       I++;
//     }

//     for(int i = 0; i < nFunctionsU; ++i) {
//       J=nFunctionsU + nFunctionsV +nFunctionsL;
//       I=nFunctionsU + nFunctionsV +nFunctionsL;
//       for(int j = 0; j < nFunctionsU; ++j) {
//         Ae[I][J++] +=
//           c0 * _feBu[i] * _feBv[j] * jac * w[k];
//       }
//       I++;
//     }
//     for(int i = 0; i < nFunctionsU; ++i) {
//       J=nFunctionsU + nFunctionsV +nFunctionsL + nFunctionsBu;
//       I=nFunctionsU + nFunctionsV +nFunctionsL + nFunctionsBu;
//       for(int j = 0; j < nFunctionsU; ++j) {
//         Ae[I][J++] +=
//           c0 * _feBu[j] * _feBv[i] * jac * w[k];
//       }
//       I++;
//     }
//   }
// }

// void feSysElm_1D_weakBC_Vec::computeBe(std::vector<double> &J, int numElem,
//                                         std::vector<feSpace *> &intSpace, feSpace *geoSpace,
//                                         std::vector<double> &geoCoord, double c0, double tn,
//                                         double dt, double *Be, std::vector<std::vector<double>> &sol,
//                                         std::vector<std::vector<double>> &solDot) {
//   int nG = geoSpace->getNbQuadPoints();
//   std::vector<double> w = geoSpace->getQuadratureWeights();
//   // Make the difference between different type of space
//   int nFunctionsU = intSpace[_idU]->getNbFunctions();
//   int nFunctionsV = intSpace[_idV]->getNbFunctions();
//   int nFunctionsL = intSpace[_idL]->getNbFunctions();
//   int nFunctionsBv = intSpace[_idBu]->getNbFunctions();
//   int nFunctionsBu = intSpace[_idBv]->getNbFunctions();

//   std::vector<double> x(3, 0.);
//   std::vector<double> gammaDot(2, 0.);
//   double u, l, v, bu,bv, buDot,bvDot, jac;
//   // std::cout<<"gammaDot= "<<gammaDot<< " x= "<<x[0]<< "tn= " <<tn<<std::endl;
//   // printf("%16.6e \n" , v);
//   for(int k = 0; k < nG; ++k) {
//     std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
//     geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
//     jac = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);
//     geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
//     if(_fct != nullptr) _fct->eval(tn, x, gammaDot);
//     u = intSpace[_idU]->interpolateFieldAtQuadNode(sol, k);
//     v = intSpace[_idV]->interpolateFieldAtQuadNode(sol, k);
//     l = intSpace[_idL]->interpolateFieldAtQuadNode(sol, k);
//     bv = intSpace[_idBu]->interpolateFieldAtQuadNode(sol, k);
//     bu = intSpace[_idBv]->interpolateFieldAtQuadNode(sol, k);
//     buDot = intSpace[_idBu]->interpolateFieldAtQuadNode(solDot, k);
//     bvDot = intSpace[_idBv]->interpolateFieldAtQuadNode(solDot, k);

//     int cnt = 0;
//     // std::cout<<"Cnt vaut  "<<cnt<<std::endl;
//     //Residu pour u
//     for(int i = 0; i < nFunctionsU; ++i) {
//       _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
//       Be[cnt++] -= _feU[i] * l * jac * w[k];
//     }
//     // std::cout<<"Cnt vaut  "<<cnt<<std::endl;
//     //Residu pour v
//     for(int i = 0; i < nFunctionsV; ++i) {
//       _feV[i] = intSpace[_idV]->getFunctionAtQuadNode(i, k);
//       Be[cnt++] -= _feV[i] * l * jac * w[k];
//     }
//     // std::cout<<"Cnt vaut  "<<cnt<<std::endl;
//     //Residu pour lambda
//     for(int i = 0; i < nFunctionsL; ++i) {
//       _feL[i] = intSpace[_idL]->getFunctionAtQuadNode(i, k);
//       Be[cnt++] -= _feL[i] * (u - bu) * jac * w[k];
//     }
//     // std::cout<<"Cnt vaut  "<<cnt<<std::endl;
//     for(int i = 0; i < nFunctionsL; ++i) {
//       _feL[i] = intSpace[_idL]->getFunctionAtQuadNode(i, k);
//       Be[cnt++] -= _feL[i] * (v - bv) * jac * w[k];
//     }
//     // std::cout<<"Cnt vaut  "<<cnt<<std::endl;
//     //Residu du a l'edo 1
//     for(int i = 0; i < nFunctionsBu; ++i) {
//       _feBu[i] = intSpace[_idBu]->getFunctionAtQuadNode(i, k);
//       Be[cnt++] -= _feBu[i] * (buDot - gammaDot[0]) * jac * w[k];
//     }
//     // std::cout<<"Cnt vaut  "<<cnt<<std::endl;
//     for(int i = 0; i < nFunctionsBv; ++i) {
//       _feBv[i] = intSpace[_idBv]->getFunctionAtQuadNode(i, k);
//       Be[cnt++] -= _feBv[i] * (bvDot - gammaDot[1]) * jac * w[k];
//     }
//   }
// }

void feSysElm_1D_Source::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
}

void feSysElm_1D_Source::computeBe(std::vector<double> &J, int numElem,
                                   std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                   std::vector<double> &geoCoord, double c0, double tn, double dt,
                                   double *Be, std::vector<std::vector<double>> &sol,
                                   std::vector<std::vector<double>> &solDot,
                                   std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
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
      // std::cout<<"feU["<<i<<"] = "<<_feU[i]<<std::endl;
      // std::cout<<"J = "<<J[nG * numElem + k]<<std::endl;
      // std::cout<<"w["<<k<<"] = "<<w[k]<<std::endl;
      // std::cout<<"Be["<<i<<"] = "<<Be[i]<<std::endl;
    }
  }
}

void feSysElm_1D_Diffusion::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
}

void feSysElm_1D_Diffusion::computeAe(std::vector<double> &J, int numElem,
                                      std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                      std::vector<double> &geoCoord, double c0, double tn,
                                      double **Ae, std::vector<std::vector<double>> &sol,
                                      std::vector<std::vector<double>> &solDot)
{
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
      for(int j = 0; j < nFunctions; ++j) {
        Ae[i][j] += _feUdx[i] * kD * _feUdx[j] * jac * w[k];
      }
    }
  }
}

void feSysElm_1D_Diffusion::computeBe(std::vector<double> &J, int numElem,
                                      std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                      std::vector<double> &geoCoord, double c0, double tn,
                                      double dt, double *Be, std::vector<std::vector<double>> &sol,
                                      std::vector<std::vector<double>> &solDot,
                                      std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double kD = _par;
  int nFunctions = intSpace[_idU]->getNbFunctions();

  double jac, dudx;
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * numElem + k];

    dudx = intSpace[_idU]->interpolateFieldAtQuadNode_rDerivative(sol[_idU], k);
    dudx /= jac;

    for(int i = 0; i < nFunctions; ++i) {
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k);
      _feUdx[i] /= jac;
      Be[i] -= _feUdx[i] * kD * dudx * jac * w[k];
    }
  }
}

void feSysElm_1D_Masse::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
}

void feSysElm_1D_Masse::computeAe(std::vector<double> &J, int numElem,
                                  std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                  std::vector<double> &geoCoord, double c0, double tn, double **Ae,
                                  std::vector<std::vector<double>> &sol, std::vector<std::vector<double>> &solDot)
{
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double rho = _par;
  int nFunctions = intSpace[_idU]->getNbFunctions();

  double jac;
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * numElem + k];

    for(int i = 0; i < nFunctions; ++i)
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);

    for(int i = 0; i < nFunctions; ++i)
      for(int j = 0; j < nFunctions; ++j)
        Ae[i][j] += _feU[i] * rho * c0 * _feU[j] * jac * w[k];
  }
}

void feSysElm_1D_Masse::computeBe(std::vector<double> &J, int numElem,
                                  std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                  std::vector<double> &geoCoord, double c0, double tn, double dt,
                                  double *Be, std::vector<std::vector<double>> &sol, std::vector<std::vector<double>> &solDot,
                                  std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double rho = _par;
  int nFunctions = intSpace[_idU]->getNbFunctions();

  double jac, uDot;
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * numElem + k];

    uDot = intSpace[_idU]->interpolateFieldAtQuadNode(solDot[_idU], k);

    for(int i = 0; i < nFunctions; ++i) {
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
      Be[i] -= _feU[i] * rho * uDot * jac * w[k];
    }
  }
}

void feSysElm_1D_NeumannBC::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
}

void feSysElm_1D_NeumannBC::computeBe(std::vector<double> &J, int numElem,
                                      std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                      std::vector<double> &geoCoord, double c0, double tn,
                                      double dt, double *Be, std::vector<std::vector<double>> &sol,
                                      std::vector<std::vector<double>> &solDot,
                                      std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  int nFunctions = intSpace[_idU]->getNbFunctions();

  double jac;
  std::vector<double> x(3, 0.0);
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * numElem + k];
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
    double h = _fct->eval(tn, x);
    for(int i = 0; i < nFunctions; ++i) {
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
      Be[i] += _feU[i] * h * jac * w[k];
    }
  }
}

void feSysElm_1D_Advection::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
}

void feSysElm_1D_Advection::computeAe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNbQuadPoints();
  int nFunctions = form->_intSpace[_idU]->getNbFunctions();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  std::vector<double> &J = form->_cnc->getJacobians();
  std::vector<double> x(3, 0.0);

  double jac;
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * form->_numElem + k];

    for(int i = 0; i < nFunctions; ++i){
      _feU[i] = form->_intSpace[_idU]->getFunctionAtQuadNode(i, k);
      _feUdx[i] = form->_intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) / jac;
    }

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, x);
    double cVelocity = _fct->eval(form->_tn, x);

    for(int i = 0; i < nFunctions; ++i) {
      for(int j = 0; j < nFunctions; ++j) {
        form->_Ae[i][j] -= _feUdx[i] * cVelocity * _feU[j] * jac * w[k];
      }
    }
  }
}

void feSysElm_1D_Advection::computeBe(std::vector<double> &J, int numElem,
                                      std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                      std::vector<double> &geoCoord, double c0, double tn,
                                      double dt, double *Be,
                                      std::vector<std::vector<double>> &sol,
                                      std::vector<std::vector<double>> &solDot,
                                      std::vector<std::vector<double>> &solPrev,
                                      std::vector<std::vector<double>> &solNext)
{
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  int nFunctions = intSpace[_idU]->getNbFunctions();

  double u, dudx, jac, cVelocity;
  std::vector<double> x(3, 0.0);
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * numElem + k];
    u = intSpace[_idU]->interpolateFieldAtQuadNode(sol[_idU], k);
    dudx = intSpace[_idU]->interpolateFieldAtQuadNode_rDerivative(sol[_idU], k) / jac;

    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
    cVelocity = _fct->eval(tn, x);
    for(int i = 0; i < nFunctions; ++i) {
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k);
      _feUdx[i] /= jac;
      Be[i] += cVelocity * u * _feUdx[i] * jac * w[k];
    }
  }
}

void feSysElm_1D_SUPGStab::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
}

void feSysElm_1D_SUPGStab::computeAe(std::vector<double> &J, int numElem,
                                      std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                      std::vector<double> &geoCoord, double c0, double tn,
                                      double **Ae, std::vector<std::vector<double>> &sol,
                                      std::vector<std::vector<double>> &solDot)
{
  // int nG = geoSpace->getNbQuadPoints();
  // int nFunctions = intSpace[_idU]->getNbFunctions();
  // std::vector<double> w = geoSpace->getQuadratureWeights();
  // std::vector<double> x(3, 0.0);

  // double jac, cVelocity;
  // for(int k = 0; k < nG; ++k) {
  //   jac = J[nG * numElem + k];

  //   for(int i = 0; i < nFunctions; ++i){
  //     _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
  //     _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) / jac;
  //   }

  //   geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
  //   cVelocity = _fct->eval(tn, x);

  //   for(int i = 0; i < nFunctions; ++i) {
  //     for(int j = 0; j < nFunctions; ++j) {
  //       Ae[i][j] -= _feUdx[i] * cVelocity * _feU[j] * jac * w[k];
  //     }
  //   }
  // }
}

void feSysElm_1D_SUPGStab::computeBe(std::vector<double> &J, int numElem,
                                      std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                      std::vector<double> &geoCoord, double c0, double tn,
                                      double dt, double *Be,
                                      std::vector<std::vector<double>> &sol,
                                      std::vector<std::vector<double>> &solDot,
                                      std::vector<std::vector<double>> &solPrev,
                                      std::vector<std::vector<double>> &solNext)
{
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  int nFunctions = intSpace[_idU]->getNbFunctions();

  double kDiffusivity = _par[0];

  double u, dudx, jac, cVelocity;
  std::vector<double> x(3, 0.0);
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * numElem + k];
    u = intSpace[_idU]->interpolateFieldAtQuadNode(sol[_idU], k);
    dudx = intSpace[_idU]->interpolateFieldAtQuadNode_rDerivative(sol[_idU], k) / jac;

    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
    cVelocity = _fct->eval(tn, x);
    for(int i = 0; i < nFunctions; ++i) {
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k);
      _feUdx[i] /= jac;

      // Compute element Peclet and optimal numerical diffusivity coth(Peh) - 1/Peh
      double he = 2.0 * jac;
      double Peh = cVelocity * he / (2. * kDiffusivity);
      double beta = cosh(Peh)/sinh(Peh) - (1./Peh);
      double tau = beta * he / (2. * cVelocity);
      double residual = cVelocity * dudx - 1.; // FIXME : completer avec d2dx2
      Be[i] -= cVelocity * _feUdx[i] * tau * residual * jac * w[k];
    }
  }
}

void feSysElm_1D_DG_Advection::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
}

void feSysElm_1D_DG_Advection::computeAe(std::vector<double> &J, int numElem,
                                      std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                      std::vector<double> &geoCoord, double c0, double tn,
                                      double **Ae, std::vector<std::vector<double>> &sol,
                                      std::vector<std::vector<double>> &solDot)
{
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  double jac;
  int nFunctions = intSpace[_idU]->getNbFunctions();
  std::vector<double> x(3, 0.0);

  // Element term: integral of c * dphi_i/dx * phi_j
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * numElem + k];

    // Initialize vectors phi (_feU) and dphi/dx = dphi/dr * dr/dx (1/J) (_feUdx)
    for(int i = 0; i < nFunctions; ++i){
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) / jac;
    }

    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
    double cVelocity = _fct->eval(tn, x);

    for(int i = 0; i < nFunctions; ++i) {
      for(int j = 0; j < nFunctions; ++j) {
        Ae[i][j] -= _feUdx[i] * cVelocity * _feU[j] * jac * w[k];
      }
    }
  }
}

void feSysElm_1D_DG_Advection::computeBe(std::vector<double> &J, int numElem,
                                      std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                      std::vector<double> &geoCoord, double c0, double tn,
                                      double dt, double *Be,
                                      std::vector<std::vector<double>> &sol,
                                      std::vector<std::vector<double>> &solDot,
                                      std::vector<std::vector<double>> &solPrev,
                                      std::vector<std::vector<double>> &solNext)
{
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  int nFunctions = intSpace[_idU]->getNbFunctions();

  double u, jac, cVelocity;
  std::vector<double> x(3, 0.0);
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * numElem + k];
    u = intSpace[_idU]->interpolateFieldAtQuadNode(sol[_idU], k);
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
    cVelocity = _fct->eval(tn, x);
    for(int i = 0; i < nFunctions; ++i) {
      _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k);
      _feUdx[i] /= jac;
      Be[i] += cVelocity * u * _feUdx[i] * jac * w[k];
      // feInfo("residual = %f", cVelocity * u * _feUdx[i] * jac * w[k]);
    }
  }

  double rLeft[3]  = {-1.0, 0., 0.};
  double rRight[3] = {+1.0, 0., 0.};
  std::vector<double> phiL(nFunctions);
  std::vector<double> phiR(nFunctions);
  intSpace[_idU]->L(rLeft,  phiL.data());
  intSpace[_idU]->L(rRight, phiR.data());

  double uMinus, uPlus;
  double uL, uR, uPrevR, uNextL;
  uL = intSpace[_idU]->interpolateField(sol[_idU], rLeft);
  uR = intSpace[_idU]->interpolateField(sol[_idU], rRight);
  uPrevR = intSpace[_idU]->interpolateField(solPrev[_idU], rRight);
  uNextL = intSpace[_idU]->interpolateField(solNext[_idU], rLeft);

  // Moyenne:
  uPlus  = (uR + uNextL)/2.;
  uMinus = (uL + uPrevR)/2.;

  // Interface (flux) term: - [(c*uh)*phi_i]_X^e-1 ^X^e
  for(int i = 0; i < nFunctions; ++i) {
    // if(cVelocity > 0){
    //   // Information comes from the left:
    //   // Interpolate solution on the right of the previous element
    //   uMinus = intSpace[_idU]->interpolateField(solPrev[_idU], rRight);
    //   uPlus  = intSpace[_idU]->interpolateField(    sol[_idU], rRight);
    // } else{
    //   // Information comes from the right:
    //   // Interpolate solution on the left of the next element
    //   uMinus = intSpace[_idU]->interpolateField(    sol[_idU], rLeft);
    //   uPlus  = intSpace[_idU]->interpolateField(solNext[_idU], rLeft);
    // }

    Be[i] -= (cVelocity * uPlus) * phiR[i] - (cVelocity * uMinus) * phiL[i];

    // feInfo("Elem %d - u+ = %f - u- = %f",
    //   uPlus,
    //   uMinus);

    // feInfo("Elem %d - Adding flux %d = %f - %f",
    //   numElem,
    //   i,
    //   (cVelocity * uPlus) * phiR[i] - (cVelocity * uMinus) * phiL[i],
    //   uPlus - uMinus);
  }
}

void feSysElm_2D_Source::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
  _feUdy.resize(space[_idU]->getNbFunctions());
  _x.resize(3);
}

void feSysElm_2D_Source::computeBe(std::vector<double> &J, int numElem,
                                   std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                   std::vector<double> &geoCoord, double c0, double tn, double dt,
                                   double *Be, std::vector<std::vector<double>> &sol,
                                   std::vector<std::vector<double>> &solDot,
                                   std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> &w = geoSpace->getQuadratureWeights();
  int nFunctions = intSpace[_idU]->getNbFunctions();
  bool globalFunctions = intSpace[_idU]->useGlobalFunctions();

  double jac;
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * numElem + k];
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, _x);

    double S = _fct->eval(tn, _x);

    for(int i = 0; i < nFunctions; ++i) {
      if(globalFunctions) {
        _feU[i] = intSpace[_idU]->getGlobalFunctionAtQuadNode(numElem, i, k);
      } else {
        _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
      }
      Be[i] -= _feU[i] * S * jac * w[k];
      // #pragma omp critical
      // {
      // printf("%d;%d;%d;",numElem,i,k);std::cout<<Be[i]<<std::endl;
      // }
    }
  }
}

void feSysElm_2D_Diffusion::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
  _feUdy.resize(space[_idU]->getNbFunctions());
  _dxdr.resize(3);
  _dxds.resize(3);
}

void feSysElm_2D_Diffusion::computeAe(std::vector<double> &Ja, int numElem,
                                      std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                      std::vector<double> &geoCoord, double c0, double tn,
                                      double **Ae, std::vector<std::vector<double>> &sol,
                                      std::vector<std::vector<double>> &solDot)
{
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> &w = geoSpace->getQuadratureWeights();
  double kD = _par;
  int nFunctions = intSpace[_idU]->getNbFunctions();
  // bool globalFunctions = intSpace[_idU]->useGlobalFunctions();

  double J;
  for(int k = 0; k < nG; ++k) {
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, _dxdr);
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, _dxds);

    // if(!globalFunctions) {
      // J = _dxdr[0]*_dxds[1]-_dxdr[1]*_dxds[0];
      J = Ja[nG * numElem + k];
      double drdx = _dxds[1] / J;
      double drdy = -_dxds[0] / J;
      double dsdx = -_dxdr[1] / J;
      double dsdy = _dxdr[0] / J;

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

    // } else {
    //   // Using global interpolation functions
    //   for(int i = 0; i < nFunctions; ++i) {
    //     _feUdx[i] = intSpace[_idU]->getdGlobalFunctiondxAtQuadNode(numElem, i, k);
    //     _feUdy[i] = intSpace[_idU]->getdGlobalFunctiondyAtQuadNode(numElem, i, k);
    //   }

    //   for(int i = 0; i < nFunctions; ++i) {
    //     for(int j = 0; j < nFunctions; ++j) {
    //       Ae[i][j] += (_feUdx[i] * _feUdx[j] + _feUdy[i] * _feUdy[j]) * kD * J * w[k];
    //     }
    //   }
    // }
  }

  // if(false) {
  //   // 1D quadrature rule to integrate over edges
  //   feQuadrature rule(20, 1, "");
  //   std::vector<double> x1D = rule.getXPoints();
  //   std::vector<double> w1D = rule.getWeights();
  //   std::vector<double> r1D(x1D.size());
  //   int n1D = rule.getNQuad();
  //   // 1D quad nodes from [-1,1] to [0,1]
  //   for(int i = 0; i < n1D; ++i) {
  //     r1D[i] = (x1D[i] + 1.0) / 2.0;
  //   }

  //   double drdt[3] = {1.0, -1.0, 0.0};
  //   double dsdt[3] = {0.0, 1.0, -1.0};

  //   // Il faut recalculer les fonctions globales aux noeuds de quadrature 1D
  //   // Pas tres efficace : elles ont ete precalculees mais seulement aux noeuds 2D
  //   // TODO : stocker les matrices de coeff
  //   std::vector<double> l(nFunctions, 0.0);
  //   std::vector<double> dldx(nFunctions, 0.0);
  //   std::vector<double> dldy(nFunctions, 0.0);

  //   std::string name = "elem" + std::to_string(numElem) + ".pos";
  //   FILE *f = fopen(name.c_str(), "w");
  //   fprintf(f, "View \"elem\" {\n");

  //   // Integrale sur les bords (devra devenir une autre forme faible si ça fonctionne)
  //   for(int iEdge = 0; iEdge < 3; ++iEdge) {
  //     double longueur = 0.0;
  //     double longueur2 = 0.0;
  //     for(int k = 0; k < n1D; ++k) {
  //       std::vector<double> x(3, 0.0);
  //       std::vector<double> _dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
  //       std::vector<double> _dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]

  //       double jac1D;

  //       switch(iEdge) {
  //         case 0: {
  //           double r[3] = {r1D[k], 0., 0.};
  //           // These interpolations use local functions because geometry is defined with Lagrange
  //           // local polynomials
  //           geoSpace->interpolateVectorField(geoCoord, r, x);
  //           geoSpace->interpolateVectorField_rDerivative(geoCoord, r, _dxdr);
  //           geoSpace->interpolateVectorField_sDerivative(geoCoord, r, _dxds);

  //           double tx = _dxdr[0] * drdt[iEdge] + _dxds[0] * dsdt[iEdge];
  //           double ty = _dxdr[1] * drdt[iEdge] + _dxds[1] * dsdt[iEdge];

  //           // printf("Tangential vector on edge 12 at phys point (%+-4.4f, %+-4.4f) - ref point
  //           // (%4.4f, %4.4f) =
  //           // (%+-4.4f, %+-4.4f)\n", x[0], x[1], r[0], r[1], tx, ty);
  //           // fprintf(f,"VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n",x[0],x[1],0.,tx,ty,0.0);
  //           // fprintf(f,"VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n",x[0],x[1],0.,ty,-tx,0.0);
  //           break;
  //         }
  //         case 1: {
  //           double r[3] = {r1D[k], 1.0 - r1D[k], 0.};
  //           geoSpace->interpolateVectorField(geoCoord, r, x);
  //           geoSpace->interpolateVectorField_rDerivative(geoCoord, r, _dxdr);
  //           geoSpace->interpolateVectorField_sDerivative(geoCoord, r, _dxds);
  //           double tx = _dxdr[0] * drdt[iEdge] + _dxds[0] * dsdt[iEdge];
  //           double ty = _dxdr[1] * drdt[iEdge] + _dxds[1] * dsdt[iEdge];

  //           // printf("Tangential vector on edge 23 at phys point (%+-4.4f, %+-4.4f) - ref point
  //           // (%4.4f, %4.4f) =
  //           // (%+-4.4f, %+-4.4f)\n",
  //           //   x[0], x[1], r[0], r[1], tx, ty);
  //           // fprintf(f,"VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n",x[0],x[1],0.,tx,ty,0.0);
  //           // fprintf(f,"VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n",x[0],x[1],0.,ty,-tx,0.0);
  //           break;
  //         }
  //         case 2: {
  //           double r[3] = {0.0, r1D[k], 0.};
  //           geoSpace->interpolateVectorField(geoCoord, r, x);
  //           geoSpace->interpolateVectorField_rDerivative(geoCoord, r, _dxdr);
  //           geoSpace->interpolateVectorField_sDerivative(geoCoord, r, _dxds);
  //           double tx = _dxdr[0] * drdt[iEdge] + _dxds[0] * dsdt[iEdge];
  //           double ty = _dxdr[1] * drdt[iEdge] + _dxds[1] * dsdt[iEdge];

  //           // printf("Tangential vector on edge 31 at phys point (%+-4.4f, %+-4.4f) - ref point
  //           // (%4.4f, %4.4f) =
  //           // (%+-4.4f, %+-4.4f)\n",
  //           //   x[0], x[1], r[0], r[1], tx, ty);
  //           // fprintf(f,"VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n",x[0],x[1],0.,tx,ty,0.0);
  //           // fprintf(f,"VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n",x[0],x[1],0.,ty,-tx,0.0);
  //           break;
  //         }
  //       }

  //       double dxdt = _dxdr[0] * drdt[iEdge] + _dxds[0] * dsdt[iEdge];
  //       double dydt = _dxdr[1] * drdt[iEdge] + _dxds[1] * dsdt[iEdge];

  //       // Le 1/2 vient de la transformation de [-1,1] à [0,1], qui est ensuite envoyé dans l'espace
  //       // physique
  //       jac1D = sqrt(dxdt * dxdt + dydt * dydt) / 2.0;

  //       longueur += jac1D * w1D[k];
  //       longueur2 += J * w1D[k];

  //       double nx = dydt;
  //       double ny = -dxdt;
  //       double N = 2.0 * jac1D;
  //       nx /= N;
  //       ny /= N;
  //       fprintf(f, "VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n", x[0], x[1], 0., -ny, nx, 0.0);
  //       fprintf(f, "VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n", x[0], x[1], 0., nx, ny, 0.0);

  //       intSpace[_idU]->Lphys(numElem, x, l, dldx, dldy);

  //       for(int i = 0; i < nFunctions; ++i) {
  //         _feU[i] = l[i];
  //         _feUdx[i] = dldx[i];
  //         _feUdy[i] = dldy[i];
  //       }

  //       for(int i = 0; i < nFunctions; ++i) {
  //         for(int j = 0; j < nFunctions; ++j) {
  //           Ae[i][j] -= (_feU[i] * (_feUdx[j] * nx + _feUdy[j] * ny)) * kD * jac1D * w1D[k];
  //         }
  //       }
  //     }
  //     // printf("Longueur = %4.4f et %4.4f\n", longueur, longueur2);
  //   }

  //   fprintf(f, "};");
  //   fclose(f);
  // }
}

void feSysElm_2D_Diffusion::computeBe(std::vector<double> &Ja, int numElem,
                                      std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                      std::vector<double> &geoCoord, double c0, double tn,
                                      double dt, double *Be, std::vector<std::vector<double>> &sol,
                                      std::vector<std::vector<double>> &solDot,
                                      std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> &w = geoSpace->getQuadratureWeights();
  double kD = _par;
  int nFunctions = intSpace[_idU]->getNbFunctions();
  bool globalFunctions = intSpace[_idU]->useGlobalFunctions();

  double J, dudx, dudy;
  for(int k = 0; k < nG; ++k) {

    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, _dxdr);
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, _dxds);
    J = _dxdr[0]*_dxds[1]-_dxdr[1]*_dxds[0];

    if(!globalFunctions) {
      double drdx = _dxds[1] / J;
      double drdy = -_dxds[0] / J;
      double dsdx = -_dxdr[1] / J;
      double dsdy = _dxdr[0] / J;

      dudx = intSpace[_idU]->interpolateFieldAtQuadNode_rDerivative(sol[_idU], k) * drdx +
             intSpace[_idU]->interpolateFieldAtQuadNode_sDerivative(sol[_idU], k) * dsdx;
      dudy = intSpace[_idU]->interpolateFieldAtQuadNode_rDerivative(sol[_idU], k) * drdy +
             intSpace[_idU]->interpolateFieldAtQuadNode_sDerivative(sol[_idU], k) * dsdy;

      for(int i = 0; i < nFunctions; ++i) {
        _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                    intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
        _feUdy[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                    intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
        Be[i] -= (_feUdx[i] * dudx + _feUdy[i] * dudy) * kD * J * w[k];
      }
    } else {

      feWarning("FORME FAIBLE DOIT ETRE CORRIGEE POUR LES FONCTIONS GLOBALES");
      feWarning("FORME FAIBLE DOIT ETRE CORRIGEE POUR LES FONCTIONS GLOBALES");
      feWarning("FORME FAIBLE DOIT ETRE CORRIGEE POUR LES FONCTIONS GLOBALES");
      feWarning("FORME FAIBLE DOIT ETRE CORRIGEE POUR LES FONCTIONS GLOBALES");
    //   dudx = intSpace[_idU]->interpolateFieldAtQuadNode_xDerivative(sol, numElem, k);
    //   dudy = intSpace[_idU]->interpolateFieldAtQuadNode_yDerivative(sol, numElem, k);

    //   for(int i = 0; i < nFunctions; ++i) {
    //     _feUdx[i] = intSpace[_idU]->getdGlobalFunctiondxAtQuadNode(numElem, i, k);
    //     _feUdy[i] = intSpace[_idU]->getdGlobalFunctiondyAtQuadNode(numElem, i, k);
    //     Be[i] -= (_feUdx[i] * dudx + _feUdy[i] * dudy) * kD * J * w[k];
      // }
    }
  }

  // if(false) {
  //   // 1D quadrature rule to integrate over edges
  //   feQuadrature rule(30, 1, "");
  //   std::vector<double> x1D = rule.getXPoints();
  //   std::vector<double> w1D = rule.getWeights();
  //   std::vector<double> r1D(x1D.size());
  //   int n1D = rule.getNQuad();
  //   // 1D quad nodes from [-1,1] to [0,1]
  //   for(int i = 0; i < n1D; ++i) {
  //     r1D[i] = (x1D[i] + 1.0) / 2.0;
  //   }

  //   double drdt[3] = {1.0, -1.0, 0.0};
  //   double dsdt[3] = {0.0, 1.0, -1.0};

  //   // std::string name = "elem" + std::to_string(numElem) + ".pos";
  //   // FILE *f = fopen(name.c_str(), "w");
  //   // fprintf(f, "View \"elem\" {\n");

  //   // Il faut recalculer les fonctions globales aux noeuds de quadrature 1D
  //   // Pas tres efficace : elles ont ete precalculees mais seulement aux noeuds 2D
  //   // TODO : stocker les matrices de coeff
  //   std::vector<double> l(nFunctions, 0.0);
  //   std::vector<double> dldx(nFunctions, 0.0);
  //   std::vector<double> dldy(nFunctions, 0.0);

  //   // std::vector<double> xc(3, 0.0);
  //   // double rc[3] = {1. / 3., 1. / 3., 1. / 3.};
  //   // geoSpace->interpolateVectorField(geoCoord, rc, xc);

  //   // Integrale sur les bords (devra devenir une autre forme faible si ça fonctionne)
  //   for(int iEdge = 0; iEdge < 3; ++iEdge) {
  //     double longueur = 0.0;
  //     double longueur2 = 0.0;

  //     double intEdge = 0.0;

  //     for(int k = 0; k < n1D; ++k) {
  //       std::vector<double> x(3, 0.0);
  //       // std::vector<double> _dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
  //       // std::vector<double> _dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]

  //       double jac1D;

  //       switch(iEdge) {
  //         case 0: {
  //           double r[3] = {r1D[k], 0., 0.};
  //           // These interpolations use local functions because geometry is defined with Lagrange
  //           // local polynomials
  //           geoSpace->interpolateVectorField(geoCoord, r, x);
  //           geoSpace->interpolateVectorField_rDerivative(geoCoord, r, _dxdr);
  //           geoSpace->interpolateVectorField_sDerivative(geoCoord, r, _dxds);
  //           break;
  //         }
  //         case 1: {
  //           double r[3] = {1.0 - r1D[k], r1D[k], 0.};
  //           geoSpace->interpolateVectorField(geoCoord, r, x);
  //           geoSpace->interpolateVectorField_rDerivative(geoCoord, r, _dxdr);
  //           geoSpace->interpolateVectorField_sDerivative(geoCoord, r, _dxds);
  //           break;
  //         }
  //         case 2: {
  //           double r[3] = {0.0, 1.0 - r1D[k], 0.};
  //           geoSpace->interpolateVectorField(geoCoord, r, x);
  //           geoSpace->interpolateVectorField_rDerivative(geoCoord, r, _dxdr);
  //           geoSpace->interpolateVectorField_sDerivative(geoCoord, r, _dxds);
  //           break;
  //         }
  //       }

  //       double dxdt = _dxdr[0] * drdt[iEdge] + _dxds[0] * dsdt[iEdge];
  //       double dydt = _dxdr[1] * drdt[iEdge] + _dxds[1] * dsdt[iEdge];

  //       // Le 1/2 vient de la transformation de [-1,1] à [0,1], qui est ensuite envoyé dans l'espace
  //       // physique
  //       jac1D = sqrt(dxdt * dxdt + dydt * dydt);

  //       longueur += jac1D * w1D[k];
  //       // longueur2 += J * w1D[k];

  //       double nx = dydt;
  //       double ny = -dxdt;
  //       double N = sqrt(dxdt * dxdt + dydt * dydt);
  //       nx /= N;
  //       ny /= N;
  //       // fprintf(f, "VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n", x[0], x[1], 0., -ny, nx, 0.0);
  //       // fprintf(f, "VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n", x[0], x[1], 0., nx, ny, 0.0);

  //       // x[0] -= xc[0];
  //       // x[1] -= xc[1];
  //       // x[2] -= xc[2];

  //       intSpace[_idU]->Lphys(numElem, x, l, dldx, dldy);

  //       dudx = intSpace[_idU]->interpolateField_xDerivative(sol, numElem, x);
  //       dudy = intSpace[_idU]->interpolateField_yDerivative(sol, numElem, x);

  //       // double sumPhi = 0.0;
  //       for(int i = 0; i < nFunctions; ++i) {
  //         _feU[i] = l[i];
  //         Be[i] += _feU[i] * (dudx * nx + dudy * ny) * kD * jac1D * w1D[k];
  //         intEdge += (_feU[i] * (dudx * nx + dudy * ny)) * kD * jac1D * w1D[k];
  //         // sumPhi += _feU[i];
  //         // sumPhi += dldx[i];
  //         // sumPhi += dldy[i];
  //         // std::cout<<_feU[i]<<" - "<<dudx<<" - "<<dudy<<" - "<<nx<<" - "<<ny<<" - "<<jac1D<<" -
  //         // "<<w1D[k]<<std::endl;
  //       }
  //       // printf("Sumphi = %+-10.10e\n", sumPhi);
  //     }
  //     // printf("Integrale sur l'arete %d de l'elm %d = %+-10.10e\n", iEdge, numElem, intEdge);
  //     // printf("Longueur = %10.16e\n", longueur);
  //   }
  //   // fprintf(f, "};");
  //   // fclose(f);
  // }
}

void feSysElm_2D_Masse::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
  _feUdy.resize(space[_idU]->getNbFunctions());
  _dxdr.resize(3);
  _dxds.resize(3);
}

void feSysElm_2D_Masse::computeAe(std::vector<double> &J, int numElem,
                                  std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                  std::vector<double> &geoCoord, double c0, double tn, double **Ae,
                                  std::vector<std::vector<double>> &sol, std::vector<std::vector<double>> &solDot)
{
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> &w = geoSpace->getQuadratureWeights();
  double rho = _par;
  int nFunctionsU = intSpace[_idU]->getNbFunctions();

  double jac, u, dudt;
  for(int k = 0; k < nG; ++k) {
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, _dxdr);
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, _dxds);
    jac = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    u = intSpace[_idU]->interpolateFieldAtQuadNode(sol[_idU], k);
    dudt = intSpace[_idU]->interpolateFieldAtQuadNode(solDot[_idU], k);

    for(int i = 0; i < nFunctionsU; ++i) _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);

    for(int i = 0; i < nFunctionsU; ++i) {
      for(int j = 0; j < nFunctionsU; ++j) {
        Ae[i][j] += _feU[i] * rho * c0 * _feU[j] * jac * w[k];
      }
    }
  }
}

void feSysElm_2D_Masse::computeBe(std::vector<double> &J, int numElem,
                                  std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                  std::vector<double> &geoCoord, double c0, double tn, double dt,
                                  double *Be, std::vector<std::vector<double>> &sol, std::vector<std::vector<double>> &solDot,
                                  std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> &w = geoSpace->getQuadratureWeights();
  double rho = _par;
  int nFunctionsU = intSpace[_idU]->getNbFunctions();

  double jac, dudt;
  for(int k = 0; k < nG; ++k) {
    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, _dxdr);
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, _dxds);
    jac = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    dudt = intSpace[_idU]->interpolateFieldAtQuadNode(solDot[_idU], k);

    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = intSpace[_idU]->getFunctionAtQuadNode(i, k);
      Be[i] -= _feU[i] * rho * dudt * jac * w[k];
    }
  }
}

void feSysElm_2D_Advection::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
  _feUdy.resize(space[_idU]->getNbFunctions());
  _dxdr.resize(3);
  _dxds.resize(3);
}

void feSysElm_2D_Advection::computeAe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double **Ae, std::vector<std::vector<double>> &sol, std::vector<std::vector<double>> &solDot)
{
  // Calculee par differences finies

  // int nG = geoSpace->getNbQuadPoints();
  // std::vector<double> w = geoSpace->getQuadratureWeights();
  // double kD = _par;
  // int nFunctions = intSpace[_idU]->getNbFunctions();

  // double J;
  // for(int k = 0; k < nG; ++k) {
  //   std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
  //   std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
  //   geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
  //   geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
  //   J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];

  //   double drdx = dxds[1] / J;
  //   double drdy = -dxds[0] / J;
  //   double dsdx = -dxdr[1] / J;
  //   double dsdy = dxdr[0] / J;

  //   for(int i = 0; i < nFunctions; ++i) {
  //     _feUdx[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
  //                 intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
  //     _feUdy[i] = intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
  //                 intSpace[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
  //   }

  //   for(int i = 0; i < nFunctions; ++i) {
  //     for(int j = 0; j < nFunctions; ++j) {
  //       Ae[i][j] += (_feUdx[i] * _feUdx[j] + _feUdy[i] * _feUdy[j]) * kD * J * w[k];
  //     }
  //   }
  // }
}

void feSysElm_2D_Advection::computeBe(std::vector<double> &J, int numElem, std::vector<feSpace *> &intSpace,
                         feSpace *geoSpace, std::vector<double> &geoCoord, double c0, double tn,
                         double dt, double *Be, std::vector<std::vector<double>> &sol,
                         std::vector<std::vector<double>> &solDot,
                         std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> w = geoSpace->getQuadratureWeights();
  int nFunctions = intSpace[_idU]->getNbFunctions();

  double Jac, u, dudt, dudx, dudy;
  for(int k = 0; k < nG; ++k) {

    // Evaluate the exterior velocity field
    std::vector<double> v(2, 0.0);
    std::vector<double> x(3, 0.0);
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
    _fct->eval(tn, x, v);

    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, _dxdr);
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, _dxds);
    Jac = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    double drdx = _dxds[1] / Jac;
    double drdy = -_dxds[0] / Jac;
    double dsdx = -_dxdr[1] / Jac;
    double dsdy = _dxdr[0] / Jac;

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
      u = intSpace[_idU]->interpolateFieldAtQuadNode(sol[_idU], k);
      dudt = intSpace[_idU]->interpolateFieldAtQuadNode(solDot[_idU],k);
      dudx = intSpace[_idU]->interpolateFieldAtQuadNode_rDerivative(sol[_idU], k) * drdx +
             intSpace[_idU]->interpolateFieldAtQuadNode_sDerivative(sol[_idU], k) * dsdx;
      dudy = intSpace[_idU]->interpolateFieldAtQuadNode_rDerivative(sol[_idU], k) * drdy +
             intSpace[_idU]->interpolateFieldAtQuadNode_sDerivative(sol[_idU], k) * dsdy;

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
      // Be[i] -= Jac * w[k] * ((v[0] * dudx + v[1] * dudy) * _feU[i] + deltaSUPG * (v[0] * _feUdx[i] + v[1] * _feUdy[i]) * (dudt + v[0] * dudx + v[1] * dudy));
      Be[i] -= Jac * w[k] * ((v[0] * dudx + v[1] * dudy) * _feU[i] + tau       * (v[0] * _feUdx[i] + v[1] * _feUdy[i]) * (dudt + v[0] * dudx + v[1] * dudy));
    }
  }
}

void feSysElm_2D_Stokes::createElementarySystem(std::vector<feSpace *> &space)
{
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
  _dxdr.resize(3); // [dx/dr, dy/dr, dz/dr]
  _dxds.resize(3); // [dx/ds, dy/ds, dz/ds]
}

void feSysElm_2D_Stokes::computeBe(std::vector<double> &Ja, int numElem,
                                   std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                   std::vector<double> &geoCoord, double c0, double tn, double dt,
                                   double *Be, std::vector<std::vector<double>> &sol,
                                   std::vector<std::vector<double>> &solDot,
                                   std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  // TODO : Verifier que les règles d'intégration soient identiques ou adapter
  // On prend les poids de geoSpace, ok pour le jacobien seulement ?
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> &w = geoSpace->getQuadratureWeights();
  // double             rho = _par[0];
  double mu = _par[1];
  int nFunctionsU = intSpace[_idU]->getNbFunctions();
  int nFunctionsV = intSpace[_idV]->getNbFunctions();
  int nFunctionsP = intSpace[_idP]->getNbFunctions();

  double J, u, v, p, dudx, dudy, dvdx, dvdy, Sxx, Sxy, Syx, Syy;
  std::vector<double> x(3);
  std::vector<double> f(2); // Forces volumiques

  for(int k = 0; k < nG; ++k) {
    // Forces volumiques
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
    if(_fct != nullptr) _fct->eval(tn, x, f);

    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, _dxdr);
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, _dxds);
    J = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    double drdx = _dxds[1] / J;
    double drdy = -_dxds[0] / J;
    double dsdx = -_dxdr[1] / J;
    double dsdy = _dxdr[0] / J;

    u = intSpace[_idU]->interpolateFieldAtQuadNode(sol[_idU], k);
    v = intSpace[_idV]->interpolateFieldAtQuadNode(sol[_idV], k);
    p = intSpace[_idP]->interpolateFieldAtQuadNode(sol[_idP], k);

    dudx = intSpace[_idU]->interpolateFieldAtQuadNode_rDerivative(sol[_idU], k) * drdx +
           intSpace[_idU]->interpolateFieldAtQuadNode_sDerivative(sol[_idU], k) * dsdx;
    dudy = intSpace[_idU]->interpolateFieldAtQuadNode_rDerivative(sol[_idU], k) * drdy +
           intSpace[_idU]->interpolateFieldAtQuadNode_sDerivative(sol[_idU], k) * dsdy;
    dvdx = intSpace[_idV]->interpolateFieldAtQuadNode_rDerivative(sol[_idV], k) * drdx +
           intSpace[_idV]->interpolateFieldAtQuadNode_sDerivative(sol[_idV], k) * dsdx;
    dvdy = intSpace[_idV]->interpolateFieldAtQuadNode_rDerivative(sol[_idV], k) * drdy +
           intSpace[_idV]->interpolateFieldAtQuadNode_sDerivative(sol[_idV], k) * dsdy;

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
                                   std::vector<double> &geoCoord, double c0, double tn, double **Ae,
                                   std::vector<std::vector<double>> &sol, std::vector<std::vector<double>> &solDot)
{
  // TODO : Verifier que les règles d'intégration soient identiques ou adapter
  // On prend les poids de geoSpace, ok pour le jacobien seulement ?
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> &w = geoSpace->getQuadratureWeights();
  // double             rho = _par[0];
  double mu = _par[1];
  int nFunctionsU = intSpace[_idU]->getNbFunctions();
  int nFunctionsV = intSpace[_idV]->getNbFunctions();
  int nFunctionsP = intSpace[_idP]->getNbFunctions();

  double jac, dudx, dudy, dvdx, dvdy;
  std::vector<double> x(3);
  std::vector<double> f(2, 0.); // Body forces

  for(int k = 0; k < nG; ++k) {
    // Body forces
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
    if(_fct != nullptr) _fct->eval(tn, x, f);

    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, _dxdr);
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, _dxds);
    jac = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    double drdx = _dxds[1] / jac;
    double drdy = -_dxds[0] / jac;
    double dsdx = -_dxdr[1] / jac;
    double dsdy = _dxdr[0] / jac;

    dudx = intSpace[_idU]->interpolateFieldAtQuadNode_rDerivative(sol[_idU], k) * drdx +
           intSpace[_idU]->interpolateFieldAtQuadNode_sDerivative(sol[_idU], k) * dsdx;
    dudy = intSpace[_idU]->interpolateFieldAtQuadNode_rDerivative(sol[_idU], k) * drdy +
           intSpace[_idU]->interpolateFieldAtQuadNode_sDerivative(sol[_idU], k) * dsdy;
    dvdx = intSpace[_idV]->interpolateFieldAtQuadNode_rDerivative(sol[_idV], k) * drdx +
           intSpace[_idV]->interpolateFieldAtQuadNode_sDerivative(sol[_idV], k) * dsdx;
    dvdy = intSpace[_idV]->interpolateFieldAtQuadNode_rDerivative(sol[_idV], k) * drdy +
           intSpace[_idV]->interpolateFieldAtQuadNode_sDerivative(sol[_idV], k) * dsdy;

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
    for(int i = 0; i < nFunctionsP; ++i) {
      _feP[i] = intSpace[_idP]->getFunctionAtQuadNode(i, k);
    }

    // Équations pour u
    for(int i = 0; i < nFunctionsU; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) {
        Ae[I][J++] += jac * w[k] * (_feUdx[i] * mu * 2. * _feUdx[j] + _feUdy[i] * mu * _feUdy[j]);
      }
      for(int j = 0; j < nFunctionsV; ++j) {
        Ae[I][J++] += jac * w[k] * (_feUdy[i] * mu * _feVdx[j]);
      }
      for(int j = 0; j < nFunctionsP; ++j) {
        Ae[I][J++] -= jac * w[k] * (_feUdx[i] * _feP[j]);
      }
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
      for(int j = 0; j < nFunctionsP; ++j) {
        Ae[I][J++] -= jac * w[k] * (_feVdy[i] * _feP[j]);
      }
      ++I;
    }

    // Équations pour p
    for(int i = 0; i < nFunctionsP; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) {
        Ae[I][J++] -= jac * w[k] * (_feP[i] * _feUdx[j]);
      }
      for(int j = 0; j < nFunctionsV; ++j) {
        Ae[I][J++] -= jac * w[k] * (_feP[i] * _feVdy[j]);
      }
      for(int j = 0; j < nFunctionsP; ++j) {
        Ae[I][J++] += 0.0;
      }
      ++I;
    }
  }
}

void feSysElm_2D_NavierStokes::createElementarySystem(std::vector<feSpace *> &space)
{
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
  _dxdr.resize(3); // [dx/dr, dy/dr, dz/dr]
  _dxds.resize(3); // [dx/ds, dy/ds, dz/ds]
}

void feSysElm_2D_NavierStokes::computeBe(std::vector<double> &Ja, int numElem,
                                         std::vector<feSpace *> &intSpace, feSpace *geoSpace,
                                         std::vector<double> &geoCoord, double c0, double tn,
                                         double dt, double *Be, std::vector<std::vector<double>> &sol,
                                         std::vector<std::vector<double>> &solDot,
                                         std::vector<std::vector<double>> &solPrev, std::vector<std::vector<double>> &solNext)
{
  // TODO : Verifier que les règles d'intégration soient identiques ou adapter
  // On prend les poids de geoSpace, ok pour le jacobien seulement ?
  int nG = geoSpace->getNbQuadPoints();
  std::vector<double> &w = geoSpace->getQuadratureWeights();
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

    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, _dxdr);
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, _dxds);
    J = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    double drdx = _dxds[1] / J;
    double drdy = -_dxds[0] / J;
    double dsdx = -_dxdr[1] / J;
    double dsdy = _dxdr[0] / J;

    u = intSpace[_idU]->interpolateFieldAtQuadNode(sol[_idU], k);
    v = intSpace[_idV]->interpolateFieldAtQuadNode(sol[_idV], k);
    p = intSpace[_idP]->interpolateFieldAtQuadNode(sol[_idP], k);

    dudt = intSpace[_idU]->interpolateFieldAtQuadNode(solDot[_idU], k);
    dvdt = intSpace[_idV]->interpolateFieldAtQuadNode(solDot[_idV], k);

    dudx = intSpace[_idU]->interpolateFieldAtQuadNode_rDerivative(sol[_idU], k) * drdx +
           intSpace[_idU]->interpolateFieldAtQuadNode_sDerivative(sol[_idU], k) * dsdx;
    dudy = intSpace[_idU]->interpolateFieldAtQuadNode_rDerivative(sol[_idU], k) * drdy +
           intSpace[_idU]->interpolateFieldAtQuadNode_sDerivative(sol[_idU], k) * dsdy;
    dvdx = intSpace[_idV]->interpolateFieldAtQuadNode_rDerivative(sol[_idV], k) * drdx +
           intSpace[_idV]->interpolateFieldAtQuadNode_sDerivative(sol[_idV], k) * dsdx;
    dvdy = intSpace[_idV]->interpolateFieldAtQuadNode_rDerivative(sol[_idV], k) * drdy +
           intSpace[_idV]->interpolateFieldAtQuadNode_sDerivative(sol[_idV], k) * dsdy;

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
                                         double **Ae, std::vector<std::vector<double>> &sol,
                                         std::vector<std::vector<double>> &solDot)
{
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
  std::vector<double> x(3);
  std::vector<double> f(2, 0.); // Forces volumiques

  for(int k = 0; k < nG; ++k) {
    // Forces volumiques
    geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
    if(_fct != nullptr) _fct->eval(tn, x, f);

    geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, _dxdr);
    geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, _dxds);
    jac = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    double drdx = _dxds[1] / jac;
    double drdy = -_dxds[0] / jac;
    double dsdx = -_dxdr[1] / jac;
    double dsdy = _dxdr[0] / jac;

    u = intSpace[_idU]->interpolateFieldAtQuadNode(sol[_idU], k);
    v = intSpace[_idV]->interpolateFieldAtQuadNode(sol[_idV], k);

    dudt = intSpace[_idU]->interpolateFieldAtQuadNode(solDot[_idU], k);
    dvdt = intSpace[_idV]->interpolateFieldAtQuadNode(solDot[_idV], k);

    dudx = intSpace[_idU]->interpolateFieldAtQuadNode_rDerivative(sol[_idU], k) * drdx +
           intSpace[_idU]->interpolateFieldAtQuadNode_sDerivative(sol[_idU], k) * dsdx;
    dudy = intSpace[_idU]->interpolateFieldAtQuadNode_rDerivative(sol[_idU], k) * drdy +
           intSpace[_idU]->interpolateFieldAtQuadNode_sDerivative(sol[_idU], k) * dsdy;
    dvdx = intSpace[_idV]->interpolateFieldAtQuadNode_rDerivative(sol[_idV], k) * drdx +
           intSpace[_idV]->interpolateFieldAtQuadNode_sDerivative(sol[_idV], k) * dsdx;
    dvdy = intSpace[_idV]->interpolateFieldAtQuadNode_rDerivative(sol[_idV], k) * drdy +
           intSpace[_idV]->interpolateFieldAtQuadNode_sDerivative(sol[_idV], k) * dsdy;

    int I = 0, J = 0;

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
    for(int i = 0; i < nFunctionsP; ++i) {
      _feP[i] = intSpace[_idP]->getFunctionAtQuadNode(i, k);
    }

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
      for(int j = 0; j < nFunctionsP; ++j) {
        Ae[I][J++] -= jac * w[k] * (_feUdx[i] * _feP[j]);
      }
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
      for(int j = 0; j < nFunctionsP; ++j) {
        Ae[I][J++] -= jac * w[k] * (_feVdy[i] * _feP[j]);
      }
      ++I;
    }

    // Équations pour p
    for(int i = 0; i < nFunctionsP; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) {
        Ae[I][J++] -= jac * w[k] * (_feP[i] * _feUdx[j]);
      }
      for(int j = 0; j < nFunctionsV; ++j) {
        Ae[I][J++] -= jac * w[k] * (_feP[i] * _feVdy[j]);
      }
      for(int j = 0; j < nFunctionsP; ++j) {
        Ae[I][J++] += 0.0;
      }
      ++I;
    }
  }
}