#include "feSysElm.h"
#include "feBilinearForm.h"

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

void feSysElm_1D_weakBC_edo1::computeAe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  int nFunctionsU = form->_intSpace[_idU]->getNbFunctions();
  int nFunctionsL = form->_intSpace[_idL]->getNbFunctions();
  double jac;

  for(int k = 0; k < nG; ++k) {
    std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, dxdr);
    jac = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);
    // jac = J[nG * numElem + k];
    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = form->_intSpace[_idU]->getFunctionAtQuadNode(i, k);
      _feV[i] = form->_intSpace[_idV]->getFunctionAtQuadNode(i, k);
    }

    for(int i = 0; i < nFunctionsL; ++i) {
      _feL[i] = form->_intSpace[_idL]->getFunctionAtQuadNode(i, k);
    }
    for(int i = 0; i < nFunctionsU; ++i) {
      for(int j = 0; j < nFunctionsL; ++j) {
        form->_Ae[i][nFunctionsU + j] += _feU[i] * _feL[j] * jac * w[k];
      }
    }
    for(int i = 0; i < nFunctionsL; ++i) {
      for(int j = 0; j < nFunctionsU; ++j) {
        form->_Ae[nFunctionsU + i][j] += _feU[j] * _feL[i] * jac * w[k];
      }
    }
    for(int i = 0; i < nFunctionsL; ++i) {
      for(int j = 0; j < nFunctionsU; ++j) {
        form->_Ae[nFunctionsU + i][(nFunctionsU + nFunctionsL) + j] += -_feL[i] * _feV[j] * jac * w[k];
      }
    }

    for(int i = 0; i < nFunctionsU; ++i) {
      for(int j = 0; j < nFunctionsU; ++j) {
        form->_Ae[(nFunctionsU + nFunctionsL) + i][(nFunctionsU + nFunctionsL) + j] +=
          form->_c0 * _feV[i] * _feV[j] * jac * w[k];
      }
    }
  }
}

void feSysElm_1D_weakBC_edo1::computeBe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  // Make the difference between different type of space
  int nFunctionsU = form->_intSpace[_idU]->getNbFunctions();
  int nFunctionsL = form->_intSpace[_idL]->getNbFunctions();

  std::vector<double> x(3, 0.0);
  double u, l, v, vDot, jac;
  // std::cout<<"gammaDot= "<<gammaDot<< " x= "<<x[0]<< "tn= " <<tn<<std::endl;
  // printf("%16.6e \n" , v);
  for(int k = 0; k < nG; ++k) {
    std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, dxdr);
    jac = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, x);
    double gammaDot = _fct->eval(form->_tn, x);
    u = form->_intSpace[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    v = form->_intSpace[_idV]->interpolateFieldAtQuadNode(form->_sol[_idV], k);
    l = form->_intSpace[_idL]->interpolateFieldAtQuadNode(form->_sol[_idL], k);
    vDot = form->_intSpace[_idV]->interpolateFieldAtQuadNode(form->_solDot[_idV], k);
    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = form->_intSpace[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= _feU[i] * l * jac * w[k];
    }
    for(int i = 0; i < nFunctionsL; ++i) {
      _feL[i] = form->_intSpace[_idL]->getFunctionAtQuadNode(i, k);
      form->_Be[nFunctionsU + i] -= _feL[i] * (u - v) * jac * w[k];
    }
    for(int i = 0; i < nFunctionsU; ++i) {
      _feV[i] = form->_intSpace[_idV]->getFunctionAtQuadNode(i, k);
      form->_Be[(nFunctionsU + nFunctionsL) + i] -= _feV[i] * (vDot - gammaDot) * jac * w[k];
    }
  }
}

void feSysElm_1D_Source::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _iVar[0] = _idU;
  _jVar[0] = _idU;
  _feU.resize(space[_idU]->getNbFunctions());
  _feUdx.resize(space[_idU]->getNbFunctions());
}

void feSysElm_1D_Source::computeBe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  int nFunctions = form->_intSpace[_idU]->getNbFunctions();
  std::vector<double> &J = form->_cnc->getJacobians();

  std::vector<double> x(3, 0.0);

  for(int k = 0; k < nG; ++k) {
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, x);
    double S = _fct->eval(form->_tn, x);
    for(int i = 0; i < nFunctions; ++i) {
      _feU[i] = form->_intSpace[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= _feU[i] * S * J[nG * form->_numElem + k] * w[k];
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

void feSysElm_1D_Diffusion::computeAe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  double kD = _par;
  int nFunctions = form->_intSpace[_idU]->getNbFunctions();
  std::vector<double> &J = form->_cnc->getJacobians();

  double jac;
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * form->_numElem + k];

    for(int i = 0; i < nFunctions; ++i)
      _feUdx[i] = form->_intSpace[_idU]->getdFunctiondrAtQuadNode(i, k) / jac;

    for(int i = 0; i < nFunctions; ++i) {
      for(int j = 0; j < nFunctions; ++j) {
        form->_Ae[i][j] += _feUdx[i] * kD * _feUdx[j] * jac * w[k];
      }
    }
  }
}

void feSysElm_1D_Diffusion::computeBe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  double kD = _par;
  int nFunctions = form->_intSpace[_idU]->getNbFunctions();
  std::vector<double> &J = form->_cnc->getJacobians();

  double jac, dudx;
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * form->_numElem + k];

    dudx = form->_intSpace[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k);
    dudx /= jac;

    for(int i = 0; i < nFunctions; ++i) {
      _feUdx[i] = form->_intSpace[_idU]->getdFunctiondrAtQuadNode(i, k);
      _feUdx[i] /= jac;
      form->_Be[i] -= _feUdx[i] * kD * dudx * jac * w[k];
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

void feSysElm_1D_Masse::computeAe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  double rho = _par;
  int nFunctions = form->_intSpace[_idU]->getNbFunctions();
  std::vector<double> &J = form->_cnc->getJacobians();

  double jac;
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * form->_numElem + k];

    for(int i = 0; i < nFunctions; ++i)
      _feU[i] = form->_intSpace[_idU]->getFunctionAtQuadNode(i, k);

    for(int i = 0; i < nFunctions; ++i)
      for(int j = 0; j < nFunctions; ++j)
        form->_Ae[i][j] += _feU[i] * rho * form->_c0 * _feU[j] * jac * w[k];
  }
}

void feSysElm_1D_Masse::computeBe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  double rho = _par;
  int nFunctions = form->_intSpace[_idU]->getNbFunctions();
  std::vector<double> &J = form->_cnc->getJacobians();

  double jac, uDot;
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * form->_numElem + k];

    uDot = form->_intSpace[_idU]->interpolateFieldAtQuadNode(form->_solDot[_idU], k);

    for(int i = 0; i < nFunctions; ++i) {
      _feU[i] = form->_intSpace[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= _feU[i] * rho * uDot * jac * w[k];
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

void feSysElm_1D_NeumannBC::computeBe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  int nFunctions = form->_intSpace[_idU]->getNbFunctions();
  std::vector<double> &J = form->_cnc->getJacobians();

  double jac;
  std::vector<double> x(3, 0.0);
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * form->_numElem + k];
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, x);
    double h = _fct->eval(form->_tn, x);
    for(int i = 0; i < nFunctions; ++i) {
      _feU[i] = form->_intSpace[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] += _feU[i] * h * jac * w[k];
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

void feSysElm_1D_Advection::computeBe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  int nFunctions = form->_intSpace[_idU]->getNbFunctions();
  std::vector<double> &J = form->_cnc->getJacobians();

  double u, dudx, jac, cVelocity;
  std::vector<double> x(3, 0.0);
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * form->_numElem + k];
    u = form->_intSpace[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    dudx = form->_intSpace[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) / jac;

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, x);
    cVelocity = _fct->eval(form->_tn, x);
    for(int i = 0; i < nFunctions; ++i) {
      _feUdx[i] = form->_intSpace[_idU]->getdFunctiondrAtQuadNode(i, k);
      _feUdx[i] /= jac;
      form->_Be[i] += cVelocity * u * _feUdx[i] * jac * w[k];
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

void feSysElm_1D_SUPGStab::computeAe(feBilinearForm *form)
{
  // Finite differences for now
}

void feSysElm_1D_SUPGStab::computeBe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  int nFunctions = form->_intSpace[_idU]->getNbFunctions();
  std::vector<double> &J = form->_cnc->getJacobians();

  double kDiffusivity = _par[0];

  double u, dudx, jac, cVelocity;
  std::vector<double> x(3, 0.0);
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * form->_numElem + k];
    u    = form->_intSpace[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    dudx = form->_intSpace[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) / jac;

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, x);
    cVelocity = _fct->eval(form->_tn, x);
    for(int i = 0; i < nFunctions; ++i) {
      _feUdx[i] = form->_intSpace[_idU]->getdFunctiondrAtQuadNode(i, k);
      _feUdx[i] /= jac;

      // Compute element Peclet and optimal numerical diffusivity coth(Peh) - 1/Peh
      double he = 2.0 * jac;
      double Peh = cVelocity * he / (2. * kDiffusivity);
      double beta = cosh(Peh)/sinh(Peh) - (1./Peh);
      double tau = beta * he / (2. * cVelocity);
      double residual = cVelocity * dudx - 1.; // FIXME : completer avec d2dx2
      form->_Be[i] -= cVelocity * _feUdx[i] * tau * residual * jac * w[k];
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

void feSysElm_1D_DG_Advection::computeAe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  double jac;
  int nFunctions = form->_intSpace[_idU]->getNbFunctions();
  std::vector<double> &J = form->_cnc->getJacobians();
  std::vector<double> x(3, 0.0);

  // Element term: integral of c * dphi_i/dx * phi_j
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * form->_numElem + k];

    // Initialize vectors phi (_feU) and dphi/dx = dphi/dr * dr/dx (1/J) (_feUdx)
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

void feSysElm_1D_DG_Advection::computeBe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNbQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  int nFunctions = form->_intSpace[_idU]->getNbFunctions();
  std::vector<double> &J = form->_cnc->getJacobians();

  double u, jac, cVelocity;
  std::vector<double> x(3, 0.0);
  for(int k = 0; k < nG; ++k) {
    jac = J[nG * form->_numElem + k];
    u = form->_intSpace[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, x);
    cVelocity = _fct->eval(form->_tn, x);
    for(int i = 0; i < nFunctions; ++i) {
      _feUdx[i] = form->_intSpace[_idU]->getdFunctiondrAtQuadNode(i, k);
      _feUdx[i] /= jac;
      form->_Be[i] += cVelocity * u * _feUdx[i] * jac * w[k];
      // feInfo("residual = %f", cVelocity * u * _feUdx[i] * jac * w[k]);
    }
  }

  double rLeft[3]  = {-1.0, 0., 0.};
  double rRight[3] = {+1.0, 0., 0.};
  std::vector<double> phiL(nFunctions);
  std::vector<double> phiR(nFunctions);
  form->_intSpace[_idU]->L(rLeft,  phiL.data());
  form->_intSpace[_idU]->L(rRight, phiR.data());

  double uMinus, uPlus;
  double uL, uR, uPrevR, uNextL;
  uL = form->_intSpace[_idU]->interpolateField(form->_sol[_idU], rLeft);
  uR = form->_intSpace[_idU]->interpolateField(form->_sol[_idU], rRight);
  uPrevR = form->_intSpace[_idU]->interpolateField(form->_solPrev[_idU], rRight);
  uNextL = form->_intSpace[_idU]->interpolateField(form->_solNext[_idU], rLeft);

  // Moyenne:
  uPlus  = (uR + uNextL)/2.;
  uMinus = (uL + uPrevR)/2.;

  // Interface (flux) term: - [(c*uh)*phi_i]_X^e-1 ^X^e
  for(int i = 0; i < nFunctions; ++i) {
    // if(cVelocity > 0){
    //   // Information comes from the left:
    //   // Interpolate solution on the right of the previous element
    //   uMinus = form->_intSpace[_idU]->interpolateField(form->_solPrev[_idU], rRight);
    //   uPlus  = form->_intSpace[_idU]->interpolateField(    sol[_idU], rRight);
    // } else{
    //   // Information comes from the right:
    //   // Interpolate solution on the left of the next element
    //   uMinus = form->_intSpace[_idU]->interpolateField(    sol[_idU], rLeft);
    //   uPlus  = form->_intSpace[_idU]->interpolateField(form->_solNext[_idU], rLeft);
    // }

    form->_Be[i] -= (cVelocity * uPlus) * phiR[i] - (cVelocity * uMinus) * phiL[i];

    // feInfo("Elem %d - u+ = %f - u- = %f",
    //   uPlus,
    //   uMinus);

    // feInfo("Elem %d - Adding flux %d = %f - %f",
    //   form->_numElem,
    //   i,
    //   (cVelocity * uPlus) * phiR[i] - (cVelocity * uMinus) * phiL[i],
    //   uPlus - uMinus);
  }
}