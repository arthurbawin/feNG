#include "feSysElm.h"
#include "feSysElm_Ext.h"
#include "feBilinearForm.h"

// -----------------------------------------------------------------------------
// 1D Neumann boundary condition (= source)
// -----------------------------------------------------------------------------
void feSysElm_1D_NeumannBC::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _feU.resize(_nFunctions);
}

void feSysElm_1D_NeumannBC::computeBe(feBilinearForm *form)
{
  for(int k = 0; k < _nQuad; ++k) {
    double jac = form->_J[_nQuad * form->_numElem + k];
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    double h = _neumannBC->eval(form->_tn, _pos);
    for(int i = 0; i < _nFunctions; ++i) {
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= _feU[i] * h * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// 1D Advection
// -----------------------------------------------------------------------------
void feSysElm_1D_Advection::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _feU.resize(_nFunctions);
  _feUdx.resize(_nFunctions);
}

void feSysElm_1D_Advection::computeAe(feBilinearForm *form)
{
  double jac, cVelocity;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    for(int i = 0; i < _nFunctions; ++i){
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) / jac;
    }

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    cVelocity = _cVelocity->eval(form->_tn, _pos);

    for(int i = 0; i < _nFunctions; ++i) {
      for(int j = 0; j < _nFunctions; ++j) {
        form->_Ae[i][j] -= cVelocity * _feUdx[i] * _feU[j] * jac * _wQuad[k];
      }
    }
  }
}

void feSysElm_1D_Advection::computeBe(feBilinearForm *form)
{
  double u, jac, cVelocity;
  std::vector<double> x(3, 0.0);
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    cVelocity = _cVelocity->eval(form->_tn, _pos);
    for(int i = 0; i < _nFunctions; ++i) {
      _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k);
      _feUdx[i] /= jac;
      form->_Be[i] += cVelocity * u * _feUdx[i] * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// 1D SUPG Stabilization
// -----------------------------------------------------------------------------
void feSysElm_1D_SUPGStab::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _feU.resize(_nFunctions);
  _feUdx.resize(_nFunctions);
}

void feSysElm_1D_SUPGStab::computeAe(feBilinearForm *form)
{
  // Finite differences for now
}

void feSysElm_1D_SUPGStab::computeBe(feBilinearForm *form)
{
  double u, dudx, jac, cVelocity;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    u    = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    dudx = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k);
    dudx /= jac;

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    cVelocity = (*_velocity)(form->_tn, _pos);

    // Compute element Peclet and optimal numerical diffusivity coth(Peh) - 1/Peh
    double he = 2.0 * jac;
    double Peh = cVelocity * he / (2. * _diffusivity);
    double beta = cosh(Peh)/sinh(Peh) - (1./Peh);
    double tau = beta * he / (2. * cVelocity);
    // FIXME : completer avec d2dx2
    // FIXME : choisir le signe du terme source
    double residual = cVelocity * dudx - (*_source)(form->_tn, _pos); 

    for(int i = 0; i < _nFunctions; ++i) {

      _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k);
      _feUdx[i] /= jac;

      form->_Be[i] -= residual * tau * cVelocity * _feUdx[i] * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// 1D DG Advection
// -----------------------------------------------------------------------------
void feSysElm_1D_DG_Advection::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
}

void feSysElm_1D_DG_Advection::computeAe(feBilinearForm *form)
{
  // double rLeft[3]  = {-1.0, 0., 0.};
  // double rRight[3] = {+1.0, 0., 0.};

  // std::vector<double> vLeft(3, 0.0);
  // std::vector<double> vRight(3, 0.0);
  // // Evaluate the imposed velocity field
  // form->_geoSpace->interpolateVectorField(form->_geoCoord, rLeft, _pos);
  // (*_velocity)(form->_tn, _pos, vLeft);
  // form->_geoSpace->interpolateVectorField(form->_geoCoord, rRight, _pos);
  // (*_velocity)(form->_tn, _pos, vRight);

  form->_Ae[0][0] = 1.;
  form->_Ae[0][1] = 1.;
  form->_Ae[1][0] = 1.;
  form->_Ae[1][1] = 1.;
}

static double LFflux(double uM, double uP, feFunction *f, feFunction *df, std::vector<double> &pos)
{
  double alpha = fmax(fabs((*df)(uP, pos)), fabs((*df)(uM, pos)));
  return 0.5 * ((*f)(uM, pos) + (*f)(uP, pos) - alpha * (uP - uM));
}

void feSysElm_1D_DG_Advection::computeBe(feBilinearForm *form)
{
  // DG correction
  double rLeft[3]  = {-1.0, 0., 0.};
  double rRight[3] = {+1.0, 0., 0.};

  std::vector<double> vLeft(3, 0.0);
  std::vector<double> vRight(3, 0.0);
  // Evaluate the imposed velocity field
  form->_geoSpace->interpolateVectorField(form->_geoCoord, rLeft, _pos);
  (*_velocity)(form->_tn, _pos, vLeft);
  form->_geoSpace->interpolateVectorField(form->_geoCoord, rRight, _pos);
  (*_velocity)(form->_tn, _pos, vRight);

  std::vector<double> phiPlus(_nFunctions);  // Phi+(x j-1/2) (left)
  std::vector<double> phiMinus(_nFunctions); // Phi-(x j+1/2) (right)
  form->_intSpaces[_idU]->L(rLeft,  phiPlus.data());
  form->_intSpaces[_idU]->L(rRight, phiMinus.data());

  double uL, uR, uPrevR, uNextL;
  uL     = form->_intSpaces[_idU]->interpolateField(form->_sol[_idU]    , rLeft);  // u+(j-1/2)
  uR     = form->_intSpaces[_idU]->interpolateField(form->_sol[_idU]    , rRight); // u-(j+1/2)
  uPrevR = form->_intSpaces[_idU]->interpolateField(form->_solPrev[_idU], rRight); // u-(j-1/2)
  uNextL = form->_intSpaces[_idU]->interpolateField(form->_solNext[_idU], rLeft);  // u+(j+1/2)

  // Flux is upwinding cu- (assuming f(u) = cu)
  // FIXME: add generic flux (Lax-Friedrichs or other)
  double fMinus = LFflux(uR, uNextL, _pdeFlux, _pdedFlux, _pos);
  double  fPlus = LFflux(uPrevR, uL, _pdeFlux, _pdedFlux, _pos);
  // fMinus = vRight[0] * uR;
  // fPlus = vLeft[0] * uPrevR;
  // fMinus = uR*uR;
  // fPlus = uPrevR * uPrevR;

  for(int i = 0; i < _nFunctions; ++i) {
    form->_Be[i] -= fMinus * phiMinus[i] - fPlus * phiPlus[i];
  }
}

// -----------------------------------------------------------------------------
// 1D Euler-Bernoulli
// -----------------------------------------------------------------------------
void feSysElm_1D_Beam::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _feU.resize(_nFunctions);
  _feUdx2.resize(_nFunctions);
}

void feSysElm_1D_Beam::computeAe(feBilinearForm *form)
{
  /* ... */
}

void feSysElm_1D_Beam::computeBe(feBilinearForm *form)
{
  double jac, d2udx2;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    d2udx2 = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rrDerivative(form->_sol[_idU], k);
    d2udx2 /= jac*jac;

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    double EI = _par;

    for(int i = 0; i < _nFunctions; ++i) {
      _feUdx2[i] = form->_intSpaces[_idU]->getd2Functiondr2AtQuadNode(i, k);
      _feUdx2[i] /= jac*jac;
      form->_Be[i] -= _feUdx2[i] * EI * d2udx2 * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Undocumented
// -----------------------------------------------------------------------------
void feSysElm_1D_weakBC_edo1::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idL = 1;
  _idV = 2;
  _fieldsLayoutI.resize(space.size());
  _fieldsLayoutJ.resize(space.size());
  _feU.resize(space[_idU]->getNumFunctions());
  _feL.resize(space[_idL]->getNumFunctions());
  _feV.resize(space[_idV]->getNumFunctions());
  for(int i = 0; i < space.size(); i++) {
    _fieldsLayoutI[i] = i;
    _fieldsLayoutJ[i] = i;
  }
}

void feSysElm_1D_weakBC_edo1::computeAe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNumQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  int nFunctionsU = form->_intSpaces[_idU]->getNumFunctions();
  int nFunctionsL = form->_intSpaces[_idL]->getNumFunctions();
  double jac;

  for(int k = 0; k < nG; ++k) {
    std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, dxdr);
    jac = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);
    // jac = J[nG * numElem + k];
    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      _feV[i] = form->_intSpaces[_idV]->getFunctionAtQuadNode(i, k);
    }

    for(int i = 0; i < nFunctionsL; ++i) {
      _feL[i] = form->_intSpaces[_idL]->getFunctionAtQuadNode(i, k);
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
  int nG = form->_geoSpace->getNumQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  // Make the difference between different type of space
  int nFunctionsU = form->_intSpaces[_idU]->getNumFunctions();
  int nFunctionsL = form->_intSpaces[_idL]->getNumFunctions();

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
    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    v = form->_intSpaces[_idV]->interpolateFieldAtQuadNode(form->_sol[_idV], k);
    l = form->_intSpaces[_idL]->interpolateFieldAtQuadNode(form->_sol[_idL], k);
    vDot = form->_intSpaces[_idV]->interpolateFieldAtQuadNode(form->_solDot[_idV], k);
    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= _feU[i] * l * jac * w[k];
    }
    for(int i = 0; i < nFunctionsL; ++i) {
      _feL[i] = form->_intSpaces[_idL]->getFunctionAtQuadNode(i, k);
      form->_Be[nFunctionsU + i] -= _feL[i] * (u - v) * jac * w[k];
    }
    for(int i = 0; i < nFunctionsU; ++i) {
      _feV[i] = form->_intSpaces[_idV]->getFunctionAtQuadNode(i, k);
      form->_Be[(nFunctionsU + nFunctionsL) + i] -= _feV[i] * (vDot - gammaDot) * jac * w[k];
    }
  }
}

void feSysElm_1D_Euler::createElementarySystem(std::vector<feSpace *> &space)
{
  _idr  = 0;
  _idru = 1;
  _idre = 2;
  _fieldsLayoutI = {_idr, _idru, _idre};
  _fieldsLayoutJ = {_idr, _idru, _idre};
  _nFunctionsr = space[_idr]->getNumFunctions();
  _nFunctionsru = space[_idru]->getNumFunctions();
  _nFunctionsre = space[_idre]->getNumFunctions();
  _phir.resize(_nFunctionsr);
  _gradphir.resize(_nFunctionsr);
  _phiru.resize(_nFunctionsru);
  _gradphiru.resize(_nFunctionsru);
  _phire.resize(_nFunctionsre);
  _gradphire.resize(_nFunctionsre);
}

void feSysElm_1D_Euler::computeAe(feBilinearForm *form)
{
  // Finite differences
}

static double area(double x)
{
  return 1. + 2.2 * (x - 1.5) * (x - 1.5); // Anderson
  // return 1.398 + 0.347 * tanh(0.8*x-4.);
}

static double dadx(double x)
{
  return 2. * 2.2 * (x - 1.5);
  // return 0.2776/(cosh(4. - 0.8 * x)*cosh(4. - 0.8 * x));
}

void feSysElm_1D_Euler::computeBe(feBilinearForm *form)
{
  double gamma = 1.4;
  double jac, rhoA, rhouA, eA, rho, u, e, p;

  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    double A = area(_pos[0]);
    double dAdx = dadx(_pos[0]);
    
    // Compute current fields and gradient of test functions
    rhoA = form->_intSpaces[_idr ]->interpolateFieldAtQuadNode(form->_sol[_idr ], k);
    rhouA = form->_intSpaces[_idru]->interpolateFieldAtQuadNode(form->_sol[_idru], k);
    eA = form->_intSpaces[_idre]->interpolateFieldAtQuadNode(form->_sol[_idre], k);

    rho = rhoA/A;
    u = rhouA/(rhoA);
    e = eA/A;
    p = (e - rho*u*u/2.) * (gamma - 1.);

    form->_intSpaces[_idru]->getFunctionsAtQuadNode(k, _phiru);

    form->_intSpaces[_idr]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, 
      form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, _gradphir.data());
    form->_intSpaces[_idru]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, 
      form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, _gradphiru.data());
    form->_intSpaces[_idre]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, 
      form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, _gradphire.data());

    double flux[3];
    flux[0] = rhouA;
    flux[1] = (rho*u*u+p)*A;
    flux[2] = u*A*(e+p);

    int cnt = 0;
    // Residu pour rhoA
    for(int i = 0; i < _nFunctionsr; ++i) {
      form->_Be[cnt++] += (flux[0] * _gradphir[i]) * jac * _wQuad[k];
    }
    // Residu pour rho*u
    for(int i = 0; i < _nFunctionsru; ++i) {
      form->_Be[cnt++] += (flux[1] * _gradphiru[i]) * jac * _wQuad[k];
      form->_Be[cnt++] += p * dAdx * _phiru[i] * jac * _wQuad[k]; // Source
    }
    // Residu pour rho*e
    for(int i = 0; i < _nFunctionsre; ++i) {
      form->_Be[cnt++] += (flux[2] * _gradphire[i]) * jac * _wQuad[k];
    }
  }
}

void feSysElm_0D_EulerFlux::createElementarySystem(std::vector<feSpace *> &space)
{
  _idrhoA_b  = 0;
  _idrhouA_b = 1;
  _ideA_b = 2;

  // Fields on the domain, to extrapolate boundary conditions
  _idrhoA_dom = 3;
  _idrhouA_dom = 4;
  _ideA_dom = 5;

  // Auxiliary FE spaces dont count for the local matrix/residual size
  _fieldsLayoutI = {_idrhoA_b, _idrhouA_b, _ideA_b};
  _fieldsLayoutJ = {_idrhoA_b, _idrhouA_b, _ideA_b};

  _nFunctionsr = space[_idrhoA_b]->getNumFunctions();
  _nFunctionsru = space[_idrhouA_b]->getNumFunctions();
  _nFunctionsre = space[_ideA_b]->getNumFunctions();

  _phir.resize(_nFunctionsr);
  _gradphir.resize(_nFunctionsr);
  _phiru.resize(_nFunctionsru);
  _gradphiru.resize(_nFunctionsru);
  _phire.resize(_nFunctionsre);
  _gradphire.resize(_nFunctionsre);
}

void feSysElm_0D_EulerFlux::computeAe(feBilinearForm *form)
{
  // Finite differences
}

void feSysElm_0D_EulerFlux::computeBe(feBilinearForm *form)
{
  double gamma = 1.4;
  double jac, rhoA, rhouA, eA;
  double cs;
  double rhob, ub, pb, eb;
  double flux[3];

  form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, 0, _pos);
  double A = area(_pos[0]);
  
  if(_boundaryType == 1)
  {
    // Inlet

    //  speed of sound from Riemann invariant
    //  temperature from isentropic relation
    //  pressure from isentropic relation  
    //  density from gas equation          
    //  velocity from energy equation

    // Extrapolate values from solNext: 2nd element on interior domain
    // Initialize solution at vertices
    double rLeft[3]  = {-1.0, 0., 0.};
    double rRight[3] = {+1.0, 0., 0.};

    double rhoA_l = form->_intSpaces[_idrhoA_dom]->interpolateField(form->_solNext[_idrhoA_dom], rLeft);
    double rhouA_l = form->_intSpaces[_idrhouA_dom]->interpolateField(form->_solNext[_idrhouA_dom], rLeft);
    double eA_l = form->_intSpaces[_ideA_dom]->interpolateField(form->_solNext[_ideA_dom], rLeft);

    double rhoA_r = form->_intSpaces[_idrhoA_dom]->interpolateField(form->_solNext[_idrhoA_dom], rRight);
    double rhouA_r = form->_intSpaces[_idrhouA_dom]->interpolateField(form->_solNext[_idrhouA_dom], rRight);
    double eA_r = form->_intSpaces[_ideA_dom]->interpolateField(form->_solNext[_ideA_dom], rRight);

    // Extrapolate values assuming constant spatial step h:
    // Get first jacobian on the domain
    double jac_domain = form->_intSpaces[_idrhoA_dom]->getCncGeo()->getJacobians()[0];
    double h = 2. * jac_domain;

    double Anext = area(_pos[0] + h);
    
    // double rhoA_ext = 2.*rhoA_l - rhoA_r;
    // double rhouA_ext = 2.*rhouA_l - rhouA_r;
    // double eA_ext = 2.*eA_l - eA_r;

    double rhoA_ext = rhoA_r;
    double rhouA_ext = rhouA_r;
    double eA_ext = eA_r;

    feInfo("%f - %f - %f", rhoA_ext, rhouA_ext, eA_ext);
    feInfo("%f - %f - %f", rhoA_ext/Anext, rhouA_ext/Anext, eA_ext/Anext);
    feInfo("%f - %f - %f", rhoA_r, rhouA_r, eA_r);

    // FIXME: peut-être juste prendre les valeurs _l de l'interieur
    double rho_e = rhoA_ext/Anext;
    double u_e = rhouA_ext / rhoA_ext;
    double e_e = eA_ext/Anext;
    double p_e = (e_e - rho_e*u_e*u_e/2.) * (gamma - 1.);

    feInfo("inlet : rho_e = %f - p_e = %f - u_e = %f - e_e = %f", rho_e, p_e, u_e, e_e);

    double cpgas = 1005.;
    double rgas = (gamma-1.) * cpgas/gamma;

    double t01 = 288.;
    double p01 = 0.72e5;

    double cs2  = gamma*p_e/rho_e;
    double c02  = cs2 + 0.5*(gamma - 1.)*u_e*u_e;
    double rinv = u_e - 2.0*sqrt(cs2)/(gamma - 1.);
    double dis  = (gamma + 1.)*c02/((gamma - 1.)*rinv*rinv) - 0.5*(gamma - 1.);

    feInfo("inlet : rinv = %f - cs2 = %f - c02 = %f - dis = %f", rinv, cs2, c02, dis);
    
    if(dis<0) {
      dis=1e-15;       
    }

    double cb   = -rinv*((gamma - 1.)/(gamma + 1.))*(1.0+sqrt(dis));
    double cc02 = fmin(cb*cb/c02, 1.0);
    double tb   = cc02*t01;
    pb   = p01*pow(tb/t01, gamma/(gamma-1.));
    rhob = pb/(rgas*tb);
    ub   = sqrt(2.0*cpgas*(t01-tb));

    feInfo("inlet  velocity = %f - cs = %f - mach = %f - rho = %f", ub, cb, ub/cb, rhob);

    flux[0] = rhob*ub*Anext;
    flux[1] = (rhob*ub*ub+pb)*Anext;
    flux[2] = ub*Anext*(pb/(gamma-1.) + rhob*ub*ub/2. + pb);
    feInfo("inlet  flux = %f %f %f (Anext = %f)", flux[0], flux[1], flux[2], Anext);

    int cnt = 0;
    // Flux pour rhoA
    form->_Be[cnt++] -= flux[0];
    // Flux pour rho*u
    form->_Be[cnt++] -= flux[1];
    // Flux pour rho*e
    form->_Be[cnt++] -= flux[2];
  } else
  {
    // Outlet

    // Extrapolate values from solNext: 2nd element on interior domain
    // Initialize solution at vertices
    double rLeft[3]  = {-1.0, 0., 0.};
    double rRight[3] = {+1.0, 0., 0.};

    double rhoA_l = form->_intSpaces[_idrhoA_dom]->interpolateField(form->_solPrev[_idrhoA_dom], rLeft);
    double rhouA_l = form->_intSpaces[_idrhouA_dom]->interpolateField(form->_solPrev[_idrhouA_dom], rLeft);
    double eA_l = form->_intSpaces[_ideA_dom]->interpolateField(form->_solPrev[_ideA_dom], rLeft);

    double rhoA_r = form->_intSpaces[_idrhoA_dom]->interpolateField(form->_solPrev[_idrhoA_dom], rRight);
    double rhouA_r = form->_intSpaces[_idrhouA_dom]->interpolateField(form->_solPrev[_idrhouA_dom], rRight);
    double eA_r = form->_intSpaces[_ideA_dom]->interpolateField(form->_solPrev[_ideA_dom], rRight);

    // Extrapolate values assuming constant spatial step h:
    // Get first jacobian on the domain
    double jac_domain = form->_intSpaces[_idrhoA_dom]->getCncGeo()->getJacobians()[0];
    double h = 2. * jac_domain;

    double Aprev = area(_pos[0] - h);
    
    // double rhoA_ext = 2.*rhoA_r - rhoA_l;
    // double rhouA_ext = 2.*rhouA_r - rhouA_l;
    // double eA_ext = 2.*eA_r - eA_l;

    double rhoA_ext = rhoA_r;
    double rhouA_ext = rhouA_r;
    double eA_ext = eA_r;

    feInfo("outlet %f - %f - %f", rhoA_ext, rhouA_ext, eA_ext);
    feInfo("outlet %f - %f - %f", rhoA_ext/Aprev, rhouA_ext/Aprev, eA_ext/Aprev);
    feInfo("outlet %f - %f - %f", rhoA_r, rhouA_r, eA_r);

    // FIXME: peut-être juste prendre les valeurs _l de l'interieur
    double rho_e = rhoA_ext/Aprev; // FIXME: prendre A_ext?
    double u_e = rhouA_ext / rhoA_ext;
    double e_e = eA_ext/Aprev;
    double p_e = (e_e - rho_e*u_e*u_e/2.) * (gamma - 1.);

    feInfo("outlet : rho_e = %f - p_e = %f - u_e = %f - e_e = %f", rho_e, p_e, u_e, e_e);

    cs  = sqrt(gamma * p_e/rho_e);

    feInfo("outlet velocity = %f - cs = %f - mach = %f - rho = %f", u_e, cs, u_e/cs, rho_e);
    feInfo("outlet e_e = %f - p_e = %f", e_e, p_e);

    if(u_e > cs){
      // Supersonic
      // Les variables _e remplacent les valeurs interpolees sans indices
      pb   = p_e;
      rhob = rho_e;
      ub   = u_e;
      feInfo("supersonic outlet: ub = %f - pb = %f - rhob = %f", ub, pb, rhob); 
    }
    else{
      // Subsonic
      double p_b = 0.7e5;
      pb   = p_b;
      rhob = rho_e + (p_b-p_e)/(cs*cs);
      ub   = u_e   - (p_b-p_e)/(cs*rho_e);
      feInfo("subsonic outlet: ub = %f - pb = %f - rhob = %f - p_e = %f", ub, pb, rhob, p_e);  
    }

    flux[0] = rhob*ub*Aprev;
    flux[1] = (rhob*ub*ub+pb)*Aprev;
    flux[2] = ub*Aprev*(pb/(gamma-1.) + rhob*ub*ub/2. + pb);
    feInfo("outlet flux = %f %f %f", flux[0], flux[1], flux[2]);

    int cnt = 0;
    // Flux pour rhoA
    form->_Be[cnt++] -= flux[0];
    // Flux pour rho*u
    form->_Be[cnt++] -= flux[1];
    // Flux pour rho*e
    form->_Be[cnt++] -= flux[2];
  }

  // int cnt = 0;
  // // Flux pour rhoA
  // form->_Be[cnt++] += flux[0];
  // // Flux pour rho*u
  // form->_Be[cnt++] += flux[1];
  // // Flux pour rho*e
  // form->_Be[cnt++] += flux[2];
}

void feSysElm_1D_Coupled::createElementarySystem(std::vector<feSpace *> &space)
{
  _idu = 0;
  _idv = 1;
  _fieldsLayoutI = {_idu, _idv};
  _fieldsLayoutJ = {_idu, _idv};
  _nFunctionsu = space[_idu]->getNumFunctions();
  _nFunctionsv = space[_idv]->getNumFunctions();
  _phiu.resize(_nFunctionsu);
  _gradphiu.resize(_nFunctionsu);
  _phiv.resize(_nFunctionsv);
  _gradphiv.resize(_nFunctionsv);
}

void feSysElm_1D_Coupled::computeAe(feBilinearForm *form)
{
  double jac;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    
    form->_intSpaces[_idu]->getFunctionsAtQuadNode(k, _phiu);
    form->_intSpaces[_idv]->getFunctionsAtQuadNode(k, _phiv);

    form->_intSpaces[_idu]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, 
      form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, _gradphiu.data());
    form->_intSpaces[_idv]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, 
      form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, _gradphiv.data());

    int I = 0, J = 0;
    // Residu pour u
    for(int i = 0; i < _nFunctionsu; ++i) {
      J = 0;
      for(int j = 0; j < _nFunctionsu; ++j) {
        form->_Ae[I][J++] -= (_phiu[j] * _gradphiu[i]) * jac * _wQuad[k];
      }
      for(int j = 0; j < _nFunctionsv; ++j) {
        form->_Ae[I][J++] -= (_phiv[j] * _phiu[i]) * jac * _wQuad[k];
      }
      ++I;
    }
    // Residu pour v
    for(int i = 0; i < _nFunctionsv; ++i) {
      J = 0;
      for(int j = 0; j < _nFunctionsu; ++j) {
        form->_Ae[I][J++] -= (_phiu[j] * _phiv[i]) * jac * _wQuad[k];
      }
      for(int j = 0; j < _nFunctionsv; ++j) {
        form->_Ae[I][J++] -= (_phiv[j] * _gradphiv[i]) * jac * _wQuad[k];
      }
      ++I;
    }
  }
}

void feSysElm_1D_Coupled::computeBe(feBilinearForm *form)
{
  double jac;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    
    // Compute current fields and gradient of test functions
    double u = form->_intSpaces[_idu]->interpolateFieldAtQuadNode(form->_sol[_idu], k);
    double v = form->_intSpaces[_idv]->interpolateFieldAtQuadNode(form->_sol[_idv], k);

    form->_intSpaces[_idu]->getFunctionsAtQuadNode(k, _phiu);
    form->_intSpaces[_idv]->getFunctionsAtQuadNode(k, _phiv);

    form->_intSpaces[_idu]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, 
      form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, _gradphiu.data());
    form->_intSpaces[_idv]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, 
      form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, _gradphiv.data());

    int cnt = 0;
    // Residu pour u
    for(int i = 0; i < _nFunctionsu; ++i) {
      form->_Be[cnt++] += (u * _gradphiu[i] + v * _phiu[i]) * jac * _wQuad[k];
    }
    // Residu pour v
    for(int i = 0; i < _nFunctionsv; ++i) {
      form->_Be[cnt++] += (v * _gradphiv[i] + u * _phiv[i]) * jac * _wQuad[k];
    }
  }
}