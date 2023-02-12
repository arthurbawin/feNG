#include "feSysElm.h"
#include "feBilinearForm.h"


void feSysElm_2D_Advection::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _feU.resize(space[_idU]->getNumFunctions());
  _feUdx.resize(space[_idU]->getNumFunctions());
  _feUdy.resize(space[_idU]->getNumFunctions());
  _dxdr.resize(3);
  _dxds.resize(3);
}

void feSysElm_2D_Advection::computeAe(feBilinearForm *form)
{
  // Calculee par differences finies
}

void feSysElm_2D_Advection::computeBe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNumQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  int nFunctions = form->_intSpaces[_idU]->getNumFunctions();

  double Jac, u, dudt, dudx, dudy;
  for(int k = 0; k < nG; ++k) {

    // Evaluate the exterior velocity field
    std::vector<double> v(2, 0.0);
    std::vector<double> x(3, 0.0);
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, x);
    (*_velocity)(form->_tn, x, v);

    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, _dxdr);
    form->_geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(form->_geoCoord, k, _dxds);
    Jac = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    double drdx = _dxds[1] / Jac;
    double drdy = -_dxds[0] / Jac;
    double dsdx = -_dxdr[1] / Jac;
    double dsdy = _dxdr[0] / Jac;

    // Compute SUPG parameter
    // double tau1 = 0.0, tau2 = 0.0, tau3 = 0.0, c = 0.0, cT = 0.0, kT = 0.0;
    // for(int i = 0; i < nFunctions; ++i) {
    //   dudt = form->_intSpaces[_idU]->interpolateSolutionDotAtQuadNode(k);
    //   dudx = form->_intSpaces[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdx +
    //          form->_intSpaces[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdx;
    //   dudy = form->_intSpaces[_idU]->interpolateSolutionAtQuadNode_rDerivative(k) * drdy +
    //          form->_intSpaces[_idU]->interpolateSolutionAtQuadNode_sDerivative(k) * dsdy;

    //   _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
    //   _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
    //               form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
    //   _feUdy[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
    //               form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
    //   c  += J * w[k] * _feU[i] * (v[0] * dudx + v[1] * dudy);
    //   cT += J * w[k] * (v[0] * _feUdx[i] + v[1] * _feUdy[i]) * dudt;
    //   kT += J * w[k] * (v[0] * _feUdx[i] + v[1] * _feUdy[i]) * (v[0] * dudx + v[1] * dudy);
    // }
    double dx01 = fabs(form->_geoCoord[3]-form->_geoCoord[0]);
    double dx02 = fabs(form->_geoCoord[6]-form->_geoCoord[0]);
    double dx12 = fabs(form->_geoCoord[6]-form->_geoCoord[3]);
    double dy01 = fabs(form->_geoCoord[4]-form->_geoCoord[1]);
    double dy02 = fabs(form->_geoCoord[7]-form->_geoCoord[1]);
    double dy12 = fabs(form->_geoCoord[7]-form->_geoCoord[4]);
    double h01 = sqrt(dx01*dx01 + dy01*dy01);
    double h02 = sqrt(dx02*dx02 + dy02*dy02);
    double h12 = sqrt(dx12*dx12 + dy12*dy12);
    double h = fmax(h12,fmax(h01,h02));
    double normeU = sqrt(v[0]*v[0]+v[1]*v[1]);
    // double Peh = normeU * h / 2.0;

    // double tau = 1.0/( sqrt( 4.0/(form->_dt*form->_dt) + 4.0*normeU*normeU/h/h + 9.0*16.0/(h*h*h*h)) );
    double tau = 1.0/( sqrt( 4.0/(form->_dt*form->_dt) + 4.0*normeU*normeU/h/h ) );
    // tau = 0.01;
    // double tau = h/2.0/sqrt(v[0]*v[0]+v[1]*v[1]) * (cosh(Peh)/sinh(Peh) - 1.0/Peh);
    // double tau = (cosh(Peh)/sinh(Peh) - 1.0/Peh);
    // double tau = 0.01;
    // std::cout<<tau<<std::endl;
    // double deltaSUPG = c/kT;
    // double Re = (v[0]*v[0]+v[1]*v[1]) * c/kT;
    // printf("Re = %+-4.4e - c = %+-4.4e - cT = %+-4.4e - kT = %+-4.4e - delta = %+-4.4e\n", Re, c, cT, kT, deltaSUPG);

    for(int i = 0; i < nFunctions; ++i) {
      u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
      dudt = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_solDot[_idU],k);
      dudx = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdx +
             form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdx;
      dudy = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdy +
             form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdy;

      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
      // form->_Be[i] -= (v[0] * _feUdx[i] + v[1] * _feUdy[i]) * u * J * w[k];
      // SUPG ?
      // double delta = form->_geoCoord[3]-form->_geoCoord[0];
      // double deltaSUPG = 0.01;
      // form->_Be[i] -= ((v[0] * dudx + v[1] * dudy) * _feU[i] + delta * (v[0] * dudx + v[1] * dudy) * (v[0] * _feUdx[i] + v[1] * _feUdy[i])) * J * w[k];
      // form->_Be[i] -= Jac * w[k] * ((v[0] * dudx + v[1] * dudy) * _feU[i] + deltaSUPG * (v[0] * _feUdx[i] + v[1] * _feUdy[i]) * (dudt + v[0] * dudx + v[1] * dudy));
      form->_Be[i] -= Jac * w[k] * ((v[0] * dudx + v[1] * dudy) * _feU[i] + tau       * (v[0] * _feUdx[i] + v[1] * _feUdy[i]) * (dudt + v[0] * dudx + v[1] * dudy));
    }
  }
}

// -----------------------------------------------------------------------------
// 2D SUPG Stabilization
// -----------------------------------------------------------------------------
void feSysElm_2D_SUPGStab::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _phiU.resize(_nFunctions);
  _gradPhiU.resize(2 * _nFunctions);
}

void feSysElm_2D_SUPGStab::computeAe(feBilinearForm *form)
{
  // Finite differences for now
}

static double tau2DTri(int dim, double velocity[3], double diffusivity, int nFunctions, std::vector<double> &gradphi)
{
  // Compute hTau from (13.17) in Fortin & Garon
  double normV = sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]);

  double res = 0., dotprod;
  for(int i = 0; i < nFunctions; ++i){
    dotprod = 0.;
    for(int iDim = 0; iDim < dim; ++iDim){
      dotprod += velocity[iDim]/normV * gradphi[i*dim + iDim];
    }
    res += dotprod*dotprod;
  }

  double hTau = sqrt(2.) / sqrt(res);

  double rho = 1.;
  double cp = 1.;
  return 1.0 / sqrt( 4.*normV*normV/(hTau*hTau) + 144.*diffusivity*diffusivity/((rho*cp*hTau*hTau)*(rho*cp*hTau*hTau)) );
}

void feSysElm_2D_SUPGStab::computeBe(feBilinearForm *form)
{
  double jac, u, grad_u[3] = {0., 0., 0.}, tau, diffusivity, source, reaction;
  std::vector<double> c(2, 0.0);
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    // Evaluate the imposed velocity field, diffusivity, reaction coefficient and source term
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    (*_velocity)(form->_tn, _pos, c);
    diffusivity = (*_diffusivity)(form->_tn, _pos);
    reaction = (*_reactionCoeff)(form->_tn, _pos);
    source = (*_source)(form->_tn, _pos);

    // Compute u
    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);

    // Compute grad(u)
    form->_intSpaces[_idU]->interpolateFieldAtQuadNode_physicalGradient(form->_sol[_idU], k, jac, form->_geoCoord, 
      form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, grad_u);

    // Compute grad(phi)
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, 
      form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, _gradPhiU.data());

    // Compute stabilization parameter tau
    tau = tau2DTri(2, c.data(), diffusivity, _nFunctions, _gradPhiU);
    
    // Compute PDE residual:
    // Reaction term
    double residual = reaction * u;
    // Convection term 
    for(int iDim = 0; iDim < 2; ++iDim){
      residual += c[iDim] * grad_u[iDim];
    }
    // Diffusion term (0 for linear elements)
    // ...
    /// Source term
    residual += source;

    for(int i = 0; i < _nFunctions; ++i) {
      for(int iDim = 0; iDim < 2; ++iDim){
        form->_Be[i] -= residual * tau * c[iDim] * _gradPhiU[i*2 + iDim] * jac * _wQuad[k];
      }
    }
  }
}

void feSysElm_2D_Stokes::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idV = 1;
  _idP = 2;
  _fieldsLayoutI = {_idU, _idV, _idP};
  _fieldsLayoutJ = {_idU, _idV, _idP};
  _feU.resize(space[_idU]->getNumFunctions());
  _feUdx.resize(space[_idU]->getNumFunctions());
  _feUdy.resize(space[_idU]->getNumFunctions());
  _feV.resize(space[_idV]->getNumFunctions());
  _feVdx.resize(space[_idV]->getNumFunctions());
  _feVdy.resize(space[_idV]->getNumFunctions());
  _feP.resize(space[_idP]->getNumFunctions());
  _dxdr.resize(3); // [dx/dr, dy/dr, dz/dr]
  _dxds.resize(3); // [dx/ds, dy/ds, dz/ds]
}

void feSysElm_2D_Stokes::computeBe(feBilinearForm *form)
{
  // TODO : Verifier que les règles d'intégration soient identiques ou adapter
  // On prend les poids de form->_geoSpace, ok pour le jacobien seulement ?
  int nG = form->_geoSpace->getNumQuadPoints();
  std::vector<double> &w = form->_geoSpace->getQuadratureWeights();
  // double             rho = _par[0];
  double mu = _par[1];
  int nFunctionsU = form->_intSpaces[_idU]->getNumFunctions();
  int nFunctionsV = form->_intSpaces[_idV]->getNumFunctions();
  int nFunctionsP = form->_intSpaces[_idP]->getNumFunctions();

  double J, u, v, p, dudx, dudy, dvdx, dvdy, Sxx, Sxy, Syx, Syy;
  std::vector<double> x(3);
  std::vector<double> f(2); // Forces volumiques

  for(int k = 0; k < nG; ++k) {
    // Forces volumiques
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, x);
    if(_fct != nullptr) _fct->eval(form->_tn, x, f);

    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, _dxdr);
    form->_geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(form->_geoCoord, k, _dxds);
    J = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    double drdx = _dxds[1] / J;
    double drdy = -_dxds[0] / J;
    double dsdx = -_dxdr[1] / J;
    double dsdy = _dxdr[0] / J;

    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    v = form->_intSpaces[_idV]->interpolateFieldAtQuadNode(form->_sol[_idV], k);
    p = form->_intSpaces[_idP]->interpolateFieldAtQuadNode(form->_sol[_idP], k);

    dudx = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdx +
           form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdx;
    dudy = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdy +
           form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdy;
    dvdx = form->_intSpaces[_idV]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idV], k) * drdx +
           form->_intSpaces[_idV]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idV], k) * dsdx;
    dvdy = form->_intSpaces[_idV]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idV], k) * drdy +
           form->_intSpaces[_idV]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idV], k) * dsdy;

    Sxx = -p + 2. * mu * dudx;
    Sxy = mu * (dudy + dvdx);
    Syx = mu * (dudy + dvdx);
    Syy = -p + 2. * mu * dvdy;

    int cnt = 0;
    // Residu pour u
    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(
        i, k); // Fonction test de u : uniquement pour les forces volumiques
      _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
      form->_Be[cnt++] -= (Sxx * _feUdx[i] + Sxy * _feUdy[i] - f[0] * _feU[i]) * J * w[k];
    }
    // Residu pour v
    for(int i = 0; i < nFunctionsV; ++i) {
      _feV[i] = form->_intSpaces[_idV]->getFunctionAtQuadNode(
        i, k); // Fonction test de v : uniquement pour les forces volumiques
      _feVdx[i] = form->_intSpaces[_idV]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idV]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feVdy[i] = form->_intSpaces[_idV]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idV]->getdFunctiondsAtQuadNode(i, k) * dsdy;
      form->_Be[cnt++] -= (Syx * _feVdx[i] + Syy * _feVdy[i] - f[1] * _feV[i]) * J * w[k];
    }
    // Residu pour p
    for(int i = 0; i < nFunctionsP; ++i) {
      _feP[i] = form->_intSpaces[_idP]->getFunctionAtQuadNode(i, k);
      form->_Be[cnt++] += _feP[i] * (dudx + dvdy) * J * w[k];
    }
  }
}

void feSysElm_2D_Stokes::computeAe(feBilinearForm *form)
{
  // TODO : Verifier que les règles d'intégration soient identiques ou adapter
  // On prend les poids de form->_geoSpace, ok pour le jacobien seulement ?
  int nG = form->_geoSpace->getNumQuadPoints();
  std::vector<double> &w = form->_geoSpace->getQuadratureWeights();
  // double             rho = _par[0];
  double mu = _par[1];
  int nFunctionsU = form->_intSpaces[_idU]->getNumFunctions();
  int nFunctionsV = form->_intSpaces[_idV]->getNumFunctions();
  int nFunctionsP = form->_intSpaces[_idP]->getNumFunctions();

  double jac, dudx, dudy, dvdx, dvdy;
  std::vector<double> x(3);
  std::vector<double> f(2, 0.); // Body forces

  for(int k = 0; k < nG; ++k) {
    // Body forces
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, x);
    if(_fct != nullptr) _fct->eval(form->_tn, x, f);

    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, _dxdr);
    form->_geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(form->_geoCoord, k, _dxds);
    jac = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    double drdx = _dxds[1] / jac;
    double drdy = -_dxds[0] / jac;
    double dsdx = -_dxdr[1] / jac;
    double dsdy = _dxdr[0] / jac;

    dudx = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdx +
           form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdx;
    dudy = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdy +
           form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdy;
    dvdx = form->_intSpaces[_idV]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idV], k) * drdx +
           form->_intSpaces[_idV]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idV], k) * dsdx;
    dvdy = form->_intSpaces[_idV]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idV], k) * drdy +
           form->_intSpaces[_idV]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idV], k) * dsdy;

    int I = 0, J;

    // Set up interpolation functions
    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
    }
    for(int i = 0; i < nFunctionsV; ++i) {
      _feV[i] = form->_intSpaces[_idV]->getFunctionAtQuadNode(i, k);
      _feVdx[i] = form->_intSpaces[_idV]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feVdy[i] = form->_intSpaces[_idV]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
    }
    for(int i = 0; i < nFunctionsP; ++i) {
      _feP[i] = form->_intSpaces[_idP]->getFunctionAtQuadNode(i, k);
    }

    // Équations pour u
    for(int i = 0; i < nFunctionsU; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) {
        form->_Ae[I][J++] += jac * w[k] * (_feUdx[i] * mu * 2. * _feUdx[j] + _feUdy[i] * mu * _feUdy[j]);
      }
      for(int j = 0; j < nFunctionsV; ++j) {
        form->_Ae[I][J++] += jac * w[k] * (_feUdy[i] * mu * _feVdx[j]);
      }
      for(int j = 0; j < nFunctionsP; ++j) {
        form->_Ae[I][J++] -= jac * w[k] * (_feUdx[i] * _feP[j]);
      }
      ++I;
    }

    // Équations pour v
    for(int i = 0; i < nFunctionsV; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) {
        form->_Ae[I][J++] += jac * w[k] * (_feVdx[i] * mu * _feUdy[j]);
      }
      for(int j = 0; j < nFunctionsV; ++j) {
        form->_Ae[I][J++] += jac * w[k] * (_feVdy[i] * mu * 2. * _feVdy[j] + _feVdx[i] * mu * _feVdx[j]);
      }
      for(int j = 0; j < nFunctionsP; ++j) {
        form->_Ae[I][J++] -= jac * w[k] * (_feVdy[i] * _feP[j]);
      }
      ++I;
    }

    // Équations pour p
    for(int i = 0; i < nFunctionsP; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) {
        form->_Ae[I][J++] -= jac * w[k] * (_feP[i] * _feUdx[j]);
      }
      for(int j = 0; j < nFunctionsV; ++j) {
        form->_Ae[I][J++] -= jac * w[k] * (_feP[i] * _feVdy[j]);
      }
      for(int j = 0; j < nFunctionsP; ++j) {
        form->_Ae[I][J++] += 0.0;
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
  _fieldsLayoutI = {_idU, _idV, _idP};
  _fieldsLayoutJ = {_idU, _idV, _idP};
  _feU.resize(space[_idU]->getNumFunctions());
  _feUdx.resize(space[_idU]->getNumFunctions());
  _feUdy.resize(space[_idU]->getNumFunctions());
  _feV.resize(space[_idV]->getNumFunctions());
  _feVdx.resize(space[_idV]->getNumFunctions());
  _feVdy.resize(space[_idV]->getNumFunctions());
  _feP.resize(space[_idP]->getNumFunctions());
  _dxdr.resize(3); // [dx/dr, dy/dr, dz/dr]
  _dxds.resize(3); // [dx/ds, dy/ds, dz/ds]
}

void feSysElm_2D_NavierStokes::computeBe(feBilinearForm *form)
{
  // TODO : Verifier que les règles d'intégration soient identiques ou adapter
  // On prend les poids de form->_geoSpace, ok pour le jacobien seulement ?
  int nG = form->_geoSpace->getNumQuadPoints();
  std::vector<double> &w = form->_geoSpace->getQuadratureWeights();
  double rho = _par[0];
  double mu = _par[1];
  int nFunctionsU = form->_intSpaces[_idU]->getNumFunctions();
  int nFunctionsV = form->_intSpaces[_idV]->getNumFunctions();
  int nFunctionsP = form->_intSpaces[_idP]->getNumFunctions();

  double J, u, v, p, dudt, dvdt, dudx, dudy, dvdx, dvdy, Sxx, Sxy, Syx, Syy;
  std::vector<double> x(3, 0.);
  std::vector<double> f(2, 0.); // Forces volumiques

  for(int k = 0; k < nG; ++k) {
    // Forces volumiques
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, x);
    if(_fct != nullptr) _fct->eval(form->_tn, x, f);

    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, _dxdr);
    form->_geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(form->_geoCoord, k, _dxds);
    J = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    double drdx = _dxds[1] / J;
    double drdy = -_dxds[0] / J;
    double dsdx = -_dxdr[1] / J;
    double dsdy = _dxdr[0] / J;

    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    v = form->_intSpaces[_idV]->interpolateFieldAtQuadNode(form->_sol[_idV], k);
    p = form->_intSpaces[_idP]->interpolateFieldAtQuadNode(form->_sol[_idP], k);

    dudt = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_solDot[_idU], k);
    dvdt = form->_intSpaces[_idV]->interpolateFieldAtQuadNode(form->_solDot[_idV], k);

    dudx = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdx +
           form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdx;
    dudy = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdy +
           form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdy;
    dvdx = form->_intSpaces[_idV]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idV], k) * drdx +
           form->_intSpaces[_idV]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idV], k) * dsdx;
    dvdy = form->_intSpaces[_idV]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idV], k) * drdy +
           form->_intSpaces[_idV]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idV], k) * dsdy;

    Sxx = -p + 2. * mu * dudx;
    Sxy = mu * (dudy + dvdx);
    Syx = mu * (dudy + dvdx);
    Syy = -p + 2. * mu * dvdy;

    int cnt = 0;
    // Residu pour u
    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
      form->_Be[cnt++] -= (_feU[i] * rho * (dudt + u * dudx + v * dudy) + Sxx * _feUdx[i] +
                    Sxy * _feUdy[i] - f[0] * _feU[i]) *
                   J * w[k];
    }
    // Residu pour v
    for(int i = 0; i < nFunctionsV; ++i) {
      _feV[i] = form->_intSpaces[_idV]->getFunctionAtQuadNode(i, k);
      _feVdx[i] = form->_intSpaces[_idV]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idV]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feVdy[i] = form->_intSpaces[_idV]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idV]->getdFunctiondsAtQuadNode(i, k) * dsdy;
      form->_Be[cnt++] -= (_feV[i] * rho * (dvdt + u * dvdx + v * dvdy) + Syx * _feVdx[i] +
                    Syy * _feVdy[i] - f[1] * _feV[i]) *
                   J * w[k];
    }
    // Residu pour p
    for(int i = 0; i < nFunctionsP; ++i) {
      _feP[i] = form->_intSpaces[_idP]->getFunctionAtQuadNode(i, k);
      form->_Be[cnt++] += _feP[i] * (dudx + dvdy) * J * w[k];
    }
  }
}

void feSysElm_2D_NavierStokes::computeAe(feBilinearForm *form)
{
  // TODO : Verifier que les règles d'intégration soient identiques ou adapter
  // On prend les poids de form->_geoSpace, ok pour le jacobien seulement ?
  int nG = form->_geoSpace->getNumQuadPoints();
  std::vector<double> w = form->_geoSpace->getQuadratureWeights();
  double rho = _par[0];
  double mu = _par[1];
  int nFunctionsU = form->_intSpaces[_idU]->getNumFunctions();
  int nFunctionsV = form->_intSpaces[_idV]->getNumFunctions();
  int nFunctionsP = form->_intSpaces[_idP]->getNumFunctions();

  double jac, u, v, dudt, dvdt, dudx, dudy, dvdx, dvdy;
  std::vector<double> x(3);
  std::vector<double> f(2, 0.); // Forces volumiques

  for(int k = 0; k < nG; ++k) {
    // Forces volumiques
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, x);
    if(_fct != nullptr) _fct->eval(form->_tn, x, f);

    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, _dxdr);
    form->_geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(form->_geoCoord, k, _dxds);
    jac = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    double drdx = _dxds[1] / jac;
    double drdy = -_dxds[0] / jac;
    double dsdx = -_dxdr[1] / jac;
    double dsdy = _dxdr[0] / jac;

    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    v = form->_intSpaces[_idV]->interpolateFieldAtQuadNode(form->_sol[_idV], k);

    dudt = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_solDot[_idU], k);
    dvdt = form->_intSpaces[_idV]->interpolateFieldAtQuadNode(form->_solDot[_idV], k);

    dudx = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdx +
           form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdx;
    dudy = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdy +
           form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdy;
    dvdx = form->_intSpaces[_idV]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idV], k) * drdx +
           form->_intSpaces[_idV]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idV], k) * dsdx;
    dvdy = form->_intSpaces[_idV]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idV], k) * drdy +
           form->_intSpaces[_idV]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idV], k) * dsdy;

    int I = 0, J = 0;

    // Set up interpolation functions
    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
    }
    for(int i = 0; i < nFunctionsV; ++i) {
      _feV[i] = form->_intSpaces[_idV]->getFunctionAtQuadNode(i, k);
      _feVdx[i] = form->_intSpaces[_idV]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feVdy[i] = form->_intSpaces[_idV]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
    }
    for(int i = 0; i < nFunctionsP; ++i) {
      _feP[i] = form->_intSpaces[_idP]->getFunctionAtQuadNode(i, k);
    }

    // Équations pour u
    for(int i = 0; i < nFunctionsU; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) {
        form->_Ae[I][J++] +=
          jac * w[k] *
          (rho * _feU[i] * (form->_c0 * _feU[j] + u * _feUdx[j] + v * _feUdy[j] + _feU[j] * dudx) +
           _feUdx[i] * mu * 2. * _feUdx[j] + _feUdy[i] * mu * _feUdy[j]);
      }
      for(int j = 0; j < nFunctionsV; ++j) {
        form->_Ae[I][J++] += jac * w[k] * (rho * _feU[i] * (_feV[j] * dudy) + _feUdy[i] * mu * _feVdx[j]);
      }
      for(int j = 0; j < nFunctionsP; ++j) {
        form->_Ae[I][J++] -= jac * w[k] * (_feUdx[i] * _feP[j]);
      }
      ++I;
    }

    // Équations pour v
    for(int i = 0; i < nFunctionsV; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) {
        form->_Ae[I][J++] += jac * w[k] * (rho * _feV[i] * (_feU[j] * dvdx) + _feVdx[i] * mu * _feUdy[j]);
      }
      for(int j = 0; j < nFunctionsV; ++j) {
        form->_Ae[I][J++] +=
          jac * w[k] *
          (rho * _feV[i] * (form->_c0 * _feV[j] + u * _feVdx[j] + v * _feVdy[j] + _feV[j] * dvdy) +
           _feVdy[i] * mu * 2. * _feVdy[j] + _feVdx[i] * mu * _feVdx[j]);
      }
      for(int j = 0; j < nFunctionsP; ++j) {
        form->_Ae[I][J++] -= jac * w[k] * (_feVdy[i] * _feP[j]);
      }
      ++I;
    }

    // Équations pour p
    for(int i = 0; i < nFunctionsP; ++i) {
      J = 0;
      for(int j = 0; j < nFunctionsU; ++j) {
        form->_Ae[I][J++] -= jac * w[k] * (_feP[i] * _feUdx[j]);
      }
      for(int j = 0; j < nFunctionsV; ++j) {
        form->_Ae[I][J++] -= jac * w[k] * (_feP[i] * _feVdy[j]);
      }
      for(int j = 0; j < nFunctionsP; ++j) {
        form->_Ae[I][J++] += 0.0;
      }
      ++I;
    }
  }
}

void feSysElm_2D_Euler::createElementarySystem(std::vector<feSpace *> &space)
{
  _idr  = 0;
  _idru = 1;
  _idrv = 2;
  _idre = 3;
  _fieldsLayoutI = {_idr, _idru, _idrv, _idre};
  _fieldsLayoutJ = {_idr, _idru, _idrv, _idre};
  _nFunctionsr = space[_idr]->getNumFunctions();
  _nFunctionsru = space[_idru]->getNumFunctions();
  _nFunctionsrv = space[_idrv]->getNumFunctions();
  _nFunctionsre = space[_idre]->getNumFunctions();
  _phir.resize(_nFunctionsr);
  _gradphir.resize(2 * _nFunctionsr);
  _phiru.resize(_nFunctionsru);
  _gradphiru.resize(2 * _nFunctionsru);
  _phirv.resize(_nFunctionsrv);
  _gradphirv.resize(2 * _nFunctionsrv);
  _phire.resize(_nFunctionsre);
  _gradphire.resize(2 * _nFunctionsre);
}

void feSysElm_2D_Euler::computeAe(feBilinearForm *form)
{
  // Finite differences
}

void feSysElm_2D_Euler::computeBe(feBilinearForm *form)
{
  double gamma = 1.4;
  double jac, r, u, ru, v, rv, re, p, T;

  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    
    // Compute current fields and gradient of test functions
     r = form->_intSpaces[_idr ]->interpolateFieldAtQuadNode(form->_sol[_idr ], k);
    ru = form->_intSpaces[_idru]->interpolateFieldAtQuadNode(form->_sol[_idru], k);
    rv = form->_intSpaces[_idrv]->interpolateFieldAtQuadNode(form->_sol[_idrv], k);
    re = form->_intSpaces[_idre]->interpolateFieldAtQuadNode(form->_sol[_idre], k);
    u = ru/r;
    v = rv/r;
    T = (re/r) - (u*u+v*v)/2.;
    p = (gamma - 1.) * r * T;

    form->_intSpaces[_idr]->getFunctionsAtQuadNode(k, _phir);
    form->_intSpaces[_idru]->getFunctionsAtQuadNode(k, _phiru);
    form->_intSpaces[_idrv]->getFunctionsAtQuadNode(k, _phirv);
    form->_intSpaces[_idre]->getFunctionsAtQuadNode(k, _phire);

    form->_intSpaces[_idr]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, 
      form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, _gradphir.data());
    form->_intSpaces[_idru]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, 
      form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, _gradphiru.data());
    form->_intSpaces[_idrv]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, 
      form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, _gradphirv.data());
    form->_intSpaces[_idre]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, 
      form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, _gradphire.data());

     // + 1.0 * _phiru[i]
     // + 1.0 * _phirv[i]

    int cnt = 0;
    // Residu pour rho
    for(int i = 0; i < _nFunctionsr; ++i) {
      form->_Be[cnt++] += (ru * _gradphir[2*i+0] + rv * _gradphir[2*i+1]) * jac * _wQuad[k];
    }
    // Residu pour rho*u
    for(int i = 0; i < _nFunctionsru; ++i) {
      form->_Be[cnt++] += ((ru*u+p) * _gradphiru[2*i+0] + ru*v * _gradphiru[2*i+1]) * jac * _wQuad[k];
    }
    // Residu pour rho*v
    for(int i = 0; i < _nFunctionsrv; ++i) {
      form->_Be[cnt++] += (ru*v * _gradphirv[2*i+0] + (rv*v+p) * _gradphirv[2*i+1]) * jac * _wQuad[k];
    }
    // Residu pour rho*e
    for(int i = 0; i < _nFunctionsre; ++i) {
      form->_Be[cnt++] += (u*(re + p) * _gradphire[2*i+0] + v*(re + p) * _gradphire[2*i+1]) * jac * _wQuad[k];
    }
  }

  // for(int i = 0; i < _nFunctionsr + _nFunctionsru + _nFunctionsrv + _nFunctionsre; ++i)
  //   feInfo("Computed %f", form->_Be[i]);
}

void feSysElm_1D_EulerBord::createElementarySystem(std::vector<feSpace *> &space)
{
  _idr  = 0;
  _idru = 1;
  _idrv = 2;
  _idre = 3;
  _fieldsLayoutI = {_idr, _idru, _idrv, _idre};
  _fieldsLayoutJ = {_idr, _idru, _idrv, _idre};
  _nFunctionsr = space[_idr]->getNumFunctions();
  _nFunctionsru = space[_idru]->getNumFunctions();
  _nFunctionsrv = space[_idrv]->getNumFunctions();
  _nFunctionsre = space[_idre]->getNumFunctions();
  _phir.resize(_nFunctionsr);
  _gradphir.resize(2 * _nFunctionsr);
  _phiru.resize(_nFunctionsru);
  _gradphiru.resize(2 * _nFunctionsru);
  _phirv.resize(_nFunctionsrv);
  _gradphirv.resize(2 * _nFunctionsrv);
  _phire.resize(_nFunctionsre);
  _gradphire.resize(2 * _nFunctionsre);
}

void feSysElm_1D_EulerBord::computeBe(feBilinearForm *form)
{
  double gamma = 1.4;
  double jac, r, u, ru, v, rv, re, p, T;
  double Fn[4] = {0., 0., 0., 0.};

  // std::string name = "elem" + std::to_string(form->_numElem) + ".pos";
  // FILE *f = fopen(name.c_str(), "w");
  // fprintf(f, "View \"elem\" {\n");

  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    
    // Compute current fields and gradient of test functions
     r = form->_intSpaces[_idr ]->interpolateFieldAtQuadNode(form->_sol[_idr ], k);
    ru = form->_intSpaces[_idru]->interpolateFieldAtQuadNode(form->_sol[_idru], k);
    rv = form->_intSpaces[_idrv]->interpolateFieldAtQuadNode(form->_sol[_idrv], k);
    re = form->_intSpaces[_idre]->interpolateFieldAtQuadNode(form->_sol[_idre], k);
    u = ru/r;
    v = rv/r;
    T = (re/r) - (u*u+v*v)/2.;
    p = (gamma - 1.) * r * T;

    form->_intSpaces[_idr]->getFunctionsAtQuadNode(k, _phir);
    form->_intSpaces[_idru]->getFunctionsAtQuadNode(k, _phiru);
    form->_intSpaces[_idrv]->getFunctionsAtQuadNode(k, _phirv);
    form->_intSpaces[_idre]->getFunctionsAtQuadNode(k, _phire);

    // form->_intSpaces[_idr]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, 
    //   form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, _gradphir.data());
    // form->_intSpaces[_idru]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, 
    //   form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, _gradphiru.data());
    // form->_intSpaces[_idrv]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, 
    //   form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, _gradphirv.data());
    // form->_intSpaces[_idre]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, 
    //   form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, _gradphire.data());

    double x0 = form->_geoCoord[0];
    double y0 = form->_geoCoord[1];
    double x1 = form->_geoCoord[3];
    double y1 = form->_geoCoord[4];
    double nx = y1-y0;
    double ny = x0-x1;
    double N = sqrt(nx*nx + ny*ny);
    nx /= N;
    ny /= N;
    if(fabs(N) < 1e-3){
      feInfo("ICI");
      exit(-1);
    }

    // fprintf(f, "VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n", _pos[0], _pos[1], 0., nx, ny, 0.0);

    if(_boundaryType == 1)
    {
      // Wall
      Fn[0] = 0.;
      Fn[1] = p * nx;
      Fn[2] = p * ny;
      Fn[3] = 0.;
    }
    else if(_boundaryType == 2)
    {
      // // Inflow: imposed u, v, p and T ?
      // double R = 1.;
      // double Ti = 1.;
      // double pi = 1.;
      // double rhoi = pi/(R*Ti);
      // double ui = 0.01;
      // double vi = 0.;
      // double Ei = Ti + (ui*ui+vi*vi)/2.;
      // double uin = ui*nx + vi*ny;
      // // Fn[0] = rhoi*un;
      // // Fn[1] = rhoi*ui*un + pi*nx;
      // // Fn[2] = rhoi*vi*un + pi*ny;
      // // Fn[3] = (rhoi*Ei+pi)*un;
      // // feInfo("Inlet flux: %f - %f - %f - %f", Fn[0], Fn[1], Fn[2], Fn[3]);

      // double un = u*nx + v*ny;
      // double cin = gamma * p/r;
      // double cinf = gamma * pi/rhoi;
      // double Rplus  = un + 2*cin/(gamma-1.);
      // double Rminus = uin + 2*cinf/(gamma-1.);
      // double Vb = (Rplus + Rminus)/2.;
      // double cb = (gamma-1.)/4. * (Rplus - Rminus);
      // double ub = ui + (Vb - uin)*nx;
      // double vb = vi + (Vb - uin)*ny;
      // double sb = cinf*cinf / (gamma * pow(rhoi, gamma-1.));
      // double rhob = cb*cb/(gamma * sb);
      // double pb = rhob*cb*cb/gamma;
      // double ubn = ub*nx + vb*ny;
      // double Eb = pb/(gamma-1.) + (ub*ub+vb*vb)/2.;
      // Fn[0] = rhob*ubn;
      // Fn[1] = rhob*ub*ubn + pb*nx;
      // Fn[2] = rhob*vb*ubn + pb*ny;
      // Fn[3] = (Eb+pb)*ubn;
      double pt_set = 102100.;
      double Tt_set = 288.6;
      double R = 287.;

      double Ht = p/r * gamma/(gamma-1.) + 0.5*(u*u+v*v);
      double ci = sqrt(gamma*p/r);
      double Rplus = -(u*nx+v*ny) - 2.*ci/(gamma-1.);
      double a = 1. + 2./(gamma-1.);
      double b = 2*Rplus;
      double c = (gamma-1.)/2. * (Rplus*Rplus - 2*Ht);
      // feInfo("%f - %f - %f", gamma, p, r);
      // feInfo("%f - %f - %f", ci, u, v);
      // feInfo("%f - %f - %f", r, Ht, Rplus);
      // feInfo("%f - %f - %f", b, c, b*b - 4.*a*c);
      double cb1 = -b/(2.*a) + sqrt(b*b - 4.*a*c)/(2.*a);
      double cb2 = -b/(2.*a) - sqrt(b*b - 4.*a*c)/(2.*a);
      if(isnan(cb1) || isnan(cb2)){
        feInfo("%f - %f", cb1, cb2);
        exit(-1);
      }
      double cb = fmax(cb1, cb2);
      double U = 2*cb/(gamma-1.);
      double Mb = U/cb;
      double pb = pt_set * pow((1. + (gamma-1.)/2. * Mb*Mb), gamma/(gamma-1.));
      double Tb = Tt_set * pow( pb/pt_set , (gamma-1.)/gamma);
      double rhob = pb/(R*Tb);
      double ub = U*nx;
      double vb = U*ny;
      double Eb = pb/(gamma-1.) + (ub*ub+vb*vb)/2.;
      double ubn = ub*nx + vb*ny;
      Fn[0] = rhob*ubn;
      Fn[1] = rhob*ub*ubn + pb*nx;
      Fn[2] = rhob*vb*ubn + pb*ny;
      Fn[3] = (Eb+pb)*ubn;
    }
    else if(_boundaryType == 3)
    {
      // Outflow: imposed p ?
      // double po = 0.98;
      // // double un = u*nx + v*ny;
      // // Fn[0] = r*un;
      // // Fn[1] = r*u*un + po*nx;
      // // Fn[2] = r*v*un + po*ny;
      // // Fn[3] = (re+po)*un;
      // double R = 1.;
      // double Ti = 1.;
      // double pi = 1.;
      // double rhoi = pi/(R*Ti);
      // double ui = 0.01;
      // double vi = 0.;
      // double uin = ui*nx + vi*ny;
      // double un = u*nx + v*ny;
      // double cin = gamma * p/r;
      // double cinf = gamma * po/rhoi;
      // double Rplus  = un + 2*cin/(gamma-1.);
      // double Rminus = uin + 2*cinf/(gamma-1.);

      // double Vb = (Rplus + Rminus)/2.;
      // double cb = (gamma-1.)/4. * (Rplus - Rminus);
      // double ub = u + (Vb - un)*nx;
      // double vb = v + (Vb - un)*ny;
      // double sb = cinf*cinf / (gamma * pow(rhoi, gamma-1.));
      // double rhob = cb*cb/(gamma * sb);
      // double pb = rhob*cb*cb/gamma;
      // double ubn = ub*nx + vb*ny;
      // double Eb = pb/(gamma-1.) + (ub*ub+vb*vb)/2.;
      // Fn[0] = rhob*ubn;
      // Fn[1] = rhob*ub*ubn + pb*nx;
      // Fn[2] = rhob*vb*ubn + pb*ny;
      // Fn[3] = (rhob*Eb+pb)*ubn;

      // Inflow/Outflow Boundary Conditions with Application to FUN3D, Jan-Renee Carlson
      double p_out = 1.; // = p_set
      double pb = p_out;
      double ub = u;
      double vb = v;
      // double Ti = gamma * p / r; // Definition?
      double Ti = (gamma - 1.) * p / r; // Definition?
      double rhob = gamma * pb / Ti;
      double ubn = ub*nx + vb*ny;
      double Eb = pb/(gamma-1.) + (ub*ub+vb*vb)/2.;
      // The state vector at the boundary is qb = (rhob, ub, vb, pb)
      // The flux F(qb) is:
      Fn[0] = rhob*ubn;
      Fn[1] = rhob*ub*ubn + pb*nx;
      Fn[2] = rhob*vb*ubn + pb*ny;
      Fn[3] = (Eb+pb)*ubn;
    }
  
    int cnt = 0;
    // Residu pour rho
    for(int i = 0; i < _nFunctionsr; ++i) {
      form->_Be[cnt++] -= Fn[0] * _phir[i] * jac * _wQuad[k];
    }
    // Residu pour rho*u
    for(int i = 0; i < _nFunctionsru; ++i) {
      form->_Be[cnt++] -= Fn[1] * _phiru[i] * jac * _wQuad[k];
    }
    // Residu pour rho*v
    for(int i = 0; i < _nFunctionsrv; ++i) {
      form->_Be[cnt++] -= Fn[2] * _phirv[i] * jac * _wQuad[k];
    }
    // Residu pour rho*e
    for(int i = 0; i < _nFunctionsre; ++i) {
      form->_Be[cnt++] -= Fn[3] * _phire[i] * jac * _wQuad[k];
    }
  }

  // for(int i = 0; i < _nFunctionsr + _nFunctionsru + _nFunctionsrv + _nFunctionsre; ++i)
  //   feInfo("Computed %f", form->_Be[i]);
}