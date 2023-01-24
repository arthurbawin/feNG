#include "feSysElm.h"
#include "feBilinearForm.h"


// -----------------------------------------------------------------------------
// 2D Source
// -----------------------------------------------------------------------------
void feSysElm_2D_Source::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _feU.resize(_nFunctions);
}

void feSysElm_2D_Source::computeBe(feBilinearForm *form)
{
  for(int k = 0; k < _nQuad; ++k) {
    double jac = form->_J[_nQuad * form->_numElem + k];
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    double S = _sourceFun->eval(form->_tn, _pos);
    for(int i = 0; i < _nFunctions; ++i) {
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= _feU[i] * S * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// 2D Diffusion (stiffness matrix)
// -----------------------------------------------------------------------------
void feSysElm_2D_Diffusion::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _feU.resize(_nFunctions);
  _feUdx.resize(_nFunctions);
  _feUdy.resize(_nFunctions);
  _dxdr.resize(3);
  _dxds.resize(3);
}

void feSysElm_2D_Diffusion::computeAe(feBilinearForm *form)
{
  double jac;
  for(int k = 0; k < _nQuad; ++k) {
    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, _dxdr);
    form->_geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(form->_geoCoord, k, _dxds);

    jac = form->_J[_nQuad * form->_numElem + k];
    double drdx = _dxds[1] / jac;
    double drdy = -_dxds[0] / jac;
    double dsdx = -_dxdr[1] / jac;
    double dsdy = _dxdr[0] / jac;

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    double kD = (*_diffusivity)(form->_tn, _pos);

    for(int i = 0; i < _nFunctions; ++i) {
      _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
    }

    // double dphidx[3], dphidy[3], dphidz[3];

    // form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, dphidx, dphidy, dphidz);

    // for(int ii = 0; ii < 3; ++ii){
    //   feInfo("Computed gradient dphidx[%d] = %f vs _feudx[%d] = %f", ii, dphidx[ii], ii, _feUdx[ii]);
    //   feInfo("Computed gradient dphidy[%d] = %f vs _feudy[%d] = %f", ii, dphidy[ii], ii, _feUdy[ii]);
    //   feInfo("Computed gradient dphidz[%d] = %f", ii, dphidz[ii]);
    // }

    // exit(-1);

    for(int i = 0; i < _nFunctions; ++i) {
      for(int j = 0; j < _nFunctions; ++j) {
        form->_Ae[i][j] += (_feUdx[i] * _feUdx[j] + _feUdy[i] * _feUdy[j]) * kD * jac * _wQuad[k];
      }
    }
  }
}

void feSysElm_2D_Diffusion::computeBe(feBilinearForm *form)
{
  double jac, dudx, dudy;
  for(int k = 0; k < _nQuad; ++k) {

    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, _dxdr);
    form->_geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(form->_geoCoord, k, _dxds);
    jac = form->_J[_nQuad * form->_numElem + k];

    double drdx = _dxds[1] / jac;
    double drdy = -_dxds[0] / jac;
    double dsdx = -_dxdr[1] / jac;
    double dsdy = _dxdr[0] / jac;

    dudx = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdx +
           form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdx;
    dudy = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k) * drdy +
           form->_intSpaces[_idU]->interpolateFieldAtQuadNode_sDerivative(form->_sol[_idU], k) * dsdy;

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    double kD = (*_diffusivity)(form->_tn, _pos);

    for(int i = 0; i < _nFunctions; ++i) {
      _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdx +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdx;
      _feUdy[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k) * drdy +
                  form->_intSpaces[_idU]->getdFunctiondsAtQuadNode(i, k) * dsdy;
      form->_Be[i] -= (_feUdx[i] * dudx + _feUdy[i] * dudy) * kD * jac * _wQuad[k];
    }
  }
}

void feSysElm_2D_Masse::createElementarySystem(std::vector<feSpace *> &space)
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

void feSysElm_2D_Masse::computeAe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNumQuadPoints();
  std::vector<double> &w = form->_geoSpace->getQuadratureWeights();
  double rho = _par;
  int nFunctionsU = form->_intSpaces[_idU]->getNumFunctions();

  double jac, u, dudt;
  for(int k = 0; k < nG; ++k) {
    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, _dxdr);
    form->_geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(form->_geoCoord, k, _dxds);
    jac = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    dudt = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_solDot[_idU], k);

    for(int i = 0; i < nFunctionsU; ++i)
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);

    for(int i = 0; i < nFunctionsU; ++i) {
      for(int j = 0; j < nFunctionsU; ++j) {
        form->_Ae[i][j] += _feU[i] * rho * form->_c0 * _feU[j] * jac * w[k];
      }
    }
  }
}

void feSysElm_2D_Masse::computeBe(feBilinearForm *form)
{
  int nG = form->_geoSpace->getNumQuadPoints();
  std::vector<double> &w = form->_geoSpace->getQuadratureWeights();
  double rho = _par;
  int nFunctionsU = form->_intSpaces[_idU]->getNumFunctions();

  double jac, dudt;
  for(int k = 0; k < nG; ++k) {
    form->_geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(form->_geoCoord, k, _dxdr);
    form->_geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(form->_geoCoord, k, _dxds);
    jac = _dxdr[0] * _dxds[1] - _dxdr[1] * _dxds[0];

    dudt = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_solDot[_idU], k);

    for(int i = 0; i < nFunctionsU; ++i) {
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= _feU[i] * rho * dudt * jac * w[k];
    }
  }
}

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
    _fct->eval(form->_tn, x, v);

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