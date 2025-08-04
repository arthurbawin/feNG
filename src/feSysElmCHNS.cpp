
#include "feBilinearForm.h"
#include "feSysElm.h"

template <int dim>
void CHNS_Abels<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU   = 0;
  _idP   = 1;
  _idPhi = 2;
  _idMu  = 3;

  _fieldsLayoutI = {_idU, _idP, _idPhi, _idMu};
  _fieldsLayoutJ = {_idU, _idP, _idPhi, _idMu};

  _nFunctionsU   = space[_idU]->getNumFunctions();
  _nFunctionsP   = space[_idP]->getNumFunctions();
  _nFunctionsPhi = space[_idPhi]->getNumFunctions();
  _nFunctionsMu  = space[_idMu]->getNumFunctions();

  // Volume force and source term
  _f.resize(dim);
  _Su.resize(dim);

  // Current fields and gradients
  _u.resize(dim);
  _dudt.resize(dim);
  _gradu.resize(dim * dim);
  _symmetricGradu.resize(dim * dim);
  _uDotGradu.resize(dim);
  _gradp.resize(dim);
  _gradphi.resize(dim);
  _gradmu.resize(dim);
  _gradmuDotGradu.resize(dim);

  // Test functions and gradients
  _gradPhiU.resize(dim * dim * _nFunctionsU);
  _phiP.resize(_nFunctionsP);
  _gradPhiP.resize(dim * _nFunctionsP);
  _phiPhi.resize(_nFunctionsPhi);
  _gradPhiPhi.resize(dim * _nFunctionsPhi);
  _phiMu.resize(_nFunctionsMu);
  _gradPhiMu.resize(dim * _nFunctionsMu);

  // Contractions
  _uDotPhiU.resize(_nFunctionsU);
  _dudtDotPhiU.resize(_nFunctionsU);
  _uDotGraduDotPhiU.resize(_nFunctionsU);
  _gradMuDotgradUdotphiU.resize(_nFunctionsU);
  _fDotPhiU.resize(_nFunctionsU);
  _SuDotPhiU.resize(_nFunctionsU);
  _gradPhiDotphiU.resize(_nFunctionsU);
  _gradMuDotphiU.resize(_nFunctionsU);
  _divPhiU.resize(_nFunctionsU);
  _doubleContraction.resize(_nFunctionsU);
}

template <int dim>
void CHNS_Abels<dim>::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Computing with FD");
  exit(-1);
}

template <int dim>
void CHNS_Abels<dim>::computeBe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace =
    static_cast<feVectorSpace<dim> *>(form->_intSpaces[_idU]);
  const feSpace         *pSpace  = form->_intSpaces[_idP];
  const feSpace         *fSpace  = form->_intSpaces[_idPhi];
  const feSpace         *mSpace  = form->_intSpaces[_idMu];
  std::vector<double>   &uSol    = form->_sol[_idU];
  std::vector<double>   &pSol    = form->_sol[_idP];
  std::vector<double>   &fSol    = form->_sol[_idPhi];
  std::vector<double>   &mSol    = form->_sol[_idMu];
  std::vector<double>   &uSolDot = form->_solDot[_idU];
  std::vector<double>   &fSolDot = form->_solDot[_idPhi];
  ElementTransformation &T       = form->_transformation;

  // const double hsize = hCircumInscr(form->_geoCoord);

  double jac, jacw, rho, drhodphi, eta, p, phi, mu, dphidt, div_u;
  double M, Sp, Sphi, Smu;
  double uDotGradPhi;

  for (int k = 0; k < _nQuad; ++k)
  {
    jac  = form->_J[_nQuad * form->_numElem + k];
    jacw = jac * _wQuad[k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, T);

    //
    // Get current fields and their derivatives
    //
    uSpace->interpolateVectorFieldAtQuadNode(uSol, k, _u, dim);
    p   = pSpace->interpolateFieldAtQuadNode(pSol, k);
    phi = fSpace->interpolateFieldAtQuadNode(fSol, k);
    mu  = mSpace->interpolateFieldAtQuadNode(mSol, k);

    //
    // Current time derivatives
    //
    uSpace->interpolateVectorFieldAtQuadNode(uSolDot, k, _dudt, dim);
    dphidt = fSpace->interpolateFieldAtQuadNode(fSolDot, k);

    // Gradients
    uSpace->interpolateVectorFieldAtQuadNode_physicalGradient(
      uSol, dim, k, T, _gradu.data());
    pSpace->interpolateFieldAtQuadNode_physicalGradient(pSol,
                                                        k,
                                                        T,
                                                        _gradp.data());
    fSpace->interpolateFieldAtQuadNode_physicalGradient(fSol,
                                                        k,
                                                        T,
                                                        _gradphi.data());
    mSpace->interpolateFieldAtQuadNode_physicalGradient(mSol,
                                                        k,
                                                        T,
                                                        _gradmu.data());

    // Set phi as current solution value in function arguments
    form->_args.u = phi;

    //
    // Evaluate all scalar coefficients
    //
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord,
                                                      k,
                                                      form->_args.pos);
    rho      = (*_density)(form->_args);
    drhodphi = (*_drhodphi)(form->_args);
    eta      = (*_viscosity)(form->_args);
    M        = (*_mobility)(form->_args);
    (*_volumeForce)(form->_args, _f);

    (*_sourceU)(form->_args, _Su);
    Sp   = (*_sourceP)(form->_args);
    Sphi = (*_sourcePhi)(form->_args);
    Smu  = (*_sourceMu)(form->_args);

    //
    // Test functions and gradients
    // Velocity test functions are accessed through the contraction functions
    //
    pSpace->getFunctionsAtQuadNode(k, _phiP);
    fSpace->getFunctionsAtQuadNode(k, _phiPhi);
    mSpace->getFunctionsAtQuadNode(k, _phiMu);

    uSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiU);
    pSpace->getFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiP.data());
    fSpace->getFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiPhi.data());
    mSpace->getFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiMu.data());

    if constexpr (dim == 2)
    {
      _uDotGradu[0]      = _u[0] * _gradu[0] + _u[1] * _gradu[2];
      _uDotGradu[1]      = _u[0] * _gradu[1] + _u[1] * _gradu[3];
      _gradmuDotGradu[0] = _gradmu[0] * _gradu[0] + _gradmu[1] * _gradu[2];
      _gradmuDotGradu[1] = _gradmu[0] * _gradu[1] + _gradmu[1] * _gradu[3];
      _symmetricGradu[0] = _gradu[0] + _gradu[0];
      _symmetricGradu[1] = _gradu[1] + _gradu[2];
      _symmetricGradu[2] = _gradu[2] + _gradu[1];
      _symmetricGradu[3] = _gradu[3] + _gradu[3];
      div_u              = _gradu[0] + _gradu[3];
      uDotGradPhi        = _u[0] * _gradphi[0] + _u[1] * _gradphi[1];
    }

    uSpace->dotProductShapeOther(k, _u, _uDotPhiU);
    uSpace->dotProductShapeOther(k, _dudt, _dudtDotPhiU);
    uSpace->dotProductShapeOther(k, _uDotGradu, _uDotGraduDotPhiU);
    uSpace->dotProductShapeOther(k, _gradmuDotGradu, _gradMuDotgradUdotphiU);
    uSpace->dotProductShapeOther(k, _f, _fDotPhiU);
    uSpace->dotProductShapeOther(k, _Su, _SuDotPhiU);
    uSpace->dotProductShapeOther(k, _gradphi, _gradPhiDotphiU);
    uSpace->dotProductShapeOther(k, _gradmu, _gradMuDotphiU);
    uSpace->divergence(_gradPhiU, _divPhiU);
    uSpace->doubleContractionGradShapeOther(_gradPhiU,
                                            _symmetricGradu,
                                            _doubleContraction);

    //
    // Increment RHS
    //
    int I = 0;

    //
    // u test functions block : Momentum equation
    //
    for (int i = 0; i < _nFunctionsU; ++i, ++I)
    {
      form->_Be[I] -=
        jacw * (
                 // Acceleration
                 rho * (_dudtDotPhiU[i] + _uDotGraduDotPhiU[i] - _fDotPhiU[i])

                 // Diffusive flux
                 - drhodphi * M * _gradMuDotgradUdotphiU[i]

                 // Korteweg force
                 // - mu * _gradPhiDotphiU[i]
                 + phi * _gradMuDotphiU[i]

                 // Pressure gradient
                 - p * _divPhiU[i]

                 // Viscous term : Symmetric gradient
                 + eta * _doubleContraction[i]

                 // Source term if needed
                 + _SuDotPhiU[i]

               );
    }

    //
    // p test functions block : continuity
    //
    for (int i = 0; i < _nFunctionsP; ++i, ++I)
    {
      form->_Be[I] -= jacw * (div_u * _phiP[i] + Sp * _phiP[i]);
    }

    //
    // phi test functions block : tracer equation
    //
    double gradMu_dot_gradPhiPhi;
    for (int i = 0; i < _nFunctionsPhi; ++i, ++I)
    {
      if constexpr (dim == 2)
      {
        gradMu_dot_gradPhiPhi = _gradmu[0] * _gradPhiPhi[2 * i + 0] +
                                _gradmu[1] * _gradPhiPhi[2 * i + 1];
      }

      form->_Be[I] -= jacw * (
                               // Time derivative
                               dphidt * _phiPhi[i]
                               // Convective term
                               + uDotGradPhi * _phiPhi[i]
                               // Diffusive term
                               + M * gradMu_dot_gradPhiPhi
                               // Source term if needed
                               + Sphi * _phiPhi[i]);
    }

    //
    // mu test functions block : potential equation
    //
    double gradPhi_dot_gradPhiMu;
    for (int i = 0; i < _nFunctionsMu; ++i, ++I)
    {
      if constexpr (dim == 2)
      {
        gradPhi_dot_gradPhiMu = _gradphi[0] * _gradPhiMu[2 * i + 0] +
                                _gradphi[1] * _gradPhiMu[2 * i + 1];
      }

      form->_Be[I] -= jacw * (
                               // Mu term
                               mu * _phiMu[i]
                               // Phi laplacian
                               - _lambda * gradPhi_dot_gradPhiMu
                               // Double well potential
                               - _lambda / (_epsilon * _epsilon) * phi *
                                   (phi * phi - 1.) * _phiMu[i]
                               // Source term if needed
                               + Smu * _phiMu[i]);
    }
  }
}

template class CHNS_Abels<2>;

template <int dim>
void CHNS_MassAveraged<dim>::createElementarySystem(
  std::vector<feSpace *> &space)
{
  _idU   = 0;
  _idP   = 1;
  _idPhi = 2;
  _idMu  = 3;

  _fieldsLayoutI = {_idU, _idP, _idPhi, _idMu};
  _fieldsLayoutJ = {_idU, _idP, _idPhi, _idMu};

  _nFunctionsU   = space[_idU]->getNumFunctions();
  _nFunctionsP   = space[_idP]->getNumFunctions();
  _nFunctionsPhi = space[_idPhi]->getNumFunctions();
  _nFunctionsMu  = space[_idMu]->getNumFunctions();

  _f.resize(dim);
  _Su.resize(dim);

  _u.resize(dim);
  _dudt.resize(dim);
  _gradu.resize(dim * dim);
  _symmetricGradu.resize(dim * dim);
  _uDotGradu.resize(dim);
  _gradp.resize(dim);
  _gradphi.resize(dim);
  _gradmu.resize(dim);
  // _gradmuDotGradu.resize(dim);

  // phiU is the test function of U, unrelated to the phase marker Phi
  // _phiU.resize(_nFunctionsU);
  _gradPhiU.resize(dim * dim * _nFunctionsU);
  _phiP.resize(_nFunctionsP);
  _gradPhiP.resize(dim * _nFunctionsP);
  _phiPhi.resize(_nFunctionsPhi);
  _gradPhiPhi.resize(dim * _nFunctionsPhi);
  _phiMu.resize(_nFunctionsMu);
  _gradPhiMu.resize(dim * _nFunctionsMu);

  _uDotPhiU.resize(_nFunctionsU);
  _dudtDotPhiU.resize(_nFunctionsU);
  _uDotGraduDotPhiU.resize(_nFunctionsU);
  _fDotPhiU.resize(_nFunctionsU);
  _SuDotPhiU.resize(_nFunctionsU);
  _gradMuDotphiU.resize(_nFunctionsU);
  _divPhiU.resize(_nFunctionsU);
  _doubleContraction.resize(_nFunctionsU);

  // _phi_idotphi_j.resize(_nFunctionsU*_nFunctionsU);
  // _u0DotGradPhiUDotPhiU.resize(_nFunctionsU*_nFunctionsU);
  // _phiUDotGradu0DotPhiU.resize(_nFunctionsU*_nFunctionsU);
  // _doubleContractionPhiPhi.resize(_nFunctionsU*_nFunctionsU);
  // _doubleContractionPhiPhiT.resize(_nFunctionsU*_nFunctionsU);
  // _gradMu0DotgradUdotphiU.resize(_nFunctionsU*_nFunctionsU);
  // _gradMuDotgradU0DotphiU.resize(_nFunctionsU*_nFunctionsMu);
  // _gradPhi0dotPhiU.resize(_nFunctionsU);
  // _gradPhiPhiDotPhiU.resize(_nFunctionsU*_nFunctionsPhi);
  // _symGraduDDotGradPhiU.resize(_nFunctionsU);
}

template <int dim>
void CHNS_MassAveraged<dim>::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Computing with FD");
  exit(-1);
}

template <int dim>
void CHNS_MassAveraged<dim>::computeBe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace =
    static_cast<feVectorSpace<dim> *>(form->_intSpaces[_idU]);
  const feSpace         *pSpace  = form->_intSpaces[_idP];
  const feSpace         *fSpace  = form->_intSpaces[_idPhi];
  const feSpace         *mSpace  = form->_intSpaces[_idMu];
  std::vector<double>   &uSol    = form->_sol[_idU];
  std::vector<double>   &pSol    = form->_sol[_idP];
  std::vector<double>   &fSol    = form->_sol[_idPhi];
  std::vector<double>   &fSol_n  = form->_solAtTimeN[_idPhi];
  std::vector<double>   &mSol    = form->_sol[_idMu];
  std::vector<double>   &uSolDot = form->_solDot[_idU];
  std::vector<double>   &fSolDot = form->_solDot[_idPhi];
  ElementTransformation &T       = form->_transformation;

  double jac, jacw, rho, rho_n, drhodphi, eta, p, phi, phi_n, phi_avg, mu, dphidt, div_u, divRhoU,
    uDotGradRho;
  double M, Sp, Sphi, Smu;
  double f_n, f_avg, f_np1, doubleWell, rho_avg;

  for (int k = 0; k < _nQuad; ++k)
  {
    jac  = form->_J[_nQuad * form->_numElem + k];
    jacw = jac * _wQuad[k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, T);

    //
    // Get current fields and their derivatives
    //
    uSpace->interpolateVectorFieldAtQuadNode(uSol, k, _u, dim);
    p   = pSpace->interpolateFieldAtQuadNode(pSol, k);
    phi = fSpace->interpolateFieldAtQuadNode(fSol, k);
    mu  = mSpace->interpolateFieldAtQuadNode(mSol, k);

    // Phi at previous time step
    phi_n = fSpace->interpolateFieldAtQuadNode(fSol_n, k);
    phi_avg = 0.5 * (phi + phi_n);
    
    //
    // Current time derivatives
    //
    uSpace->interpolateVectorFieldAtQuadNode(uSolDot, k, _dudt, dim);
    dphidt = fSpace->interpolateFieldAtQuadNode(fSolDot, k);

    // Gradients
    uSpace->interpolateVectorFieldAtQuadNode_physicalGradient(
      uSol, dim, k, T, _gradu.data());
    pSpace->interpolateFieldAtQuadNode_physicalGradient(pSol,
                                                        k,
                                                        T,
                                                        _gradp.data());
    fSpace->interpolateFieldAtQuadNode_physicalGradient(fSol,
                                                        k,
                                                        T,
                                                        _gradphi.data());
    mSpace->interpolateFieldAtQuadNode_physicalGradient(mSol,
                                                        k,
                                                        T,
                                                        _gradmu.data());

    // Set phi as current solution value in function arguments
    form->_args.u = phi;

    //
    // Evaluate all scalar coefficients
    //
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord,
                                                      k,
                                                      form->_args.pos);
    rho      = (*_density)(form->_args);
    drhodphi = (*_drhodphi)(form->_args);
    eta      = (*_viscosity)(form->_args);
    M        = (*_mobility)(form->_args);
    (*_volumeForce)(form->_args, _f);

    (*_sourceU)(form->_args, _Su);
    Sp   = (*_sourceP)(form->_args);
    Sphi = (*_sourcePhi)(form->_args);
    Smu  = (*_sourceMu)(form->_args);

    // rho at previous time step
    form->_args.u = phi_n;
    rho_n = (*_density)(form->_args);

    // rho at averaged phi
    form->_args.u = phi_avg;
    rho_avg = (*_density)(form->_args);

    //
    // Test functions and gradients
    // Velocity test functions are accessed through the contraction functions
    //
    pSpace->getFunctionsAtQuadNode(k, _phiP);
    fSpace->getFunctionsAtQuadNode(k, _phiPhi);
    mSpace->getFunctionsAtQuadNode(k, _phiMu);

    uSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiU);
    pSpace->getFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiP.data());
    fSpace->getFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiPhi.data());
    mSpace->getFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiMu.data());

    if constexpr (dim == 2)
    {
      _uDotGradu[0]      = _u[0] * _gradu[0] + _u[1] * _gradu[2];
      _uDotGradu[1]      = _u[0] * _gradu[1] + _u[1] * _gradu[3];
      _symmetricGradu[0] = _gradu[0] + _gradu[0];
      _symmetricGradu[1] = _gradu[1] + _gradu[2];
      _symmetricGradu[2] = _gradu[2] + _gradu[1];
      _symmetricGradu[3] = _gradu[3] + _gradu[3];

      div_u       = _gradu[0] + _gradu[3];
      uDotGradRho = drhodphi * (_u[0] * _gradphi[0] + _u[1] * _gradphi[1]);
    }

    divRhoU = rho * div_u + uDotGradRho;

    // Time-averaged double well
    f_n   = phi_n   * (phi_n   * phi_n   - 1.) * _beta;
    f_avg = phi_avg * (phi_avg * phi_avg - 1.) * _beta;
    f_np1 = phi     * (phi     * phi     - 1.) * _beta;
    doubleWell = 1./6. * (f_np1 + 4. * f_avg + f_n);

    uSpace->dotProductShapeOther(k, _u, _uDotPhiU);
    uSpace->dotProductShapeOther(k, _dudt, _dudtDotPhiU);
    uSpace->dotProductShapeOther(k, _uDotGradu, _uDotGraduDotPhiU);
    uSpace->dotProductShapeOther(k, _f, _fDotPhiU);
    uSpace->dotProductShapeOther(k, _Su, _SuDotPhiU);
    uSpace->dotProductShapeOther(k, _gradmu, _gradMuDotphiU);
    uSpace->divergence(_gradPhiU, _divPhiU);
    uSpace->doubleContractionGradShapeOther(_gradPhiU,
                                            _symmetricGradu,
                                            _doubleContraction);

    //
    // Increment RHS
    //
    int I = 0;
    UNUSED(rho_n, rho_avg, doubleWell);
    //
    // u test functions block : Momentum equation
    //
    for (int i = 0; i < _nFunctionsU; ++i, ++I)
    {
      form->_Be[I] -=
        jacw *
        (
          // Acceleration
          // rho_n * _dudtDotPhiU[i] + rho * (_uDotGraduDotPhiU[i])  - rho_avg * _fDotPhiU[i]
          rho * (_dudtDotPhiU[i] + _uDotGraduDotPhiU[i] - _fDotPhiU[i])

          // Mass conservation
          + 0.5 * _uDotPhiU[i] * (drhodphi * dphidt + divRhoU)

          // Korteweg force
          + phi * _gradMuDotphiU[i]

          // Pressure gradient
          - p * _divPhiU[i]

          // Viscous term : Symmetric gradient and compressible contribution
          + eta *
              (_doubleContraction[i] - (2. / (double)dim) * div_u * _divPhiU[i])

          // Source term if needed
          + _SuDotPhiU[i]

        );
    }

    //
    // p test functions block : continuity
    //
    double gradMu_dot_gradPhiP, gradP_dot_gradPhiP;
    for (int i = 0; i < _nFunctionsP; ++i, ++I)
    {
      if constexpr (dim == 2)
      {
        gradMu_dot_gradPhiP =
          _gradmu[0] * _gradPhiP[2 * i + 0] + _gradmu[1] * _gradPhiP[2 * i + 1];
        gradP_dot_gradPhiP =
          _gradp[0] * _gradPhiP[2 * i + 0] +  _gradp[1] * _gradPhiP[2 * i + 1];
      }

      form->_Be[I] -=
        jacw *
        (
          // Div u
          div_u * _phiP[i]
          
          // alpha * Div (M (grad(mu + alpha * p))
          + _alpha * M * (gradMu_dot_gradPhiP + _alpha * gradP_dot_gradPhiP)

          // Source term if needed
          + Sp * _phiP[i]);
    }

    //
    // phi test functions block : tracer equation
    //
    double gradMu_dot_gradPhiPhi, gradP_dot_gradPhiPhi, uDotGradPhiPhi;
    for (int i = 0; i < _nFunctionsPhi; ++i, ++I)
    {
      if constexpr (dim == 2)
      {
        uDotGradPhiPhi =
          _u[0] * _gradPhiPhi[2 * i + 0] + _u[1] * _gradPhiPhi[2 * i + 1];
        gradMu_dot_gradPhiPhi = _gradmu[0] * _gradPhiPhi[2 * i + 0] +
                                _gradmu[1] * _gradPhiPhi[2 * i + 1];
        gradP_dot_gradPhiPhi = _gradp[0] * _gradPhiPhi[2 * i + 0] +
                               _gradp[1] * _gradPhiPhi[2 * i + 1];
      }

      form->_Be[I] -=
        jacw * (
                 // Time derivative
                 dphidt * _phiPhi[i]

                 // Convective term
                 - phi * uDotGradPhiPhi

                 // Diffusive term
                 + M * (gradMu_dot_gradPhiPhi + _alpha * gradP_dot_gradPhiPhi)

                 // Source term if needed
                 + Sphi * _phiPhi[i]);
    }

    //
    // mu test functions block : potential equation
    //
    double gradPhi_dot_gradPhiMu;
    for (int i = 0; i < _nFunctionsMu; ++i, ++I)
    {
      if constexpr (dim == 2)
      {
        gradPhi_dot_gradPhiMu = _gradphi[0] * _gradPhiMu[2 * i + 0] +
                                _gradphi[1] * _gradPhiMu[2 * i + 1];
      }

      form->_Be[I] -= jacw * (
                               // Mu term
                               mu * _phiMu[i]

                               // Phi laplacian
                               - _tau * gradPhi_dot_gradPhiMu

                               // Double well potential
                               // - _beta * (phi * phi * phi - phi) * _phiMu[i]
                               - doubleWell * _phiMu[i]

                               // Source term if needed
                               + Smu * _phiMu[i]);
    }
  }
}

template class CHNS_MassAveraged<2>;

template <int dim>
void CHNS_Khanwale<dim>::createElementarySystem(
  std::vector<feSpace *> &space)
{
  _idU   = 0;
  _idP   = 1;
  _idPhi = 2;
  _idMu  = 3;

  _fieldsLayoutI = {_idU, _idP, _idPhi, _idMu};
  _fieldsLayoutJ = {_idU, _idP, _idPhi, _idMu};

  _nFunctionsU   = space[_idU]->getNumFunctions();
  _nFunctionsP   = space[_idP]->getNumFunctions();
  _nFunctionsPhi = space[_idPhi]->getNumFunctions();
  _nFunctionsMu  = space[_idMu]->getNumFunctions();

  _f = {0., -1., 0.};
  _Su.resize(dim);

  _u.resize(dim);
  _u_n.resize(dim);
  _u_avg.resize(dim);

  _gradu.resize(dim * dim);
  _gradp.resize(dim);
  _gradphi.resize(dim);
  _gradmu.resize(dim);

  _gradu_n.resize(dim * dim);
  _gradphi_n.resize(dim);
  _gradmu_n.resize(dim);
  _gradu_avg.resize(dim * dim);
  _gradphi_avg.resize(dim);
  _gradmu_avg.resize(dim);

  _jflux.resize(dim);
  _dudt.resize(dim);
  _symmetricGradu.resize(dim * dim);
  _uDotGradu.resize(dim);
  _jFluxDotGradu.resize(dim);

  _gphiOtimesgphi.resize(dim * dim);

  _gradPhiU.resize(dim * dim * _nFunctionsU);
  _phiP.resize(_nFunctionsP);
  _gradPhiP.resize(dim * _nFunctionsP);
  _phiPhi.resize(_nFunctionsPhi);
  _gradPhiPhi.resize(dim * _nFunctionsPhi);
  _phiMu.resize(_nFunctionsMu);
  _gradPhiMu.resize(dim * _nFunctionsMu);

  // _uDotPhiU.resize(_nFunctionsU);
  _dudtDotPhiU.resize(_nFunctionsU);
  _uDotGraduDotPhiU.resize(_nFunctionsU);
  _jFluxDotGraduDotPhiU.resize(_nFunctionsU);
  _fDotPhiU.resize(_nFunctionsU);
  _SuDotPhiU.resize(_nFunctionsU);
  _divPhiU.resize(_nFunctionsU);
  _doubleContraction.resize(_nFunctionsU);
  _gphiOtimesgphi_times_gradPhiU.resize(_nFunctionsU);
}

template <int dim>
void CHNS_Khanwale<dim>::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Computing with FD");
  exit(-1);
}

template <int dim>
void CHNS_Khanwale<dim>::computeBe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace =
    static_cast<feVectorSpace<dim> *>(form->_intSpaces[_idU]);
  const feSpace         *pSpace  = form->_intSpaces[_idP];
  const feSpace         *fSpace  = form->_intSpaces[_idPhi];
  const feSpace         *mSpace  = form->_intSpaces[_idMu];

  std::vector<double>   &uSol    = form->_sol[_idU];
  std::vector<double>   &pSol    = form->_sol[_idP];
  std::vector<double>   &fSol    = form->_sol[_idPhi];
  std::vector<double>   &mSol    = form->_sol[_idMu];

  std::vector<double>   &uSol_n  = form->_solAtTimeN[_idU];
  std::vector<double>   &pSol_n  = form->_solAtTimeN[_idP];
  std::vector<double>   &fSol_n  = form->_solAtTimeN[_idPhi];
  std::vector<double>   &mSol_n  = form->_solAtTimeN[_idMu];

  std::vector<double>   &uSolDot = form->_solDot[_idU];
  std::vector<double>   &fSolDot = form->_solDot[_idPhi];

  ElementTransformation &T       = form->_transformation;

  for (int k = 0; k < _nQuad; ++k)
  {
    const double jac  = form->_J[_nQuad * form->_numElem + k];
    const double jacw = jac * _wQuad[k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, T);
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord,
                                                      k,
                                                      form->_args.pos);

    //
    // Current fields
    //
    uSpace->interpolateVectorFieldAtQuadNode(uSol, k, _u, dim);
    const double p   = pSpace->interpolateFieldAtQuadNode(pSol, k);
    const double phi = fSpace->interpolateFieldAtQuadNode(fSol, k);
    const double mu  = mSpace->interpolateFieldAtQuadNode(mSol, k);

    // Time averages
    uSpace->interpolateVectorFieldAtQuadNode(uSol_n, k, _u_n, dim);
    const double p_n   = pSpace->interpolateFieldAtQuadNode(pSol_n, k);
    const double phi_n = fSpace->interpolateFieldAtQuadNode(fSol_n, k);
    const double mu_n  = mSpace->interpolateFieldAtQuadNode(mSol_n, k);

    for(unsigned int d = 0; d < dim; ++d)
    {
      _u_avg[d] = 0.5 * (_u[d] + _u_n[d]);
    }
    const double p_avg   = 0.5 * (p + p_n);
    const double phi_avg = 0.5 * (phi + phi_n);
    const double mu_avg  = 0.5 * (mu + mu_n);
    
    //
    // Current time derivatives (Implicit Euler for this formulation)
    //
    uSpace->interpolateVectorFieldAtQuadNode(uSolDot, k, _dudt, dim);
    const double dphidt = fSpace->interpolateFieldAtQuadNode(fSolDot, k);

    // Current and averaged gradients
    uSpace->interpolateVectorFieldAtQuadNode_physicalGradient(uSol, dim, k, T, _gradu.data());
    pSpace->interpolateFieldAtQuadNode_physicalGradient(pSol, k, T, _gradp.data());
    fSpace->interpolateFieldAtQuadNode_physicalGradient(fSol, k, T, _gradphi.data());
    mSpace->interpolateFieldAtQuadNode_physicalGradient(mSol, k, T, _gradmu.data());

    uSpace->interpolateVectorFieldAtQuadNode_physicalGradient(uSol_n, dim, k, T, _gradu_n.data());
    fSpace->interpolateFieldAtQuadNode_physicalGradient(fSol_n, k, T, _gradphi_n.data());
    mSpace->interpolateFieldAtQuadNode_physicalGradient(mSol_n, k, T, _gradmu_n.data());

    for(unsigned int d = 0; d < dim; ++d)
    {
      for(unsigned int e = 0; e < dim; ++e)
      {
        _gradu_avg[d * dim + e] = 0.5 * (_gradu[d * dim + e] + _gradu_n[d * dim + e]);
      }
      _gradphi_avg[d] = 0.5 * (_gradphi[d] + _gradphi_n[d]);
      _gradmu_avg[d]  = 0.5 * (_gradmu[d] + _gradmu_n[d]);
    }

    // Evaluate source terms at previous time step (n)
    form->_args.u       = phi_n;
    const double rho_n  = (*_density)(form->_args);

    // Evaluate coefficients at average time
    form->_args.u         = phi_avg;
    const double rho_avg  = (*_density)(form->_args);
    const double eta_avg  = (*_viscosity)(form->_args);
    const double well_avg = phi_avg * (phi_avg * phi_avg - 1.);

    // Evaluate source terms at next time step (n+1)
    form->_args.u          = phi;
    const double rho       = (*_density)(form->_args);
    const double drhodphi  = (*_drhodphi)(form->_args);
    (*_sourceU)(form->_args, _Su);
    const double Sp   = (*_sourceP)(form->_args);
    const double Sphi = (*_sourcePhi)(form->_args);
    const double Smu  = (*_sourceMu)(form->_args);

    double div_u = 0., div_u_avg = 0.;
    if constexpr (dim == 2)
    {
      div_u              = _gradu[0] + _gradu[3];
      div_u_avg          = _gradu_avg[0] + _gradu_avg[3];

      _uDotGradu[0]      = _u_avg[0] * _gradu_avg[0] + _u_avg[1] * _gradu_avg[2];
      _uDotGradu[1]      = _u_avg[0] * _gradu_avg[1] + _u_avg[1] * _gradu_avg[3];
      _symmetricGradu[0] = _gradu_avg[0] + _gradu_avg[0];
      _symmetricGradu[1] = _gradu_avg[1] + _gradu_avg[2];
      _symmetricGradu[2] = _gradu_avg[2] + _gradu_avg[1];
      _symmetricGradu[3] = _gradu_avg[3] + _gradu_avg[3];
      _jflux[0]          = (_rhoB - _rhoA) / (2. * _rhoA * _Cn) * _gradmu_avg[0];
      _jflux[1]          = (_rhoB - _rhoA) / (2. * _rhoA * _Cn) * _gradmu_avg[1];
      _jFluxDotGradu[0]  = _jflux[0] * _gradu_avg[0] + _jflux[1] * _gradu_avg[2];
      _jFluxDotGradu[1]  = _jflux[0] * _gradu_avg[1] + _jflux[1] * _gradu_avg[3];
      _gphiOtimesgphi[0] = _gradphi_avg[0] * _gradphi_avg[0];
      _gphiOtimesgphi[1] = _gradphi_avg[0] * _gradphi_avg[1];
      _gphiOtimesgphi[2] = _gradphi_avg[1] * _gradphi_avg[0];
      _gphiOtimesgphi[3] = _gradphi_avg[1] * _gradphi_avg[1];
    }

    const double gRhoAvgDotUavg = drhodphi * (_gradphi_avg[0] * _u_avg[0] + _gradphi_avg[1] * _u_avg[1]);
    const double div_rhoAvgUavg = gRhoAvgDotUavg + rho_avg * div_u_avg;
    //
    // Test functions and gradients
    // Velocity test functions are accessed through the contraction functions
    //
    pSpace->getFunctionsAtQuadNode(k, _phiP);
    fSpace->getFunctionsAtQuadNode(k, _phiPhi);
    mSpace->getFunctionsAtQuadNode(k, _phiMu);

    uSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiU);
    pSpace->getFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiP.data());
    fSpace->getFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiPhi.data());
    mSpace->getFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiMu.data());

    uSpace->dotProductShapeOther(k, _dudt, _dudtDotPhiU);
    uSpace->dotProductShapeOther(k, _uDotGradu, _uDotGraduDotPhiU);
    uSpace->dotProductShapeOther(k, _jFluxDotGradu, _jFluxDotGraduDotPhiU);
    uSpace->dotProductShapeOther(k, _f, _fDotPhiU);
    uSpace->dotProductShapeOther(k, _Su, _SuDotPhiU);
    uSpace->divergence(_gradPhiU, _divPhiU);
    uSpace->doubleContractionGradShapeOther(_gradPhiU,
                                            _symmetricGradu,
                                            _doubleContraction);
    uSpace->doubleContractionGradShapeOther(_gradPhiU, _gphiOtimesgphi, _gphiOtimesgphi_times_gradPhiU);

    //
    // Increment RHS
    //
    int I = 0;

    //
    // u test functions block : Momentum equation
    //
    for (int i = 0; i < _nFunctionsU; ++i, ++I)
    {
      form->_Be[I] -=
        jacw *
        (
          rho_avg * (_dudtDotPhiU[i] + _uDotGraduDotPhiU[i])

          + 1./_Pe * _jFluxDotGraduDotPhiU[i]

          - _Cn/_We * _gphiOtimesgphi_times_gradPhiU[i]

          - 1./_We * p_avg * _divPhiU[i]

          + 1./_Re * eta_avg * _doubleContraction[i]

          - 1./_Fr * rho_avg * _fDotPhiU[i]

          // Source term if needed
          + _SuDotPhiU[i]
        );
    }

    //
    // p test functions block : continuity
    //
    double jFlux_dot_gradPhiP;
    for (int i = 0; i < _nFunctionsP; ++i, ++I)
    {
      if constexpr (dim == 2)
      {
        jFlux_dot_gradPhiP =
          _jflux[0] * _gradPhiP[2 * i + 0] + _jflux[1] * _gradPhiP[2 * i + 1];
      }

      form->_Be[I] -=
        jacw *
        (
          // Div u
          div_u * _phiP[i]
          
          + (rho - rho_n) / form->_dt * _phiP[i]

          + div_rhoAvgUavg * _phiP[i]

          - 1./_Pe * jFlux_dot_gradPhiP

          // Source term if needed
          + Sp * _phiP[i]);
    }

    //
    // phi test functions block : tracer equation
    //
    double gradMuAvg_dot_gradPhiPhi, uAvgDotGradPhiPhi;
    for (int i = 0; i < _nFunctionsPhi; ++i, ++I)
    {
      if constexpr (dim == 2)
      {
        uAvgDotGradPhiPhi =
          _u_avg[0] * _gradPhiPhi[2 * i + 0] + _u_avg[1] * _gradPhiPhi[2 * i + 1];
        gradMuAvg_dot_gradPhiPhi = _gradmu_avg[0] * _gradPhiPhi[2 * i + 0] +
                                   _gradmu_avg[1] * _gradPhiPhi[2 * i + 1];
      }

      form->_Be[I] -=
        jacw * (
                 // Time derivative
                 dphidt * _phiPhi[i]

                 // Convective term
                 - phi_avg * uAvgDotGradPhiPhi

                 // Diffusive term
                 + 1./(_Pe*_Cn) * gradMuAvg_dot_gradPhiPhi

                 // Source term if needed
                 + Sphi * _phiPhi[i]);
    }

    //
    // mu test functions block : potential equation
    //
    double gradPhiAvg_dot_gradPhiMu;
    for (int i = 0; i < _nFunctionsMu; ++i, ++I)
    {
      if constexpr (dim == 2)
      {
        gradPhiAvg_dot_gradPhiMu = _gradphi_avg[0] * _gradPhiMu[2 * i + 0] +
                                   _gradphi_avg[1] * _gradPhiMu[2 * i + 1];
      }

      form->_Be[I] -= jacw * (
                               // Mu term
                               mu_avg * _phiMu[i]

                               // Double well potential
                               - well_avg * _phiMu[i]

                               // Phi laplacian
                               - _Cn * _Cn * gradPhiAvg_dot_gradPhiMu

                               // Source term if needed
                               + Smu * _phiMu[i]);
    }
  }
}

template class CHNS_Khanwale<2>;

template <int dim>
void CHNS_VolumeAveragedGeneric<dim>::createElementarySystem(
  std::vector<feSpace *> &space)
{
  _idU   = 0;
  _idP   = 1;
  _idPhi = 2;
  _idMu  = 3;
  _idV   = 4;

  _fieldsLayoutI = {_idU, _idP, _idPhi, _idMu, _idV};
  _fieldsLayoutJ = {_idU, _idP, _idPhi, _idMu, _idV};

  _nFunctionsU   = space[_idU]->getNumFunctions();
  _nFunctionsP   = space[_idP]->getNumFunctions();
  _nFunctionsPhi = space[_idPhi]->getNumFunctions();
  _nFunctionsMu  = space[_idMu]->getNumFunctions();
  _nFunctionsV   = space[_idV]->getNumFunctions();

  _f.resize(dim);
  _Su.resize(dim);
  _Sv.resize(dim);

  _u.resize(dim);
  _v.resize(dim);
  _dudt.resize(dim);
  _dvdt.resize(dim);
  _gradu.resize(dim * dim);
  _gradv.resize(dim * dim);
  _symmetricGradu.resize(dim * dim);
  _uDotGradu.resize(dim);
  _gradp.resize(dim);
  _gradphi.resize(dim);
  _gradmu.resize(dim);
  // _gradmuDotGradu.resize(dim);

  // phiU is the test function of U, unrelated to the phase marker Phi
  // _phiU.resize(_nFunctionsU);
  _gradPhiU.resize(dim * dim * _nFunctionsU);
  _gradPhiV.resize(dim * dim * _nFunctionsV);
  _phiP.resize(_nFunctionsP);
  _gradPhiP.resize(dim * _nFunctionsP);
  _phiPhi.resize(_nFunctionsPhi);
  _gradPhiPhi.resize(dim * _nFunctionsPhi);
  _phiMu.resize(_nFunctionsMu);
  _gradPhiMu.resize(dim * _nFunctionsMu);

  _uDotPhiU.resize(_nFunctionsU);
  _vDotPhiU.resize(_nFunctionsU);
  _dudtDotPhiU.resize(_nFunctionsU);
  _uDotGraduDotPhiU.resize(_nFunctionsU);
  _fDotPhiU.resize(_nFunctionsU);
  _SuDotPhiU.resize(_nFunctionsU);
  _gradMuDotphiU.resize(_nFunctionsU);
  _divPhiU.resize(_nFunctionsU);
  _doubleContraction.resize(_nFunctionsU);

  _JfluxDotPhiU.resize(_nFunctionsU);

  _vDotPhiV.resize(_nFunctionsV);
  _dvdtDotPhiV.resize(_nFunctionsV);
  _SvDotPhiV.resize(_nFunctionsV);
  _fDotPhiV.resize(_nFunctionsV);
  // _uDotGraduDotPhiU.resize(_nFunctionsU);
  // _fDotPhiU.resize(_nFunctionsU);
  // _SuDotPhiU.resize(_nFunctionsU);
  // _gradMuDotphiU.resize(_nFunctionsU);
  // _divPhiU.resize(_nFunctionsU);
  // _doubleContraction.resize(_nFunctionsU);

  _Jflux.resize(dim);
  _viscousStress.resize(dim * dim);
  _SDotGradPhiV.resize(_nFunctionsV);
  _momentumTensor.resize(dim * dim);
  _TDotGradPhiV.resize(_nFunctionsV);
}

template <int dim>
void CHNS_VolumeAveragedGeneric<dim>::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Computing with FD");
  exit(-1);
}

template <int dim>
void CHNS_VolumeAveragedGeneric<dim>::computeBe(feBilinearForm *form)
{
  const feVectorSpace<dim> *uSpace =
    static_cast<feVectorSpace<dim> *>(form->_intSpaces[_idU]);
  const feVectorSpace<dim> *vSpace =
    static_cast<feVectorSpace<dim> *>(form->_intSpaces[_idV]);
  const feSpace       *pSpace = form->_intSpaces[_idP];
  const feSpace       *fSpace = form->_intSpaces[_idPhi];
  const feSpace       *mSpace = form->_intSpaces[_idMu];
  std::vector<double> &uSol   = form->_sol[_idU];
  std::vector<double> &vSol   = form->_sol[_idV];
  std::vector<double> &pSol   = form->_sol[_idP];
  std::vector<double> &fSol   = form->_sol[_idPhi];
  std::vector<double> &mSol   = form->_sol[_idMu];
  // std::vector<double> &uSolDot = form->_solDot[_idU];
  std::vector<double> &vSolDot = form->_solDot[_idV];
  std::vector<double> &fSolDot = form->_solDot[_idPhi];

  ElementTransformation &T = form->_transformation;

  double jac, jacw, rho, drhodphi, eta, p, phi, mu, dphidt, div_u, div_v, drhodt;
  double M, Sp, Sphi, Smu;

  double bulk_visc = -2. / (double) dim;

  for (int k = 0; k < _nQuad; ++k)
  {
    jac  = form->_J[_nQuad * form->_numElem + k];
    jacw = jac * _wQuad[k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, T);

    //
    // Get current fields and their derivatives
    //
    uSpace->interpolateVectorFieldAtQuadNode(uSol, k, _u, dim);
    vSpace->interpolateVectorFieldAtQuadNode(vSol, k, _v, dim);
    p   = pSpace->interpolateFieldAtQuadNode(pSol, k);
    phi = fSpace->interpolateFieldAtQuadNode(fSol, k);
    mu  = mSpace->interpolateFieldAtQuadNode(mSol, k);

    //
    // Current time derivatives
    //
    // uSpace->interpolateVectorFieldAtQuadNode(uSolDot, k, _dudt, dim);
    vSpace->interpolateVectorFieldAtQuadNode(vSolDot, k, _dvdt, dim);
    dphidt = fSpace->interpolateFieldAtQuadNode(fSolDot, k);

    // Gradients
    uSpace->interpolateVectorFieldAtQuadNode_physicalGradient(
      uSol, dim, k, T, _gradu.data());
    vSpace->interpolateVectorFieldAtQuadNode_physicalGradient(
      vSol, dim, k, T, _gradv.data());
    pSpace->interpolateFieldAtQuadNode_physicalGradient(pSol,
                                                        k,
                                                        T,
                                                        _gradp.data());
    fSpace->interpolateFieldAtQuadNode_physicalGradient(fSol,
                                                        k,
                                                        T,
                                                        _gradphi.data());
    mSpace->interpolateFieldAtQuadNode_physicalGradient(mSol,
                                                        k,
                                                        T,
                                                        _gradmu.data());

    // Set phi as current solution value in function arguments
    form->_args.u = phi;

    //
    // Evaluate all scalar coefficients
    //
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord,
                                                      k,
                                                      form->_args.pos);
    rho      = (*_density)(form->_args);
    drhodphi = (*_drhodphi)(form->_args);
    eta      = (*_viscosity)(form->_args);
    M        = (*_mobility)(form->_args);
    (*_volumeForce)(form->_args, _f);

    (*_sourceU)(form->_args, _Su);
    if(_sourceV) {
      (*_sourceV)(form->_args, _Sv);
    }
    Sp   = (*_sourceP)(form->_args);
    Sphi = (*_sourcePhi)(form->_args);
    Smu  = (*_sourceMu)(form->_args);

    //
    // Test functions and gradients
    // Velocity test functions are accessed through the contraction functions
    //
    pSpace->getFunctionsAtQuadNode(k, _phiP);
    fSpace->getFunctionsAtQuadNode(k, _phiPhi);
    mSpace->getFunctionsAtQuadNode(k, _phiMu);

    uSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiU);
    vSpace->getVectorFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiV);
    pSpace->getFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiP.data());
    fSpace->getFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiPhi.data());
    mSpace->getFunctionsPhysicalGradientAtQuadNode(k, T, _gradPhiMu.data());

    if constexpr (dim == 2)
    {
      // _uDotGradu[0]      = _u[0] * _gradu[0] + _u[1] * _gradu[2];
      // _uDotGradu[1]      = _u[0] * _gradu[1] + _u[1] * _gradu[3];

      _Jflux[0]          = - drhodphi * M * (_gradmu[0] + _alpha * _gradp[0]);
      _Jflux[1]          = - drhodphi * M * (_gradmu[1] + _alpha * _gradp[1]);

      // entries 00, 01, 10, 11
      _momentumTensor[0] = rho * _u[0]*_u[0] + _u[0]*_Jflux[0] + _Jflux[0]*_u[0] + _Jflux[0] * _Jflux[0] / rho;
      _momentumTensor[1] = rho * _u[0]*_u[1] + _u[0]*_Jflux[1] + _Jflux[0]*_u[1] + _Jflux[0] * _Jflux[1] / rho;
      _momentumTensor[2] = rho * _u[1]*_u[0] + _u[1]*_Jflux[0] + _Jflux[1]*_u[0] + _Jflux[1] * _Jflux[0] / rho;
      _momentumTensor[3] = rho * _u[1]*_u[1] + _u[1]*_Jflux[1] + _Jflux[1]*_u[1] + _Jflux[1] * _Jflux[1] / rho;

      // _symmetricGradu[0] = _gradu[0] + _gradu[0];
      // _symmetricGradu[1] = _gradu[1] + _gradu[2];
      // _symmetricGradu[2] = _gradu[2] + _gradu[1];
      // _symmetricGradu[3] = _gradu[3] + _gradu[3];

      div_u       = _gradu[0] + _gradu[3];
      div_v       = _gradv[0] + _gradv[3];

      _viscousStress[0] = eta * (_gradv[0] + _gradv[0] + bulk_visc * div_v);
      _viscousStress[1] = eta * (_gradv[1] + _gradv[2]);
      _viscousStress[2] = eta * (_gradv[2] + _gradv[1]);
      _viscousStress[3] = eta * (_gradv[3] + _gradv[3] + bulk_visc * div_v);

      // uDotGradRho = drhodphi * (_u[0] * _gradphi[0] + _u[1] * _gradphi[1]);
    }

    drhodt  = drhodphi * dphidt;

    uSpace->dotProductShapeOther(k, _u, _uDotPhiU);
    uSpace->dotProductShapeOther(k, _v, _vDotPhiU);
    uSpace->dotProductShapeOther(k, _Jflux, _JfluxDotPhiU);
    // uSpace->dotProductShapeOther(k, _dudt, _dudtDotPhiU);
    // uSpace->dotProductShapeOther(k, _uDotGradu, _uDotGraduDotPhiU);
    uSpace->dotProductShapeOther(k, _Su, _SuDotPhiU);
    // uSpace->dotProductShapeOther(k, _gradmu, _gradMuDotphiU);
    // uSpace->divergence(_gradPhiU, _divPhiU);
    // uSpace->doubleContractionGradShapeOther(_gradPhiU,
    //                                         _symmetricGradu,
    //                                         _doubleContraction);

    vSpace->dotProductShapeOther(k, _v, _vDotPhiV);
    vSpace->dotProductShapeOther(k, _dvdt, _dvdtDotPhiV);
    vSpace->dotProductShapeOther(k, _Sv, _SvDotPhiV);
    vSpace->dotProductShapeOther(k, _f, _fDotPhiV);
    vSpace->doubleContractionGradShapeOther(_gradPhiV, _momentumTensor, _TDotGradPhiV);
    vSpace->doubleContractionGradShapeOther(_gradPhiV, _viscousStress,  _SDotGradPhiV);

    //
    // Increment RHS
    //
    int I = 0;

    //
    // u test functions block : u-v coupling
    //
    for (int i = 0; i < _nFunctionsU; ++i, ++I)
    {
      form->_Be[I] -=
        jacw *
        (
          // Coupling : rho*v = rho*u + J
          rho * _vDotPhiU[i] - rho * _uDotPhiU[i] - _JfluxDotPhiU[i]

          // Source term if needed
          + _SuDotPhiU[i]
        );
    }

    //
    // p test functions block : continuity
    //
    double gradMu_dot_gradPhiP, gradP_dot_gradPhiP;
    for (int i = 0; i < _nFunctionsP; ++i, ++I)
    {
      if constexpr (dim == 2)
      {
        gradMu_dot_gradPhiP =
          _gradmu[0] * _gradPhiP[2 * i + 0] + _gradmu[1] * _gradPhiP[2 * i + 1];
        gradP_dot_gradPhiP =
          _gradp[0] * _gradPhiP[2 * i + 0] + _gradp[1] * _gradPhiP[2 * i + 1];
      }

      form->_Be[I] -=
        jacw *
        (
          // Div u
          div_u * _phiP[i]
          // Source term if needed
          + Sp * _phiP[i]);
    }

    //
    // phi test functions block : tracer equation
    //
    double gradMu_dot_gradPhiPhi, gradP_dot_gradPhiPhi, uDotGradPhiPhi;
    for (int i = 0; i < _nFunctionsPhi; ++i, ++I)
    {
      if constexpr (dim == 2)
      {
        uDotGradPhiPhi =
          _u[0] * _gradPhiPhi[2 * i + 0] + _u[1] * _gradPhiPhi[2 * i + 1];
        gradMu_dot_gradPhiPhi = _gradmu[0] * _gradPhiPhi[2 * i + 0] +
                                _gradmu[1] * _gradPhiPhi[2 * i + 1];
        gradP_dot_gradPhiPhi = _gradp[0] * _gradPhiPhi[2 * i + 0] +
                               _gradp[1] * _gradPhiPhi[2 * i + 1];
      }

      form->_Be[I] -=
        jacw * (
                 // Time derivative
                 dphidt * _phiPhi[i]
                 // Convective term
                 - phi * uDotGradPhiPhi
                 // Diffusive term
                 + M * (gradMu_dot_gradPhiPhi + _alpha * gradP_dot_gradPhiPhi)
                 // Source term if needed
                 + Sphi * _phiPhi[i]);
    }

    //
    // mu test functions block : potential equation
    //
    double gradPhi_dot_gradPhiMu;
    for (int i = 0; i < _nFunctionsMu; ++i, ++I)
    {
      if constexpr (dim == 2)
      {
        gradPhi_dot_gradPhiMu = _gradphi[0] * _gradPhiMu[2 * i + 0] +
                                _gradphi[1] * _gradPhiMu[2 * i + 1];
      }

      form->_Be[I] -= jacw * (
                               // Mu term
                               mu * _phiMu[i]
                               // Phi laplacian
                               - _lambda * gradPhi_dot_gradPhiMu
                               // Double well potential
                               - _lambda / (_epsilon * _epsilon) * phi *
                                   (phi * phi - 1.) * _phiMu[i]
                               // Source term if needed
                               + Smu * _phiMu[i]);
    }

    //
    // v test functions block : Momentum equation
    //
    for (int i = 0; i < _nFunctionsV; ++i, ++I)
    {
      form->_Be[I] -=
        jacw *
        (
          // d/dt (rho*v)
          drhodt * _vDotPhiV[i] + rho * _dvdtDotPhiV[i]

          // Acceleration
          + _TDotGradPhiV[i]

          // Gravity
          - rho * _fDotPhiV[i]

          // Pressure gradient
          - p * _divPhiU[i]

          // Viscous term : Symmetric gradient and compressible contribution
          + _SDotGradPhiV[i]

          // Korteweg force
          + phi * _gradMuDotphiU[i]

          // Source term if needed
          + _SvDotPhiV[i]

        );
    }
  }
}

template class CHNS_VolumeAveragedGeneric<2>;