#include "feSysElm.h"
#include "feBilinearForm.h"
#include "feMatrixInterface.h"

// -----------------------------------------------------------------------------
// Linear form: scalar source term
// -----------------------------------------------------------------------------
void feSysElm_Source::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _phiU.resize(_nFunctions);
}

void feSysElm_Source::computeBe(feBilinearForm *form)
{
  for(int k = 0; k < _nQuad; ++k) {
    double jac = form->_J[_nQuad * form->_numElem + k];
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    double S = (*_source)(form->_tn, _pos);
    for(int i = 0; i < _nFunctions; ++i) {
      _phiU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= _phiU[i] * S * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: mass matrix for scalar field
// -----------------------------------------------------------------------------
void feSysElm_Mass::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _phiU.resize(_nFunctions);
}

void feSysElm_Mass::computeAe(feBilinearForm *form)
{
  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    for(int i = 0; i < _nFunctions; ++i)
      _phiU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);

    for(int i = 0; i < _nFunctions; ++i)
      for(int j = 0; j < _nFunctions; ++j)
        form->_Ae[i][j] += coeff * _phiU[i] * _phiU[j] * jac * _wQuad[k];
  }
}

void feSysElm_Mass::computeBe(feBilinearForm *form)
{
  double jac, u, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);

    for(int i = 0; i < _nFunctions; ++i) {
      _phiU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= coeff * u * _phiU[i] * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: transient mass matrix
// -----------------------------------------------------------------------------
void feSysElm_TransientMass::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _feU.resize(_nFunctions);
}

void feSysElm_TransientMass::computeAe(feBilinearForm *form)
{
  double jac;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    for(int i = 0; i < _nFunctions; ++i)
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);

    for(int i = 0; i < _nFunctions; ++i)
      for(int j = 0; j < _nFunctions; ++j)
        form->_Ae[i][j] += _feU[i] * _coeff * form->_c0 * _feU[j] * jac * _wQuad[k];
  }
}

void feSysElm_TransientMass::computeBe(feBilinearForm *form)
{
  double jac, uDot;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    uDot = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_solDot[_idU], k);

    for(int i = 0; i < _nFunctions; ++i) {
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= _feU[i] * _coeff * uDot * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: diffusion (stiffness matrix) of scalar field
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_Diffusion<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _gradPhi.resize(dim * _nFunctions);
}

template <int dim> void feSysElm_Diffusion<dim>::computeAe(feBilinearForm *form)
{
  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar parameter
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    // Compute grad(phi)
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhi.data());

    for(int i = 0; i < _nFunctions; ++i) {
      for(int j = 0; j < _nFunctions; ++j) {
        for(int iDim = 0; iDim < dim; ++iDim) {
          form->_Ae[i][j] -=
            coeff * _gradPhi[i * dim + iDim] * _gradPhi[j * dim + iDim] * jac * _wQuad[k];
        }
      }
    }
  }
}

template <int dim> void feSysElm_Diffusion<dim>::computeBe(feBilinearForm *form)
{
  double jac, coeff, grad_u[3] = {0., 0., 0.};
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate diffusivity k(t,x)
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    // Compute grad(u)
    form->_intSpaces[_idU]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], k, form->_transformation, grad_u);

    // Compute grad(phi)
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhi.data());

    for(int i = 0; i < _nFunctions; ++i) {
      for(int iDim = 0; iDim < dim; ++iDim) {
        form->_Be[i] += coeff * _gradPhi[i * dim + iDim] * grad_u[iDim] * jac * _wQuad[k];
      }
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: diffusion of scalar field with nonlinear diffusivity
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_NonlinearDiffusion<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _phiU.resize(_nFunctions);
  _gradPhi.resize(dim * _nFunctions);
}

template <int dim> void feSysElm_NonlinearDiffusion<dim>::computeAe(feBilinearForm *form)
{
  double grad_u[3] = {0., 0., 0.};
  for(int k = 0; k < _nQuad; ++k) {
    double jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate diffusivity k(u) and derivative dkdu
    double u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    double ku = (*_diffusivity)(u);
    double dkdu = (*_ddiffdu)(u);

    // Compute grad(u)
    form->_intSpaces[_idU]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], k, form->_transformation, grad_u);

    // Get phi
    form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);

    // Get grad(phi)
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhi.data());

    for(int i = 0; i < _nFunctions; ++i) {
      for(int j = 0; j < _nFunctions; ++j) {
        for(int iDim = 0; iDim < dim; ++iDim) {
          // Linear part
          form->_Ae[i][j] -=
            ku * _gradPhi[i * dim + iDim] * _gradPhi[j * dim + iDim] * jac * _wQuad[k];
          // Nonlinear part
          form->_Ae[i][j] -=
            dkdu * _phiU[j] * grad_u[iDim] * _gradPhi[i * dim + iDim] * jac * _wQuad[k];
        }
      }
    }
  }
}

template <int dim> void feSysElm_NonlinearDiffusion<dim>::computeBe(feBilinearForm *form)
{
  double grad_u[3] = {0., 0., 0.};
  for(int k = 0; k < _nQuad; ++k) {
    double jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate diffusivity k(u)
    double u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    double kD = (*_diffusivity)(u);

    // Compute grad(u)
    form->_intSpaces[_idU]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], k, form->_transformation, grad_u);

    // Compute grad(phi)
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhi.data());

    for(int i = 0; i < _nFunctions; ++i) {
      for(int iDim = 0; iDim < dim; ++iDim) {
        form->_Be[i] += _gradPhi[i * dim + iDim] * grad_u[iDim] * kD * jac * _wQuad[k];
      }
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: advection of scalar field
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_Advection<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _phiU.resize(_nFunctions);
  _gradPhi.resize(dim * _nFunctions);
}

template <int dim> void feSysElm_Advection<dim>::computeAe(feBilinearForm *form)
{
  std::vector<double> c(dim, 0.0);
  for(int k = 0; k < _nQuad; ++k) {
    double jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate the imposed velocity field
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    (*_velocity)(form->_tn, _pos, c);
    feInfo("c = %f", c[0]);

    // Get phi
    form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);

    // Get grad(phi)
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhi.data());

    for(int i = 0; i < _nFunctions; ++i) {
      for(int j = 0; j < _nFunctions; ++j) {
        for(int iDim = 0; iDim < dim; ++iDim) {
          form->_Ae[i][j] -=
            c[iDim] * _gradPhi[i * dim + iDim] * _phiU[j] * jac * _wQuad[k]; // int c*u*dvdx
          // form->_Ae[i][j] += c[iDim] * _gradPhi[j*dim + iDim] * _phiU[i] * jac * _wQuad[k]; //
          // int c*dudx*v
        }
      }
    }
  }
}

template <int dim> void feSysElm_Advection<dim>::computeBe(feBilinearForm *form)
{
  std::vector<double> c(dim, 0.0);
  double grad_u[3] = {0., 0., 0.};
  for(int k = 0; k < _nQuad; ++k) {
    double jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate the imposed velocity field
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    (*_velocity)(form->_tn, _pos, c);

    // Compute grad(u)
    form->_intSpaces[_idU]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], k, form->_transformation, grad_u);

    // Compute grad(phi)
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhi.data());

    double u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    for(int i = 0; i < _nFunctions; ++i) {
      _phiU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);

      for(int iDim = 0; iDim < dim; ++iDim) {
        form->_Be[i] += c[iDim] * u * _gradPhi[i * dim + iDim] * jac * _wQuad[k]; // int c*u*dvdx
        // form->_Be[i] -= c[iDim] * grad_u[iDim] * _phiU[i] * jac * _wQuad[k]; // int c*dudx*v
      }
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: nonlinear advection of scalar field
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_NonlinearAdvection<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _phiU.resize(_nFunctions);
  _gradPhi.resize(dim * _nFunctions);
}

template <int dim> void feSysElm_NonlinearAdvection<dim>::computeAe(feBilinearForm *form)
{
  // std::vector<double> v(dim, 0.0);
  //  for(int k = 0; k < _nQuad; ++k) {
  //    double jac = form->_J[_nQuad * form->_numElem + k];

  //    // Evaluate the imposed velocity field
  //    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
  //    (*_velocity)(form->_tn, _pos, v);

  //    // Get phi
  //    form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);

  //    // Get grad(phi)
  //    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord,
  //    	form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz,
  //    _gradPhi.data());

  //    for(int i = 0; i < _nFunctions; ++i) {
  //      for(int j = 0; j < _nFunctions; ++j) {
  //        for(int iDim = 0; iDim < dim; ++iDim){
  //          form->_Ae[i][j] += v[iDim] * _gradPhi[j*dim + iDim] * _phiU[i] * jac * _wQuad[k];
  //        }
  //      }
  //    }
  //  }
}

template <int dim> void feSysElm_NonlinearAdvection<dim>::computeBe(feBilinearForm *form)
{
  std::vector<double> f(dim, 0.0);
  for(int k = 0; k < _nQuad; ++k) {
    double jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate the nonlinear flux
    double u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    (*_flux)(u, _pos, f);

    // Compute grad(phi)
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhi.data());

    for(int i = 0; i < _nFunctions; ++i) {
      for(int iDim = 0; iDim < dim; ++iDim) {
        form->_Be[i] += f[iDim] * _gradPhi[i * dim + iDim] * jac * _wQuad[k];
      }
    }
  }
}

// -----------------------------------------------------------------------------
// LDG
// -----------------------------------------------------------------------------
void feSysElm_LDG1::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idQ = 1;
  _fieldsLayoutI = {_idU, _idQ};
  _fieldsLayoutJ = {_idU, _idQ};

  // _nFunctions = space[_idU]->getNumFunctions();
  // _gradPhi.resize(_nFunctions);

  _nFunctionsQ = space[_idQ]->getNumFunctions();
  _nFunctionsU = space[_idU]->getNumFunctions();
  _phiQ.resize(_nFunctionsQ);
  _phiU.resize(_nFunctionsU);
  _gradPhiQ.resize(_nFunctionsQ);
  _gradPhiU.resize(_nFunctionsU);
}

void feSysElm_LDG1::computeAe(feBilinearForm *form)
{
  // int I = 0, J = 0;
  // for(int i = 0; i < _nFunctionsU; ++i) {
  //   J = 0;

  //   // Bloc test u-u : Skip
  //   for(int j = 0; j < _nFunctionsQ; ++j) {
  //     double hu_u = c * _phiU[j];
  //     form->_Ae[I][J++] -= hu_u * _gradPhiU[i];
  //   }

  //   // Bloc test u-q
  //   for(int j = 0; j < _nFunctionsQ; ++j) {
  //     double hu_q = - sqrt(kd) * _phiQ[j];
  //     form->_Ae[I][J++] -= hu_q * _gradPhiU[i];
  //   }

  //   I++;
  // }

  // for(int i = 0; i < _nFunctionsQ; ++i) {
  //   J = 0;

  //   // Bloc test q-u
  //   for(int j = 0; j < _nFunctionsU; ++j) {
  //     double hq = - sqrt(kd) * _phiU[j];
  //     form->_Ae[I][J++] -= hq * _gradPhiQ[i];
  //   }

  //   // Bloc test q-q
  //   for(int j = 0; j < _nFunctionsQ; ++j) {
  //     // Skip
  //     J++;
  //   }

  //   I++;
  // }
}

void feSysElm_LDG1::computeBe(feBilinearForm *form)
{
  // double d2 = 1.;

  double rLeft[3]  = {-1.0, 0., 0.};
  double rRight[3] = {+1.0, 0., 0.};

  // std::vector<double> phiL(_nFunctions); // Phi+(x j-1/2) (left)
  // std::vector<double> phiR(_nFunctions); // Phi-(x j+1/2) (right)
  // form->_intSpaces[_idU]->L(rLeft, phiL.data());
  // form->_intSpaces[_idU]->L(rRight, phiR.data());

  // form->_geoSpace->interpolateVectorField(form->_geoCoord, rLeft, _pos);
  // double fL = (*_source)(form->_tn, _pos);
  // form->_geoSpace->interpolateVectorField(form->_geoCoord, rRight, _pos);
  // double fR = (*_source)(form->_tn, _pos);

  // double uL, uR, uPrevR, uNextL;
  // uL = form->_intSpaces[_idU]->interpolateField(form->_sol[_idU], rLeft); // u+(j-1/2)
  // uR = form->_intSpaces[_idU]->interpolateField(form->_sol[_idU], rRight); // u-(j+1/2)
  // uPrevR = form->_intSpaces[_idU]->interpolateField(form->_solPrev[_idU], rRight); // u-(j-1/2)
  // uNextL = form->_intSpaces[_idU]->interpolateField(form->_solNext[_idU], rLeft); // u+(j+1/2)

  // for(int i = 0; i < _nFunctions; ++i) {
  //   feInfo("Adding %f and %f on elm %d", uR - uNextL, uPrevR - uL, form->_numElem);
  //   form->_Be[i] -= - d2 * (uR - uNextL) * phiR[i] - d2 * (uPrevR - uL) * phiL[i];
  // }

  // =================================================
  // L2 Projection
  // double src;
  // for(int k = 0; k < _nQuad; ++k) {
  //   double jac = form->_J[_nQuad * form->_numElem + k];
  //   form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

  //   // Source
  //   form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
  //   src = (*_source)(form->_tn, _pos);

  //   // Compute grad(phi)
  //   form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
  //                                                                  _gradPhi.data());

  //   for(int i = 0; i < _nFunctions; ++i) {
  //     form->_Be[i] += src * _gradPhi[i] * jac * _wQuad[k];
  //   }
  // }

  // // DG term
  // for(int i = 0; i < _nFunctions; ++i) {
  //   form->_Be[i] -= fR * phiR[i] - fL * phiL[i];
  // }
  // =================================================

  //==================================================
  // Diffusion
  double kd = 1.0;
  double c = _velocity;
  // double c12 = -sqrt(kd)/2.;
  double c12 = 0.;

  //==================================================================
  // Decoupled problem
  // // DG jumps, averages and fluxes
  // double uL, uR, uPrevR, uNextL;
  // uL = form->_intSpaces[_idU]->interpolateField(form->_sol[_idU], rLeft); // u+(j-1/2)
  // uR = form->_intSpaces[_idU]->interpolateField(form->_sol[_idU], rRight); // u-(j+1/2)
  // uPrevR = form->_intSpaces[_idU]->interpolateField(form->_solPrev[_idU], rRight); // u-(j-1/2)
  // uNextL = form->_intSpaces[_idU]->interpolateField(form->_solNext[_idU], rLeft); // u+(j+1/2)

  // std::vector<double> phiuL(_nFunctionsU); // Phi+(x j-1/2) (left)
  // std::vector<double> phiuR(_nFunctionsU); // Phi-(x j+1/2) (right)
  // std::vector<double> phiqL(_nFunctionsQ); // Phi+(x j-1/2) (left)
  // std::vector<double> phiqR(_nFunctionsQ); // Phi-(x j+1/2) (right)
  // form->_intSpaces[_idU]->L(rLeft,  phiuL.data());
  // form->_intSpaces[_idU]->L(rRight, phiuR.data());
  // form->_intSpaces[_idQ]->L(rLeft,  phiqL.data());
  // form->_intSpaces[_idQ]->L(rRight, phiqR.data());

  // double uAvgR = 0.5 * (uR + uNextL);
  // double uAvgL = 0.5 * (uPrevR + uL);
  // double uJumpR = uNextL - uR;
  // double uJumpL = uL - uPrevR;

  // // Choix du flux
  // double huHatL, huHatR, hqHatL, hqHatR;

  // // Lax Friedrichs convective flux
  // double hConvR = c*uAvgR - fabs(c)/2. * uJumpR;
  // double hConvL = c*uAvgL - fabs(c)/2. * uJumpL;

  // // huHatR = hConvR - sqrt(kd) * qAvgR - c12 * qJumpR;
  // // huHatL = hConvL - sqrt(kd) * qAvgL - c12 * qJumpL;
  // hqHatR = 0.     - sqrt(kd) * uAvgR + c12 * uJumpR;
  // hqHatL = 0.     - sqrt(kd) * uAvgL + c12 * uJumpL;

  // // Solve the local problem for q
  // SquareMatrix qMass(_nFunctionsU);
  // Vector qRHS(_nFunctionsU);
  // for(int k = 0; k < _nQuad; ++k) {
  //   double jac = form->_J[_nQuad * form->_numElem + k];
  //   form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

  //   double u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
  //   form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);
  //   form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiU.data());

  //   double hq = - sqrt(kd) * u;

  //   for(int i = 0; i < _nFunctionsU; ++i) {
  //     for(int j = 0; j < _nFunctionsU; ++j) {
  //       qMass(i,j) += _phiU[i] * _phiU[j] * jac * _wQuad[k];
  //     }
  //     qRHS(i) += hq * _gradPhiU[i] * jac * _wQuad[k];
  //   }
  // }

  // for(int i = 0; i < _nFunctionsU; ++i) {
  //   qRHS(i) += hqHatL * phiqL[i] - hqHatR * phiqR[i];
  // }

  // // Solve for q
  // Vector qLoc = qMass.inverse() * qRHS;
  // std::vector<double> qVec(_nFunctionsU);
  // for(int i = 0; i < _nFunctionsU; ++i){
  //   qVec[i] = qLoc(i);
  // }

  // // Compute jumps and averages from qLoc
  // double qL, qR, qPrevR, qNextL;
  // qL     = form->_intSpaces[_idQ]->interpolateField(qVec, rLeft); // u+(j-1/2)
  // qR     = form->_intSpaces[_idQ]->interpolateField(qVec, rRight); // u-(j+1/2)
  // qPrevR = form->_intSpaces[_idQ]->interpolateField(form->_solPrev[_idQ], rRight); // u-(j-1/2)
  // qNextL = form->_intSpaces[_idQ]->interpolateField(form->_solNext[_idQ], rLeft); // u+(j+1/2)

  // double qAvgR = 0.5 * (qR + qNextL);
  // double qAvgL = 0.5 * (qPrevR + qL);
  // double qJumpR = qNextL - qR;
  // double qJumpL = qL - qPrevR;

  // huHatR = hConvR - sqrt(kd) * qAvgR - c12 * qJumpR;
  // huHatL = hConvL - sqrt(kd) * qAvgL - c12 * qJumpL;
  // hqHatR = 0.     - sqrt(kd) * uAvgR + c12 * uJumpR;
  // hqHatL = 0.     - sqrt(kd) * uAvgL + c12 * uJumpL;

  // // Now solve for u
  // for(int k = 0; k < _nQuad; ++k) {
  //   double jac = form->_J[_nQuad * form->_numElem + k];
  //   form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

  //   // Compute q, grad(phi)
  //   double u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
  //   form->_intSpaces[_idQ]->getFunctionsAtQuadNode(k, _phiQ);
  //   form->_intSpaces[_idQ]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiQ.data());
  //   form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiU.data());

  //   // Get q from qLoc, not from the other field
  //   // double q = form->_intSpaces[_idQ]->interpolateFieldAtQuadNode(form->_sol[_idQ], k);
  //   double q = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(qVec, k);

  //   double hu = c * u - sqrt(k) * q;

  //   for(int i = 0; i < _nFunctionsU; ++i) {
  //     form->_Be[i] += hu * _gradPhiU[i] * jac * _wQuad[k];
  //   }
  // }

  // for(int i = 0; i < _nFunctionsU; ++i) {
  //   form->_Be[i] -= huHatR * phiuR[i] - huHatL * phiuL[i];
  // }
  //==================================================================

  // // Coupled problem
  // int cnt;
  // // Integrale de q * vq
  // for(int k = 0; k < _nQuad; ++k) {
  //   double jac = form->_J[_nQuad * form->_numElem + k];
  //   form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

  //   // Compute q, grad(phi)
  //   double q = form->_intSpaces[_idQ]->interpolateFieldAtQuadNode(form->_sol[_idQ], k);
  //   double u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
  //   form->_intSpaces[_idQ]->getFunctionsAtQuadNode(k, _phiQ);
  //   form->_intSpaces[_idQ]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiQ.data());
  //   form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiU.data());

  //   // form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
  //   // double c = (*_velocity)(form->_tn, _pos);

  //   double hu = c * u - sqrt(kd) * q;
  //   double hq = - sqrt(kd) * u;

  //   cnt = 0;
  //   for(int i = 0; i < _nFunctionsU; ++i) {
  //     form->_Be[cnt] += hu * _gradPhiU[i] * jac * _wQuad[k];
  //     cnt++;
  //   }

  //   for(int i = 0; i < _nFunctionsQ; ++i) {
  //     // form->_Be[cnt] -= q * _phiQ[i] * jac * _wQuad[k];
  //     form->_Be[cnt] += hq * _gradPhiQ[i] * jac * _wQuad[k];
  //     cnt++;
  //   }
  // }

  // DG jumps, averages and fluxes
  double uL, uR, uPrevR, uNextL;
  double qL, qR, qPrevR, qNextL;
  uL     = form->_intSpaces[_idU]->interpolateField(form->_sol[_idU], rLeft); // u+(j-1/2)
  uR     = form->_intSpaces[_idU]->interpolateField(form->_sol[_idU], rRight); // u-(j+1/2)
  uPrevR = form->_intSpaces[_idU]->interpolateField(form->_solPrev[_idU], rRight); // u-(j-1/2)
  uNextL = form->_intSpaces[_idU]->interpolateField(form->_solNext[_idU], rLeft); // u+(j+1/2)

  qL     = form->_intSpaces[_idQ]->interpolateField(form->_sol[_idQ], rLeft); // u+(j-1/2)
  qR     = form->_intSpaces[_idQ]->interpolateField(form->_sol[_idQ], rRight); // u-(j+1/2)
  qPrevR = form->_intSpaces[_idQ]->interpolateField(form->_solPrev[_idQ], rRight); // u-(j-1/2)
  qNextL = form->_intSpaces[_idQ]->interpolateField(form->_solNext[_idQ], rLeft); // u+(j+1/2)

  std::vector<double> phiuL(_nFunctionsU); // Phi+(x j-1/2) (left)
  std::vector<double> phiuR(_nFunctionsU); // Phi-(x j+1/2) (right)
  std::vector<double> phiqL(_nFunctionsQ); // Phi+(x j-1/2) (left)
  std::vector<double> phiqR(_nFunctionsQ); // Phi-(x j+1/2) (right)
  form->_intSpaces[_idU]->L(rLeft,  phiuL.data());
  form->_intSpaces[_idU]->L(rRight, phiuR.data());
  form->_intSpaces[_idQ]->L(rLeft,  phiqL.data());
  form->_intSpaces[_idQ]->L(rRight, phiqR.data());

  double uAvgR = 0.5 * (uR + uNextL);
  double uAvgL = 0.5 * (uPrevR + uL);
  double uJumpR = uNextL - uR;
  double uJumpL = uL - uPrevR;

  double qAvgR = 0.5 * (qR + qNextL);
  double qAvgL = 0.5 * (qPrevR + qL);
  double qJumpR = qNextL - qR;
  double qJumpL = qL - qPrevR;

  // Choix du flux
  double huHatL, huHatR, hqHatL, hqHatR;

  // Lax Friedrichs convective flux
  double hConvR = c*uAvgR - fabs(c)/2. * uJumpR;
  double hConvL = c*uAvgL - fabs(c)/2. * uJumpL;

  // huHatR = hConvR - sqrt(kd) * qAvgR - c12 * qJumpR;
  // huHatL = hConvL - sqrt(kd) * qAvgL - c12 * qJumpL;
  // hqHatR = 0.     - sqrt(kd) * uAvgR + c12 * uJumpR;
  // hqHatL = 0.     - sqrt(kd) * uAvgL + c12 * uJumpL;

  // Alternating q- / u+
  huHatR = hConvR - sqrt(kd) * qR;
  huHatL = hConvL - sqrt(kd) * qPrevR;
  hqHatR = 0.     - sqrt(kd) * uNextL;
  hqHatL = 0.     - sqrt(kd) * uL;

  // Alternating q+ / u-
  // huHatR = hConvR - sqrt(kd) * qNextL;
  // huHatL = hConvL - sqrt(kd) * qL;
  // hqHatR = 0.     - sqrt(kd) * uR;
  // hqHatL = 0.     - sqrt(kd) * uPrevR;

  int cnt = 0;
  for(int i = 0; i < _nFunctionsU; ++i) {
    form->_Be[cnt++] -= (huHatR * phiuR[i] - huHatL * phiuL[i]);
  }

  for(int i = 0; i < _nFunctionsQ; ++i) {
    form->_Be[cnt++] -= (hqHatR * phiqR[i] - hqHatL * phiqL[i]);
  }

  //-----------------------

  // // Integrale de u*wx et de qh*w
  // int cntQ = cnt;
  // for(int k = 0; k < _nQuad; ++k) {
  //   double jac = form->_J[_nQuad * form->_numElem + k];
  //   form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

  //   // Compute u, q, phiU, grad(phiU)
  //   double u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
  //   double q = form->_intSpaces[_idQ]->interpolateFieldAtQuadNode(form->_sol[_idQ], k);
  //   form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);
  //   form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiU.data());

  //   cnt = _nFunctionsQ;
  //   for(int i = 0; i < _nFunctionsU; ++i) {
  //     form->_Be[cnt] -= q * _phiU[i]     * jac * _wQuad[k];
  //     form->_Be[cnt] -= u * _gradPhiU[i] * jac * _wQuad[k];
  //     cnt++;
  //   }
  // }

  // // DG term for u
  // double uHatL, uHatR, uL, uR, uPrevR, uNextL;
  // uL = form->_intSpaces[_idU]->interpolateField(form->_sol[_idU], rLeft); // u+(j-1/2)
  // uR = form->_intSpaces[_idU]->interpolateField(form->_sol[_idU], rRight); // u-(j+1/2)
  // uPrevR = form->_intSpaces[_idU]->interpolateField(form->_solPrev[_idU], rRight); // u-(j-1/2)
  // uNextL = form->_intSpaces[_idU]->interpolateField(form->_solNext[_idU], rLeft); // u+(j+1/2)

  // // Choix du flux
  // uHatL = (uL + uPrevR)/2.;
  // uHatR = (uR + uNextL)/2.;

  // std::vector<double> phiuL(_nFunctions); // Phi+(x j-1/2) (left)
  // std::vector<double> phiuR(_nFunctions); // Phi-(x j+1/2) (right)
  // form->_intSpaces[_idU]->L(rLeft, phiuL.data());
  // form->_intSpaces[_idU]->L(rRight, phiuR.data());

  // cnt = cntQ;
  // for(int i = 0; i < _nFunctionsU; ++i) {
  //   indices.insert(cnt);
  //   form->_Be[cnt++] -= uHatR * phiuR[i] - uHatL * phiuL[i];
  // }

  //==================================================
}

void feSysElm_LDG2::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idQ = 1;
  _fieldsLayoutI = {_idU, _idQ};
  _fieldsLayoutJ = {_idU, _idQ};

  // _nFunctions = space[_idU]->getNumFunctions();
  // _gradPhi.resize(_nFunctions);

  _nFunctionsQ = space[_idQ]->getNumFunctions();
  _nFunctionsU = space[_idU]->getNumFunctions();
  _phiQ.resize(_nFunctionsQ);
  _phiU.resize(_nFunctionsU);
  _gradPhiQ.resize(_nFunctionsQ);
  _gradPhiU.resize(_nFunctionsU);
}

void feSysElm_LDG2::computeAe(feBilinearForm *form)
{
  double kd = 1.0;
  double c = _velocity;

  for(int k = 0; k < _nQuad; ++k) {
    double jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Compute q, grad(phi)
    double q = form->_intSpaces[_idQ]->interpolateFieldAtQuadNode(form->_sol[_idQ], k);
    double u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    form->_intSpaces[_idQ]->getFunctionsAtQuadNode(k, _phiQ);
    form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);
    form->_intSpaces[_idQ]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiQ.data());
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiU.data());

    // double hu = c * u - sqrt(kd) * q;
    // double hq = - sqrt(kd) * u;

    int I = 0, J = 0;
    for(int i = 0; i < _nFunctionsU; ++i) {
      J = 0;

      // Bloc test u-u : Skip
      for(int j = 0; j < _nFunctionsQ; ++j) {
        double hu_u = c * _phiU[j];
        form->_Ae[I][J++] -= hu_u * _gradPhiU[i] * jac * _wQuad[k];
      }

      // Bloc test u-q
      for(int j = 0; j < _nFunctionsQ; ++j) {
        double hu_q = - sqrt(kd) * _phiQ[j];
        form->_Ae[I][J++] -= hu_q * _gradPhiU[i] * jac * _wQuad[k];
      }

      I++;
    }

    for(int i = 0; i < _nFunctionsQ; ++i) {
      J = 0;

      // Bloc test q-u
      for(int j = 0; j < _nFunctionsU; ++j) {
        double hq = - sqrt(kd) * _phiU[j];
        form->_Ae[I][J++] -= hq * _gradPhiQ[i] * jac * _wQuad[k];
      }

      // Bloc test q-q
      for(int j = 0; j < _nFunctionsQ; ++j) {
        // Skip
        J++;
      }

      I++;
    }
  }
}

void feSysElm_LDG2::computeBe(feBilinearForm *form)
{

  double rLeft[3]  = {-1.0, 0., 0.};
  double rRight[3] = {+1.0, 0., 0.};

  //==================================================
  // Diffusion
  double kd = 1.0;
  double c = _velocity;
  // double c12 = -sqrt(kd)/2.;
  double c12 = 0.;

  // Coupled problem
  int cnt;
  for(int k = 0; k < _nQuad; ++k) {
    double jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Compute q, grad(phi)
    double q = form->_intSpaces[_idQ]->interpolateFieldAtQuadNode(form->_sol[_idQ], k);
    double u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    form->_intSpaces[_idQ]->getFunctionsAtQuadNode(k, _phiQ);
    form->_intSpaces[_idQ]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiQ.data());
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiU.data());

    double hu = c * u - sqrt(kd) * q;
    double hq = - sqrt(kd) * u;

    cnt = 0;
    for(int i = 0; i < _nFunctionsU; ++i) {
      form->_Be[cnt] += hu * _gradPhiU[i] * jac * _wQuad[k];
      cnt++;
    }

    for(int i = 0; i < _nFunctionsQ; ++i) {
      // form->_Be[cnt] -= q * _phiQ[i] * jac * _wQuad[k];
      form->_Be[cnt] += hq * _gradPhiQ[i] * jac * _wQuad[k];
      cnt++;
    }
  }
}

// Explicit instantiation of templated forms
template class feSysElm_Diffusion<0>;
template class feSysElm_Diffusion<1>;
template class feSysElm_Diffusion<2>;
template class feSysElm_NonlinearDiffusion<0>;
template class feSysElm_NonlinearDiffusion<1>;
template class feSysElm_NonlinearDiffusion<2>;
template class feSysElm_Advection<0>;
template class feSysElm_Advection<1>;
template class feSysElm_Advection<2>;
template class feSysElm_NonlinearAdvection<0>;
template class feSysElm_NonlinearAdvection<1>;
template class feSysElm_NonlinearAdvection<2>;