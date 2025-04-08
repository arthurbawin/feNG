#include "feSysElm.h"
#include "feBilinearForm.h"

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
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    double S = (*_source)(form->_args);
    for(int i = 0; i < _nFunctions; ++i) {
      _phiU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= _phiU[i] * S * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Linear form: scalar source term where the source term is the Dirac delta
// -----------------------------------------------------------------------------
void feSysElm_SourceDirac::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _phiAtX0.resize(_nFunctions);
}

void feSysElm_SourceDirac::computeBe(feBilinearForm *form)
{
  #if defined(HAVE_OMP)
  #pragma omp critical
  #endif
  {
    // Locate x0 in the element
    int elm;
    double xsi[3];
    bool isFound = form->_cnc->getMeshPtr()->locateVertex(_x0.data(), elm, xsi, 1e-12);

    if(!isFound)
      feWarning("Could not locate Dirac vertex (%+-1.4e - %+-1.4e) in the mesh!", _x0[0], _x0[1]);

    if(isFound && elm == form->_numElem) {
      // Evaluate shape functions at x0
      form->_intSpaces[_idU]->L(xsi, _phiAtX0.data());
      for(int i = 0; i < _nFunctions; ++i) {
        form->_Be[i] -= _phiAtX0[i];
      }
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

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    for(int i = 0; i < _nFunctions; ++i)
      _phiU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);

    for(int i = 0; i < _nFunctions; ++i)
      for(int j = 0; j < _nFunctions; ++j)
        form->_Ae[i][j] += coeff * _phiU[i] * _phiU[j] * jac * _wQuad[k];
  }
  // if(form->_numElem == 10) {
  //   feInfo("%s - Local matrix on elem 10 :", toString(_ID).data());
  //   for(int i = 0; i < form->_M; ++i)
  //     for(int j = 0; j < form->_N; ++j)
  //       feInfo("A[%2d][%2d] = %+-1.10e", i, j, form->_Ae[i][j]);
  // }
}

void feSysElm_Mass::computeBe(feBilinearForm *form)
{
  double jac, u, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);

    for(int i = 0; i < _nFunctions; ++i) {
      _phiU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= coeff * u * _phiU[i] * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: mass matrix for scalar field to given power
// -----------------------------------------------------------------------------
void feSysElm_MassPower::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _phiU.resize(_nFunctions);
}

void feSysElm_MassPower::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Implement matrix for MASSPOWER");
  // double jac, coeff;
  // for(int k = 0; k < _nQuad; ++k) {
  //   jac = form->_J[_nQuad * form->_numElem + k];

  //   form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
  //   coeff = (*_coeff)(form->_args);

  //   for(int i = 0; i < _nFunctions; ++i)
  //     _phiU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);

  //   for(int i = 0; i < _nFunctions; ++i)
  //     for(int j = 0; j < _nFunctions; ++j)
  //       form->_Ae[i][j] += coeff * _phiU[i] * _phiU[j] * jac * _wQuad[k];
  // }
}

void feSysElm_MassPower::computeBe(feBilinearForm *form)
{
  double jac, u, uP, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    uP = pow(u, _p);

    for(int i = 0; i < _nFunctions; ++i) {
      _phiU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= coeff * uP * _phiU[i] * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: mixed mass matrix for scalar field (unknown and tests functions
//                from different variables)
// -----------------------------------------------------------------------------
void feSysElm_MixedMass::createElementarySystem(std::vector<feSpace *> &space)
{
  _idV = 0;
  _idU = 1;
  _fieldsLayoutI = {_idV}; // Rectangular local matrix
  _fieldsLayoutJ = {_idU};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsV = space[_idV]->getNumFunctions();
  _phiU.resize(_nFunctionsU);
  _phiV.resize(_nFunctionsV);
}

void feSysElm_MixedMass::computeAe(feBilinearForm *form)
{
  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    for(int i = 0; i < _nFunctionsU; ++i)
      _phiU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
    for(int i = 0; i < _nFunctionsV; ++i)
      _phiV[i] = form->_intSpaces[_idV]->getFunctionAtQuadNode(i, k);

    for(int i = 0; i < _nFunctionsV; ++i)
      for(int j = 0; j < _nFunctionsU; ++j)
        form->_Ae[i][j] += coeff * _phiV[i] * _phiU[j] * jac * _wQuad[k];
  }
  // if(form->_numElem == 10) {
  //   feInfo("%s - Local matrix on elem 10 :", toString(_ID).data());
  //   for(int i = 0; i < form->_M; ++i)
  //     for(int j = 0; j < form->_N; ++j)
  //       feInfo("A[%2d][%2d] = %+-1.10e", i, j, form->_Ae[i][j]);
  // }
}

void feSysElm_MixedMass::computeBe(feBilinearForm *form)
{
  double jac, u, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);

    for(int i = 0; i < _nFunctionsV; ++i) {
      _phiV[i] = form->_intSpaces[_idV]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= coeff * u * _phiV[i] * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: mixed mass power
// -----------------------------------------------------------------------------
void feSysElm_MixedMassPower::createElementarySystem(std::vector<feSpace *> &space)
{
  _idV = 0;
  _idU = 1;
  _fieldsLayoutI = {_idV}; // Rectangular local matrix
  _fieldsLayoutJ = {_idU};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsV = space[_idV]->getNumFunctions();
  _phiU.resize(_nFunctionsU);
  _phiV.resize(_nFunctionsV);
}

void feSysElm_MixedMassPower::computeAe(feBilinearForm *form)
{
  double jac, coeff, u, upm1;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    upm1 = pow(u, _p - 1);
    form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);
    form->_intSpaces[_idV]->getFunctionsAtQuadNode(k, _phiV);

    for(int i = 0; i < _nFunctionsV; ++i)
      for(int j = 0; j < _nFunctionsU; ++j)
        form->_Ae[i][j] += coeff * _p * upm1 * _phiV[i] * _phiU[j] * jac * _wQuad[k];
  }
  // if(form->_numElem == 10) {
  //   feInfo("%s - Local matrix on elem 10 :", toString(_ID).data());
  //   for(int i = 0; i < form->_M; ++i)
  //     for(int j = 0; j < form->_N; ++j)
  //       feInfo("A[%2d][%2d] = %+-1.10e", i, j, form->_Ae[i][j]);
  // }
}

void feSysElm_MixedMassPower::computeBe(feBilinearForm *form)
{
  double jac, u, up, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    up = pow(u, _p);
    form->_intSpaces[_idV]->getFunctionsAtQuadNode(k, _phiV);

    for(int i = 0; i < _nFunctionsV; ++i) {
      form->_Be[i] -= coeff * up * _phiV[i] * jac * _wQuad[k];
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
  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    for(int i = 0; i < _nFunctions; ++i)
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);

    for(int i = 0; i < _nFunctions; ++i)
      for(int j = 0; j < _nFunctions; ++j)
        form->_Ae[i][j] += _feU[i] * coeff * form->_c0 * _feU[j] * jac * _wQuad[k];
  }
  // if(form->_numElem == 10) {
  //   feInfo("%s - Local matrix on elem 10 :", toString(_ID).data());
  //   for(int i = 0; i < form->_M; ++i)
  //     for(int j = 0; j < form->_N; ++j)
  //       feInfo("A[%2d][%2d] = %+-1.10e", i, j, form->_Ae[i][j]);
  // }
}

void feSysElm_TransientMass::computeBe(feBilinearForm *form)
{
  double jac, uDot, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    uDot = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_solDot[_idU], k);

    for(int i = 0; i < _nFunctions; ++i) {
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= _feU[i] * coeff * uDot * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: mixed transient mass matrix
// -----------------------------------------------------------------------------
void feSysElm_MixedTransientMass::createElementarySystem(std::vector<feSpace *> &space)
{
  _idV = 0;
  _idU = 1;
  _fieldsLayoutI = {_idV}; // Rectangular local matrix
  _fieldsLayoutJ = {_idU};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsV = space[_idV]->getNumFunctions();
  _phiU.resize(_nFunctionsU);
  _phiV.resize(_nFunctionsV);
}

void feSysElm_MixedTransientMass::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Implement matrix for MixedTransientMass");
  // double jac, coeff;
  // for(int k = 0; k < _nQuad; ++k) {
  //   jac = form->_J[_nQuad * form->_numElem + k];

  //   form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
  //   coeff = (*_coeff)(form->_args);

  //   for(int i = 0; i < _nFunctions; ++i)
  //     _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);

  //   for(int i = 0; i < _nFunctions; ++i)
  //     for(int j = 0; j < _nFunctions; ++j)
  //       form->_Ae[i][j] += _feU[i] * coeff * form->_c0 * _feU[j] * jac * _wQuad[k];
  // }
}

void feSysElm_MixedTransientMass::computeBe(feBilinearForm *form)
{
  double jac, uDot, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    // Compute uDot and get phiV
    uDot = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_solDot[_idU], k);
    form->_intSpaces[_idV]->getFunctionsAtQuadNode(k, _phiV);

    for(int i = 0; i < _nFunctionsV; ++i) {
      form->_Be[i] -= coeff * uDot * _phiV[i] * jac * _wQuad[k];
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
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar parameter
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    // Compute grad(phi)
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhi.data());

    for(int i = 0; i < _nFunctions; ++i) {
      for(int j = 0; j < _nFunctions; ++j) {
        for(int iDim = 0; iDim < dim; ++iDim) {
          form->_Ae[i][j] +=
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
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate diffusivity k(t,x)
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    // Compute grad(u)
    form->_intSpaces[_idU]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], k, form->_transformation, grad_u);

    // Compute grad(phi)
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhi.data());

    for(int i = 0; i < _nFunctions; ++i) {
      for(int iDim = 0; iDim < dim; ++iDim) {
        form->_Be[i] -= coeff * _gradPhi[i * dim + iDim] * grad_u[iDim] * jac * _wQuad[k];
      }
    }
  }
}

// Explicit instantiation
template class feSysElm_Diffusion<0>;
template class feSysElm_Diffusion<1>;
template class feSysElm_Diffusion<2>;
template class feSysElm_Diffusion<3>;

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
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate diffusivity k(u) and derivative dkdu
    double u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    form->_args.u = u;
    double ku = (*_diffusivity).eval(form->_args);
    double dkdu = (*_ddiffdu).eval(form->_args);

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
          form->_Ae[i][j] +=
            ku * _gradPhi[i * dim + iDim] * _gradPhi[j * dim + iDim] * jac * _wQuad[k];
          // Nonlinear part
          form->_Ae[i][j] +=
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
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate diffusivity k(u)
    double u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    form->_args.u = u;
    double ku = (*_diffusivity).eval(form->_args);

    // Compute grad(u)
    form->_intSpaces[_idU]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], k, form->_transformation, grad_u);

    // Compute grad(phi)
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhi.data());

    for(int i = 0; i < _nFunctions; ++i) {
      for(int iDim = 0; iDim < dim; ++iDim) {
        form->_Be[i] -= _gradPhi[i * dim + iDim] * grad_u[iDim] * ku * jac * _wQuad[k];
      }
    }
  }
}

template class feSysElm_NonlinearDiffusion<1>;
template class feSysElm_NonlinearDiffusion<2>;
template class feSysElm_NonlinearDiffusion<3>;

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
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate the imposed velocity field
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    (*_velocity)(form->_args, c);
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
          // form->_Ae[i][j] += c[iDim] * _gradPhi[j*dim + iDim] * _phiU[i] * jac * _wQuad[k]; // int c*dudx*v
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
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate the imposed velocity field
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    (*_velocity)(form->_args, c);

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
        // UNUSED(u);
        // form->_Be[i] -= c[iDim] * grad_u[iDim] * _phiU[i] * jac * _wQuad[k]; // int c*dudx*v
      }
    }
  }
}

template class feSysElm_Advection<1>;
template class feSysElm_Advection<2>;
template class feSysElm_Advection<3>;

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
  UNUSED(form);
  // To implement
}

template <int dim> void feSysElm_NonlinearAdvection<dim>::computeBe(feBilinearForm *form)
{
  std::vector<double> f(dim, 0.0);
  for(int k = 0; k < _nQuad; ++k) {
    double jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate the nonlinear flux
    form->_args.u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    (*_flux)(form->_args, f);

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

template class feSysElm_NonlinearAdvection<1>;
template class feSysElm_NonlinearAdvection<2>;
template class feSysElm_NonlinearAdvection<3>;

// -----------------------------------------------------------------------------
// Bilinear form: mixed gradient-gradient
// -----------------------------------------------------------------------------
template <int dim>
void feSysElm_MixedGradGrad<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idV = 0; // Scalar test functions
  _idU = 1; // Scalar unknown field
  _fieldsLayoutI = {_idV}; // Rectangular local matrix
  _fieldsLayoutJ = {_idU};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsV = space[_idV]->getNumFunctions();
  _gradPhiU.resize(dim * _nFunctionsU);
  _gradPhiV.resize(dim * _nFunctionsV);
}

template <int dim>
void feSysElm_MixedGradGrad<dim>::computeAe(feBilinearForm *form)
{
  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar parameter
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    // Get gradPhiU and gradPhiV
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiU.data());
    form->_intSpaces[_idV]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiV.data());

    for(int i = 0; i < _nFunctionsV; ++i) {
      for(int j = 0; j < _nFunctionsU; ++j) {
        for(int iDim = 0; iDim < dim; ++iDim) {
          form->_Ae[i][j] +=
            coeff * _gradPhiV[i * dim + iDim] * _gradPhiU[j * dim + iDim] * jac * _wQuad[k];
        }
      }
    }
  }
  // if(form->_numElem == 10) {
  //   feInfo("%s - Local matrix on elem 10 :", toString(_ID).data());
  //   for(int i = 0; i < form->_M; ++i)
  //     for(int j = 0; j < form->_N; ++j)
  //       feInfo("A[%2d][%2d] = %+-1.10e", i, j, form->_Ae[i][j]);
  // }
}

template <int dim>
void feSysElm_MixedGradGrad<dim>::computeBe(feBilinearForm *form)
{
  double jac, coeff, grad_u[3] = {0., 0., 0.};
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->getElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, form->_args.pos);
    coeff = (*_coeff)(form->_args);

    // Get gradU and gradPhiV
    form->_intSpaces[_idU]->interpolateFieldAtQuadNode_physicalGradient(form->_sol[_idU], k, form->_transformation, grad_u);
    form->_intSpaces[_idV]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiV.data());

    for(int i = 0; i < _nFunctionsV; ++i) {
      
      double dotprod = 0.;
      for(int m = 0; m < dim; ++m) {
        dotprod += grad_u[m] * _gradPhiV[i*dim + m];
      }

      form->_Be[i] -= coeff * dotprod * jac * _wQuad[k];
    }
  }
  // if(form->_numElem == 10) {
  //   feInfo("%s - Local RHS on elem 10 :", toString(_ID).data());
  //   for(int i = 0; i < form->_M; ++i)
  //     feInfo("B[%2d] = %+-1.10e", i, form->_Be[i]);
  // }
}

template class feSysElm_MixedGradGrad<1>;
template class feSysElm_MixedGradGrad<2>;
template class feSysElm_MixedGradGrad<3>;