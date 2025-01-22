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
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    double S = (*_source)(form->_tn, _pos);
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

  //   form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
  //   coeff = (*_coeff)(form->_tn, _pos);

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

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

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

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    for(int i = 0; i < _nFunctionsU; ++i)
      _phiU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
    for(int i = 0; i < _nFunctionsV; ++i)
      _phiV[i] = form->_intSpaces[_idV]->getFunctionAtQuadNode(i, k);

    for(int i = 0; i < _nFunctionsV; ++i)
      for(int j = 0; j < _nFunctionsU; ++j)
        form->_Ae[i][j] += coeff * _phiV[i] * _phiU[j] * jac * _wQuad[k];
  }
}

void feSysElm_MixedMass::computeBe(feBilinearForm *form)
{
  double jac, u, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

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
  // Attention à la dérivation avec la puissance
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Implement matrix for MIXEDMASSPOWER");
  // double jac, coeff;
  // for(int k = 0; k < _nQuad; ++k) {
  //   jac = form->_J[_nQuad * form->_numElem + k];

  //   form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
  //   coeff = (*_coeff)(form->_tn, _pos);

  //   for(int i = 0; i < _nFunctionsU; ++i)
  //     _phiU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
  //   for(int i = 0; i < _nFunctionsV; ++i)
  //     _phiV[i] = form->_intSpaces[_idV]->getFunctionAtQuadNode(i, k);

  //   for(int i = 0; i < _nFunctionsV; ++i)
  //     for(int j = 0; j < _nFunctionsU; ++j)
  //       form->_Ae[i][j] += coeff * _phiV[i] * _phiU[j] * jac * _wQuad[k];
  // }
}

void feSysElm_MixedMassPower::computeBe(feBilinearForm *form)
{
  double jac, u, up, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

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

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    for(int i = 0; i < _nFunctions; ++i)
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);

    for(int i = 0; i < _nFunctions; ++i)
      for(int j = 0; j < _nFunctions; ++j)
        form->_Ae[i][j] += _feU[i] * coeff * form->_c0 * _feU[j] * jac * _wQuad[k];
  }
}

void feSysElm_TransientMass::computeBe(feBilinearForm *form)
{
  double jac, uDot, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

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

  //   form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
  //   coeff = (*_coeff)(form->_tn, _pos);

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

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

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
        form->_Be[i] -= _gradPhi[i * dim + iDim] * grad_u[iDim] * kD * jac * _wQuad[k];
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
  _gradPhiV.resize(dim * _nFunctionsV);
}

template <int dim>
void feSysElm_MixedGradGrad<dim>::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Implement matrix for MixedGradGrad");
  // double jac, coeff, div_v;
  // for(int k = 0; k < _nQuad; ++k) {
  //   jac = form->_J[_nQuad * form->_numElem + k];
  //   form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

  //   // Evaluate scalar coefficient
  //   form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
  //   coeff = (*_coeff)(form->_tn, _pos);

  //   // Get phiU
  //   form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);

  //   // Compute compacted grad(phiV) without the trivial zeros
  //   form->_intSpaces[_idV]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
  //                                                                  _gradPhiV.data());

  //   for(int i = 0; i < _nFunctionsV; ++i) {
  //     for(int j = 0; j < _nFunctionsU; ++j) {
  //       div_v = _gradPhiV[i * _nComponents + (i % _nComponents)];
  //       form->_Ae[i][j] -= coeff * _phiU[j] * div_v * jac * _wQuad[k];
  //     }
  //   }
  // }
}

template <int dim>
void feSysElm_MixedGradGrad<dim>::computeBe(feBilinearForm *form)
{
  double jac, coeff, grad_u[3] = {0., 0., 0.};
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    // Get gradU and gradPhiV
    form->_intSpaces[_idU]->interpolateFieldAtQuadNode_physicalGradient(form->_sol[_idU], k, form->_transformation, grad_u);
    form->_intSpaces[_idV]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiV.data());

    for(int i = 0; i < _nFunctionsV; ++i) {
      for(int iDim = 0; iDim < dim; ++iDim) {
        form->_Be[i] -= coeff * _gradPhiV[i * dim + iDim] * grad_u[iDim] * jac * _wQuad[k];
      }
    }
  }
}

template class feSysElm_MixedGradGrad<1>;
template class feSysElm_MixedGradGrad<2>;
template class feSysElm_MixedGradGrad<3>;