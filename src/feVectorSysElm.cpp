#include "feSysElm.h"
#include "feBilinearForm.h"
#include "feNumeric.h"

// -----------------------------------------------------------------------------
// Linear form: vector-valued source term
// -----------------------------------------------------------------------------
void feSysElm_VectorSource::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _phiU.resize(_nFunctions);
}

void feSysElm_VectorSource::computeBe(feBilinearForm *form)
{
  for(int k = 0; k < _nQuad; ++k) {
    double jac = form->_J[_nQuad * form->_numElem + k];
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    (*_source)(form->_tn, _pos, _S);
    for(int i = 0; i < _nFunctions; ++i) {
      _phiU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= _phiU[i] * _S[i % _nComponents] * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Linear form: weak gradient of scalar source term
// -----------------------------------------------------------------------------
void feSysElm_GradSource::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _gradPhi.resize(_nComponents * _nFunctions);
}

void feSysElm_GradSource::computeBe(feBilinearForm *form)
{
  double jac, S, div_v;

  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    S = (*_source)(form->_tn, _pos);

    // Compute compacted grad(phi) without the trivial zeros
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhi.data());

    for(int i = 0; i < _nFunctions; ++i) {
      // Divergence of test function:
      div_v = _gradPhi[i * _nComponents + (i % _nComponents)];
      form->_Be[i] += S * div_v * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: mass matrix for vector-valued field
// -----------------------------------------------------------------------------
void feSysElm_VectorMass::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _u.resize(_nComponents);
  _phiU.resize(_nFunctions);
}

void feSysElm_VectorMass::computeAe(feBilinearForm *form)
{
  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    for(int i = 0; i < _nFunctions; ++i)
      _phiU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);

    // TODO: CHECK THE DOT PRODUCT, CFR MIXEDVECTORMASS BELOW
    for(int i = 0; i < _nFunctions; ++i)
      for(int j = 0; j < _nFunctions; ++j)
        form->_Ae[i][j] += coeff * _phiU[i] * _phiU[j] * jac * _wQuad[k];
  }
}

void feSysElm_VectorMass::computeBe(feBilinearForm *form)
{
  double jac, coeff;

  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);

    for(int i = 0; i < _nFunctions; ++i) {
      _phiU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= coeff * _u[i % _nComponents] * _phiU[i] * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: mixed mass matrix for two vector-valued fields
// -----------------------------------------------------------------------------
void feSysElm_MixedVectorMass::createElementarySystem(std::vector<feSpace *> &space)
{
  _idV = 0;
  _idU = 1;
  _fieldsLayoutI = {_idV}; // Rectangular local matrix
  _fieldsLayoutJ = {_idU};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsV = space[_idV]->getNumFunctions();
  _nComponentsU = space[_idU]->getNumComponents();
  _nComponentsV = space[_idV]->getNumComponents();
  _u.resize(_nComponentsU);
  _phiU.resize(_nFunctionsU);
  _phiV.resize(_nFunctionsV);
}

void feSysElm_MixedVectorMass::computeAe(feBilinearForm *form)
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

    for(int i = 0; i < _nFunctionsV; ++i) {
      for(int j = 0; j < _nFunctionsU; ++j) {
        // The dot product du * v is nonzero only when the components match
        if(i % _nComponentsU != j % _nComponentsV) continue;
        form->_Ae[i][j] += coeff * _phiV[i] * _phiU[j] * jac * _wQuad[k];
      }
    }
  }
}

void feSysElm_MixedVectorMass::computeBe(feBilinearForm *form)
{
  double jac, coeff;

  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u,
                                                             _nComponentsU);

    for(int i = 0; i < _nFunctionsV; ++i) {
      _phiV[i] = form->_intSpaces[_idV]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= coeff * _phiV[i] * _u[i % _nComponentsU] * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: mixed scalar-vector mass matrix
// -----------------------------------------------------------------------------
void feSysElm_MixedScalarVectorMass::createElementarySystem(std::vector<feSpace *> &space)
{
  _idV = 0;
  _idU = 1;
  _fieldsLayoutI = {_idV};
  _fieldsLayoutJ = {_idU};
  _nFunctionsV = space[_idV]->getNumFunctions();
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nComponentsU = space[_idU]->getNumComponents();
  _phiV.resize(_nFunctionsV);
  _phiU.resize(_nFunctionsU);
  _u.resize(_nComponentsU);
}

void feSysElm_MixedScalarVectorMass::computeAe(feBilinearForm *form)
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

void feSysElm_MixedScalarVectorMass::computeBe(feBilinearForm *form)
{
  double jac, coeff;

  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u,
                                                             _nComponentsU);

    for(int i = 0; i < _nFunctionsV; ++i) {
      _phiV[i] = form->_intSpaces[_idV]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= coeff * _phiV[i] * _u[i % _nComponentsU] * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: transient mass matrix for vector-valued field
// -----------------------------------------------------------------------------
void feSysElm_TransientVectorMass::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _phiU.resize(_nFunctions);
  _dudt.resize(_nComponents);
}

void feSysElm_TransientVectorMass::computeAe(feBilinearForm *form)
{
  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    for(int i = 0; i < _nFunctions; ++i)
      _phiU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);

    for(int i = 0; i < _nFunctions; ++i) {
      for(int j = 0; j < _nFunctions; ++j) {
        // The dot product dudt * v is nonzero only when the components match
        if(i % _nComponents != j % _nComponents) continue;
        form->_Ae[i][j] += _phiU[i] * coeff * form->_c0 * _phiU[j] * jac * _wQuad[k];
      }
    }
  }
}

void feSysElm_TransientVectorMass::computeBe(feBilinearForm *form)
{
  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_solDot[_idU], k, _dudt,
                                                             _nComponents);

    for(int i = 0; i < _nFunctions; ++i) {
      _phiU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= coeff * _dudt[i % _nComponents] * _phiU[i] * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: diffusion of vector-valued field
// -----------------------------------------------------------------------------
void feSysElm_VectorDiffusion::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _gradu.resize(_nComponents * _nComponents);
  _gradPhi.resize(_nComponents * _nFunctions);
}

void feSysElm_VectorDiffusion::computeAe(feBilinearForm *form)
{
  double jac, doubleContraction, coeff, diffusivity;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar diffusivity k(t,x)
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    diffusivity = (*_diffusivity)(form->_tn, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    // Compute compacted grad(phi) without the trivial zeros
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhi.data());

    // Increment A with -diff * grad(u) : grad(v)
    for(int i = 0; i < _nFunctions; ++i) {
      for(int j = 0; j < _nFunctions; ++j) {
        doubleContraction = 0.;
        // Only contractions of the form
        // (phi_i, 0) : (phi_j, 0) or (0, phi_i) : (0, phi_j)
        // are nonzero (mixed products are zero)
        if((i % _nComponents) == (j % _nComponents)) {
          int iOffset = i * _nComponents;
          int jOffset = j * _nComponents;
          // Only "nComponents" terms in the sum are nonzero, since the vector test
          // functions are v_i = (phi_i, 0) or (0, phi_i)
          for(int m = 0; m < _nComponents; ++m) {
            doubleContraction += _gradPhi[iOffset + m] * _gradPhi[jOffset + m];
          }
          form->_Ae[i][j] += coeff * diffusivity * doubleContraction * jac * _wQuad[k];
        }
      }
    }
  }
}

void feSysElm_VectorDiffusion::computeBe(feBilinearForm *form)
{
  double doubleContraction;
  for(int k = 0; k < _nQuad; ++k) {
    double jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar diffusivity k(t,x)
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    double diffusivity = (*_diffusivity)(form->_tn, _pos);
    double coeff = (*_coeff)(form->_tn, _pos);

    // Compute gradient of current solution:
    // gradu = [grad(u_1) ... grad(u_nComponents)]
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());

    // Compute compacted grad(phi) without the trivial zeros
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhi.data());

    for(int i = 0; i < _nFunctions; ++i) {
      doubleContraction = 0.;
      int offset = i * _nComponents;
      // Only "nComponents" terms in the sum are nonzero, since the vector test
      // functions are v_i = (phi_i, 0) or (0, phi_i)
      for(int m = 0; m < _nComponents; ++m) {
        doubleContraction += _gradu[(i % _nComponents) * _nComponents + m] * _gradPhi[offset + m];
      }
      form->_Be[i] -= coeff * diffusivity * doubleContraction * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: mixed weak gradient
// -----------------------------------------------------------------------------
void feSysElm_MixedGradient::createElementarySystem(std::vector<feSpace *> &space)
{
  _idV = 0; // Vector-valued test functions
  _idU = 1; // Scalar unknown field
  _fieldsLayoutI = {_idV}; // Rectangular local matrix
  _fieldsLayoutJ = {_idU};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsV = space[_idV]->getNumFunctions();
  _nComponents = space[_idV]->getNumComponents();
  _phiU.resize(_nFunctionsU);
  _gradPhiV.resize(_nComponents * _nFunctionsV);
}

void feSysElm_MixedGradient::computeAe(feBilinearForm *form)
{
  double jac, coeff, div_v;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    // Get phiU
    form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);

    // Compute compacted grad(phiV) without the trivial zeros
    form->_intSpaces[_idV]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiV.data());

    for(int i = 0; i < _nFunctionsV; ++i) {
      for(int j = 0; j < _nFunctionsU; ++j) {
        div_v = _gradPhiV[i * _nComponents + (i % _nComponents)];
        form->_Ae[i][j] -= coeff * _phiU[j] * div_v * jac * _wQuad[k];
      }
    }
  }
}

void feSysElm_MixedGradient::computeBe(feBilinearForm *form)
{
  double jac, coeff, u, div_v;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    // Get current scalar solution
    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);

    // Compute compacted grad(phiV) without the trivial zeros
    form->_intSpaces[_idV]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiV.data());

    for(int i = 0; i < _nFunctionsV; ++i) {
      div_v = _gradPhiV[i * _nComponents + (i % _nComponents)];
      form->_Be[i] += coeff * u * div_v * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: mixed divergence
// -----------------------------------------------------------------------------
void feSysElm_MixedDivergence::createElementarySystem(std::vector<feSpace *> &space)
{
  _idV = 0; // Scalar test functions
  _idU = 1; // Vector-valued unknown field
  _fieldsLayoutI = {_idV}; // Rectangular local matrix
  _fieldsLayoutJ = {_idU};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _nFunctionsV = space[_idV]->getNumFunctions();
  _gradu.resize(_nComponents * _nComponents);
  _gradPhiU.resize(_nComponents * _nFunctionsU);
  _phiV.resize(_nFunctionsV);
}

void feSysElm_MixedDivergence::computeAe(feBilinearForm *form)
{
  double jac, coeff, div_phiU;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    // Get _phiV
    form->_intSpaces[_idV]->getFunctionsAtQuadNode(k, _phiV);

    // Get _gradPhiU
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiU.data());

    for(int i = 0; i < _nFunctionsV; ++i) {
      for(int j = 0; j < _nFunctionsU; ++j) {
        div_phiU = _gradPhiU[j * _nComponents + (j % _nComponents)];
        form->_Ae[i][j] += coeff * _phiV[i] * div_phiU * jac * _wQuad[k];
      }
    }
  }
}

void feSysElm_MixedDivergence::computeBe(feBilinearForm *form)
{
  double jac, coeff, div_u;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    // Compute gradient of current solution:
    // gradu = [grad(u_1) ... grad(u_nComponents)]
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());

    div_u = 0.;
    for(int iComp = 0; iComp < _nComponents; ++iComp) {
      div_u += _gradu[(iComp + 1) * (iComp + 1) - 1];
    }

    // Get phiV
    form->_intSpaces[_idV]->getFunctionsAtQuadNode(k, _phiV);

    for(int i = 0; i < _nFunctionsV; ++i) {
      form->_Be[i] -= coeff * div_u * _phiV[i] * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: mixed curl
// -----------------------------------------------------------------------------
void feSysElm_MixedCurl::createElementarySystem(std::vector<feSpace *> &space)
{
  _idV = 0; // Scalar test functions
  _idU = 1; // Vector-valued unknown field
  _fieldsLayoutI = {_idV}; // Rectangular local matrix
  _fieldsLayoutJ = {_idU};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nComponents = 2;
  _nFunctionsV = space[_idV]->getNumFunctions();
  _gradu.resize(_nComponents * _nComponents);
  _phiV.resize(_nFunctionsV);
}

void feSysElm_MixedCurl::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  // ...
}

void feSysElm_MixedCurl::computeBe(feBilinearForm *form)
{
  double jac, coeff, curl_u;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    // Compute phi and gradu
    form->_intSpaces[_idV]->getFunctionsAtQuadNode(k, _phiV);
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());

    // Only for 2D fields
    curl_u = _gradu[2] - _gradu[1];

    for(int i = 0; i < _nFunctionsV; ++i) {
      form->_Be[i] -= coeff * curl_u * _phiV[i] * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: mixed weak gradient involving 2 scalar fields and the test
//                functions of a 3rd vector field
// -----------------------------------------------------------------------------
void feSysElm_TripleMixedGradient::createElementarySystem(std::vector<feSpace *> &space)
{
  _idW = 0; // Vector-valued test functions
  _idU = 1; // Scalar unknown field
  _idV = 2; // Scalar unknown field
  _fieldsLayoutI = {_idW}; // Rectangular local matrix
  _fieldsLayoutJ = {_idU, _idV};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsV = space[_idV]->getNumFunctions();
  _nFunctionsW = space[_idW]->getNumFunctions();
  _nComponents = space[_idW]->getNumComponents();
  _gradV.resize(_nComponents); // Space dimension = nComponents of vector field
  _phiW.resize(_nFunctionsW);
}

void feSysElm_TripleMixedGradient::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Implement matrix for TripleMixedGradient");
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

void feSysElm_TripleMixedGradient::computeBe(feBilinearForm *form)
{
  double jac, coeff, u, gradVDotPhiW;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    // Get u, gradV and phiW
    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    form->_intSpaces[_idV]->interpolateFieldAtQuadNode_physicalGradient(form->_sol[_idV], k, form->_transformation, _gradV.data());
    form->_intSpaces[_idW]->getFunctionsAtQuadNode(k, _phiW);

    for(int i = 0; i < _nFunctionsW; ++i)
    {
      // _phiW contains the basis functions for each component of the vector field.
      // Take the i-th vector of _phiW = [phiW_x,1 phiW_y,1 phiW_z,1 ... phiW_x,N phiW_y,N phiW_z,N]
      gradVDotPhiW = dotProdN(_nComponents, _gradV.data(), _phiW.data() + i*_nComponents);
      form->_Be[i] -= coeff * u * gradVDotPhiW * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: CHNS Momentum equation
// -----------------------------------------------------------------------------
void feSysElm_CHNS_Momentum::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU   = 0;
  _idP   = 1;
  _idPhi = 2;
  _idMu  = 3;
  _fieldsLayoutI = {_idU};
  _fieldsLayoutJ = {_idU, _idP, _idPhi, _idMu};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _f.resize(_nComponents);
  _u.resize(_nComponents);
  _dudt.resize(_nComponents);
  _gradu.resize(_nComponents * _nComponents);
  _uDotGradu.resize(_nComponents);
  _gradphi.resize(_nComponents);
  _gradmu.resize(_nComponents);
  _gradmuDotGradu.resize(_nComponents);
  _symmetricGradu.resize(_nComponents * _nComponents);
  // phiU is the test function of U, unrelated to the phase marker Phi
  _phiU.resize(_nFunctionsU);
  _gradPhiU.resize(_nComponents * _nFunctionsU);
}

void feSysElm_CHNS_Momentum::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Implement FE matrix for feSysElm_CHNS_Momentum weak form.");
  // Computed with finite differences for now
}

void feSysElm_CHNS_Momentum::computeBe(feBilinearForm *form)
{
  double jac, rho, M, drhodphi, coeffKorteweg, eta, p, phi, mu, div_v, doubleContraction;
  double dudtDotphiU, uDotGraduDotphiU, gradPhiDotphiU, fDotphiU, gradMuDotgradUdotphiU;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficients
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    phi = form->_intSpaces[_idPhi]->interpolateFieldAtQuadNode(form->_sol[_idPhi], k);
    rho = (*_density)(phi);
    M = (*_mobility)(phi);
    drhodphi = (*_drhodphi)(phi);
    eta = (*_viscosity)(phi);
    coeffKorteweg = (*_coeffKorteweg)(form->_tn, _pos);
    (*_volumeForce)(phi, _f);

    // Get all relevant values and gradients
    // dudt
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_solDot[_idU], k, _dudt, _nComponents);
    // u
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    // grad_u
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());
    // p
    p = form->_intSpaces[_idP]->interpolateFieldAtQuadNode(form->_sol[_idP], k);
    // grad_phi
    form->_intSpaces[_idPhi]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idPhi], k, form->_transformation, _gradphi.data());
    // mu
    mu = form->_intSpaces[_idMu]->interpolateFieldAtQuadNode(form->_sol[_idMu], k);
    // grad_mu
    form->_intSpaces[_idMu]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idMu], k, form->_transformation, _gradmu.data());

    // phiU
    form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);
    // grad_phiU
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiU.data());
    // Compute u dot gradu
    std::fill(_uDotGradu.begin(), _uDotGradu.end(), 0.);
    for(int m = 0; m < _nComponents; ++m) {
      for(int n = 0; n < _nComponents; ++n) {
        _uDotGradu[m] += _u[n] * _gradu[m * _nComponents + n];
      }
    }

    // Compute gradmu dot gradu
    std::fill(_gradmuDotGradu.begin(), _gradmuDotGradu.end(), 0.);
    for(int m = 0; m < _nComponents; ++m) {
      for(int n = 0; n < _nComponents; ++n) {
        _gradmuDotGradu[m] += _gradmu[n] * _gradu[m * _nComponents + n];
      }
    }

    // Compute symmetric gradient (2*d(u))
    for(int m = 0; m < _nComponents; ++m) {
      for(int n = 0; n < _nComponents; ++n) {
        _symmetricGradu[m * _nComponents + n] =
          _gradu[m * _nComponents + n] + _gradu[n * _nComponents + m];
      }
    }

    for(int i = 0; i < _nFunctionsU; ++i)
    {
      // dudt dot phiU (always only 1 component to the dot product as vector test functions have 1 nonzero component)
      dudtDotphiU = _dudt[i % _nComponents] * _phiU[i];

      // (u dot gradu) dot phiU
      uDotGraduDotphiU = _uDotGradu[i % _nComponents] * _phiU[i];

      // grad_phi dot phiU
      gradPhiDotphiU = _gradphi[i % _nComponents] * _phiU[i];

      // (grad_mu dot grad_u) dot phiU
      gradMuDotgradUdotphiU = _gradmuDotGradu[i % _nComponents] * _phiU[i];

      // Divergence of test function
      div_v = _gradPhiU[i * _nComponents + (i % _nComponents)];

      // Double contraction 2*d(u) : grad(phiU) = (gradu + gradu^T) : grad(phiU)
      doubleContraction = 0.;
      for(int m = 0; m < _nComponents; ++m) {
        doubleContraction +=
          _symmetricGradu[(i % _nComponents) * _nComponents + m] * _gradPhiU[i * _nComponents + m];
      }

      // Source term
      fDotphiU = _f[i % _nComponents] * _phiU[i];

      // Increment RHS
      form->_Be[i] -= (rho * (dudtDotphiU + uDotGraduDotphiU - fDotphiU)
                     + M * drhodphi * gradMuDotgradUdotphiU
                     + coeffKorteweg * mu * gradPhiDotphiU
                     - p * div_v
                     + eta * doubleContraction) * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: convective acceleration of vector-valued field
// -----------------------------------------------------------------------------
void feSysElm_VectorConvectiveAcceleration::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _u.resize(_nComponents);
  _uDotGradu.resize(_nComponents);
  _gradu.resize(_nComponents * _nComponents);
  _phiU.resize(_nFunctions);
  _gradPhiU.resize(_nComponents * _nFunctions);
}

void feSysElm_VectorConvectiveAcceleration::computeAe(feBilinearForm *form)
{
  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    // Get u, grad_u, phiU and gradPhiU
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());
    form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiU.data());

    for(int i = 0; i < _nFunctions; ++i) {
      for(int j = 0; j < _nFunctions; ++j) {
        // Compute (u0 dot gradu) + (u dot gradu0)

        // Generic but slow
        // std::fill(_uDotGradu.begin(), _uDotGradu.end(), 0.);

        // for(int n = 0; n < _nComponents; ++n) {
        //   _uDotGradu[j % _nComponents] += _u[n] * _gradPhiU[j * _nComponents + n];
        // }
        // for(int n = 0; n < _nComponents; ++n) {
        //   _uDotGradu[n] += _phiU[j] * _gradu[n * _nComponents + (j % _nComponents)];
        // }

        // Explicit and faster computation for 2D velocity field
        if(j % _nComponents == 0) {
          _uDotGradu[0] = _u[0] * _gradPhiU[j * _nComponents] +
                          _u[1] * _gradPhiU[j * _nComponents + 1] +
                          _phiU[j] * _gradu[0 * _nComponents + (j % _nComponents)];
          _uDotGradu[1] = _phiU[j] * _gradu[1 * _nComponents + (j % _nComponents)];
        } else {
          _uDotGradu[0] = _phiU[j] * _gradu[0 * _nComponents + (j % _nComponents)];
          _uDotGradu[1] = _u[0] * _gradPhiU[j * _nComponents] +
                          _u[1] * _gradPhiU[j * _nComponents + 1] +
                          _phiU[j] * _gradu[1 * _nComponents + (j % _nComponents)];
        }

        // Only one nonzero component to both v_i and udotgradu, hence dot product
        // ((u0 dot gradu) + (u dot gradu0)) dot v is only one term
        form->_Ae[i][j] += coeff * _phiU[i] * _uDotGradu[i % _nComponents] * jac * _wQuad[k];
      }
    }
  }
}

void feSysElm_VectorConvectiveAcceleration::computeBe(feBilinearForm *form)
{
  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    // Get u, grad_u and phiU
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());
    form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);

    // Compute (u dot gradu)
    std::fill(_uDotGradu.begin(), _uDotGradu.end(), 0.);
    for(int m = 0; m < _nComponents; ++m) {
      for(int n = 0; n < _nComponents; ++n) {
        _uDotGradu[m] += _u[n] * _gradu[m * _nComponents + n];
      }
    }

    // Increment with (u dot gradu) dot v
    for(int i = 0; i < _nFunctions; ++i) {
      form->_Be[i] -= coeff * _phiU[i] * _uDotGradu[i % _nComponents] * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------------
// Bilinear form: convection of a scalar tracer C in the resolved velocity field u
// -----------------------------------------------------------------------------------
void feSysElm_TracerConvection::createElementarySystem(std::vector<feSpace *> &space)
{
  _idC = 0;
  _idU = 1;
  _fieldsLayoutI = {_idC};
  _fieldsLayoutJ = {_idU};
  _nFunctionsC = space[_idC]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _u.resize(_nComponents);
  _gradC.resize(_nComponents); // Dimension = number of velocity components
  _phiC.resize(_nFunctionsC);
}

void feSysElm_TracerConvection::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Implement FE matrix for TRACER_CONVECTION weak form.");
  // double jac, coeff;
  // for(int k = 0; k < _nQuad; ++k) {
  //   jac = form->_J[_nQuad * form->_numElem + k];
  //   form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

  //   // Evaluate scalar coefficient
  //   form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
  //   coeff = (*_coeff)(form->_tn, _pos);

  //   // Get u, grad_u, phiU and gradPhiU
  //   form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
  //   form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(
  //     form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());
  //   form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);
  //   form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
  //                                                                  _gradPhiU.data());

  //   for(int i = 0; i < _nFunctions; ++i) {
  //     for(int j = 0; j < _nFunctions; ++j) {
  //       // Compute (u0 dot gradu) + (u dot gradu0)

  //       // Generic but slow
  //       // std::fill(_uDotGradu.begin(), _uDotGradu.end(), 0.);

  //       // for(int n = 0; n < _nComponents; ++n) {
  //       //   _uDotGradu[j % _nComponents] += _u[n] * _gradPhiU[j * _nComponents + n];
  //       // }
  //       // for(int n = 0; n < _nComponents; ++n) {
  //       //   _uDotGradu[n] += _phiU[j] * _gradu[n * _nComponents + (j % _nComponents)];
  //       // }

  //       // Explicit and faster computation for 2D velocity field
  //       if(j % _nComponents == 0) {
  //         _uDotGradu[0] = _u[0] * _gradPhiU[j * _nComponents] +
  //                         _u[1] * _gradPhiU[j * _nComponents + 1] +
  //                         _phiU[j] * _gradu[0 * _nComponents + (j % _nComponents)];
  //         _uDotGradu[1] = _phiU[j] * _gradu[1 * _nComponents + (j % _nComponents)];
  //       } else {
  //         _uDotGradu[0] = _phiU[j] * _gradu[0 * _nComponents + (j % _nComponents)];
  //         _uDotGradu[1] = _u[0] * _gradPhiU[j * _nComponents] +
  //                         _u[1] * _gradPhiU[j * _nComponents + 1] +
  //                         _phiU[j] * _gradu[1 * _nComponents + (j % _nComponents)];
  //       }

  //       // Only one nonzero component to both v_i and udotgradu, hence dot product
  //       // ((u0 dot gradu) + (u dot gradu0)) dot v is only one term
  //       form->_Ae[i][j] += coeff * _phiU[i] * _uDotGradu[i % _nComponents] * jac * _wQuad[k];
  //     }
  //   }
  // }
}

void feSysElm_TracerConvection::computeBe(feBilinearForm *form)
{
  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    // Get u, gradC and phiC
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    form->_intSpaces[_idC]->interpolateFieldAtQuadNode_physicalGradient(form->_sol[_idC], k, form->_transformation, _gradC.data());
    form->_intSpaces[_idC]->getFunctionsAtQuadNode(k, _phiC);

    // Compute u dot gradC
    double uDotGradC = 0.;
    for(int m = 0; m < _nComponents; ++m) {
      uDotGradC += _u[m] * _gradC[m];
    }

    // Increment elementwise residual with (u dot gradC) * phiC
    for(int i = 0; i < _nFunctionsC; ++i) {
      form->_Be[i] -= coeff * uDotGradC * _phiC[i] * jac * _wQuad[k];
    }
  }
}

// ------------------------------------------------------------------------------------------------------
// Bilinear form: convective acceleration of vector-valued field for the adjoint Navier-Stokes equations
// ------------------------------------------------------------------------------------------------------
void feSysElm_VectorAdjointConvectiveAcceleration::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idUad = 1;
  _fieldsLayoutI = {_idUad}; // Local matrix is only the submatrix at (uAd, phi_uAD) ?
  _fieldsLayoutJ = {_idUad};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsUad = space[_idUad]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _u.resize(_nFunctionsU);
  _uDotGraduAd.resize(_nComponents);
  _grad_uAd.resize(_nComponents * _nComponents);
  _phi_uAd.resize(_nFunctionsUad);
  _gradPhiU.resize(_nFunctionsU * _nFunctionsU);
}

void feSysElm_VectorAdjointConvectiveAcceleration::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Implement FE matrix for VectorAdjointConvectiveAcceleration weak form.");
  // double jac, coeff;
  // for(int k = 0; k < _nQuad; ++k) {
  //   jac = form->_J[_nQuad * form->_numElem + k];
  //   form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

  //   // Evaluate scalar coefficient
  //   form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
  //   coeff = (*_coeff)(form->_tn, _pos);

  //   // Get u, grad_u, phiU and gradPhiU
  //   form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
  //   form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(
  //     form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());
  //   form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);
  //   form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
  //                                                                  _gradPhiU.data());

  //   for(int i = 0; i < _nFunctions; ++i) {
  //     for(int j = 0; j < _nFunctions; ++j) {
  //       // Compute (u0 dot gradu) + (u dot gradu0)

  //       // Generic but slow
  //       // std::fill(_uDotGradu.begin(), _uDotGradu.end(), 0.);

  //       // for(int n = 0; n < _nComponents; ++n) {
  //       //   _uDotGradu[j % _nComponents] += _u[n] * _gradPhiU[j * _nComponents + n];
  //       // }
  //       // for(int n = 0; n < _nComponents; ++n) {
  //       //   _uDotGradu[n] += _phiU[j] * _gradu[n * _nComponents + (j % _nComponents)];
  //       // }

  //       // Explicit and faster computation for 2D velocity field
  //       if(j % _nComponents == 0) {
  //         _uDotGradu[0] = _u[0] * _gradPhiU[j * _nComponents] +
  //                         _u[1] * _gradPhiU[j * _nComponents + 1] +
  //                         _phiU[j] * _gradu[0 * _nComponents + (j % _nComponents)];
  //         _uDotGradu[1] = _phiU[j] * _gradu[1 * _nComponents + (j % _nComponents)];
  //       } else {
  //         _uDotGradu[0] = _phiU[j] * _gradu[0 * _nComponents + (j % _nComponents)];
  //         _uDotGradu[1] = _u[0] * _gradPhiU[j * _nComponents] +
  //                         _u[1] * _gradPhiU[j * _nComponents + 1] +
  //                         _phiU[j] * _gradu[1 * _nComponents + (j % _nComponents)];
  //       }

  //       // Only one nonzero component to both v_i and udotgradu, hence dot product
  //       // ((u0 dot gradu) + (u dot gradu0)) dot v is only one term
  //       form->_Ae[i][j] += coeff * _phiU[i] * _uDotGradu[i % _nComponents] * jac * _wQuad[k];
  //     }
  //   }
  // }
}

void feSysElm_VectorAdjointConvectiveAcceleration::computeBe(feBilinearForm *form)
{
  double jac, coeff;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    // Get u, grad_uAd and phi_uAd
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    form->_intSpaces[_idUad]->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idUad], _nComponents, k, form->_transformation, _grad_uAd.data());
    form->_intSpaces[_idUad]->getFunctionsAtQuadNode(k, _phi_uAd);

    // Compute (u dot grad_uAd)
    std::fill(_uDotGraduAd.begin(), _uDotGraduAd.end(), 0.);
    for(int m = 0; m < _nComponents; ++m) {
      for(int n = 0; n < _nComponents; ++n) {
        _uDotGraduAd[m] += _u[n] * _grad_uAd[m * _nComponents + n];
        // ESSAI : grad transpose
        _uDotGraduAd[m] += _u[n] * _grad_uAd[n * _nComponents + m];
      }
    }

    // Increment with (u dot grad_uAd) dot v
    for(int i = 0; i < _nFunctionsUad; ++i) {
      form->_Be[i] -= coeff * _phi_uAd[i] * _uDotGraduAd[i % _nComponents] * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear form: weak divergence of newtonian stress tensor (mixed u and p)
// -----------------------------------------------------------------------------
void feSysElm_DivergenceNewtonianStress::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idP = 1;
  _fieldsLayoutI = {_idU};
  _fieldsLayoutJ = {_idU, _idP};
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _gradu.resize(_nComponents * _nComponents);
  _symmetricGradu.resize(_nComponents * _nComponents);
  _gradPhiU.resize(_nComponents * _nFunctionsU);
}

void feSysElm_DivergenceNewtonianStress::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Implement FE matrix for feSysElm_DivergenceNewtonianStress weak form.");
  // Computed with finite differences for now
}

void feSysElm_DivergenceNewtonianStress::computeBe(feBilinearForm *form)
{
  double jac, coeff, mu, p, div_v, doubleContraction;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficients
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    mu = (*_viscosity)(form->_tn, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    // Get p, grad_u and gradphiU
    p = form->_intSpaces[_idP]->interpolateFieldAtQuadNode(form->_sol[_idP], k);
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiU.data());

    for(int m = 0; m < _nComponents; ++m) {
      for(int n = 0; n < _nComponents; ++n) {
        _symmetricGradu[m * _nComponents + n] =
          _gradu[m * _nComponents + n] + _gradu[n * _nComponents + m];
      }
    }

    for(int i = 0; i < _nFunctionsU; ++i) {
      // Divergence of test function
      div_v = _gradPhiU[i * _nComponents + (i % _nComponents)];
      // Double contraction d(u) : grad(v)
      doubleContraction = 0.;
      for(int m = 0; m < _nComponents; ++m) {
        doubleContraction +=
          _symmetricGradu[(i % _nComponents) * _nComponents + m] * _gradPhiU[i * _nComponents + m];
      }

      form->_Be[i] -= coeff * (-p * div_v + mu * doubleContraction) * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Bilinear forms: Stabilization for the (Navier-)Stokes equations
// -----------------------------------------------------------------------------
// static double hFortin(int dim, const double velocity[3], double normVelocity, int nFunctions,
//                       const std::vector<double> &gradphi)
// {
//   // Compute h from (13.17) in Fortin & Garon
//   double res = 0., dotprod;
//   for(int i = 0; i < nFunctions; ++i) {
//     dotprod = 0.;
//     for(int iDim = 0; iDim < dim; ++iDim) {
//       dotprod += velocity[iDim] / normVelocity * gradphi[i * dim + iDim];
//     }
//     res += dotprod * dotprod;
//   }

//   return sqrt(2.) / sqrt(res);
// }

// static double getInnerRadius(const std::vector<double> &triCoord)
// {
//   // radius of inscribed circle = 2 * Area / sum(Line_i)
//   double dist[3], k = 0.;
//   for(int i = 0; i < 3; i++) {
//     double x0 = triCoord[3 * i + 0];
//     double y0 = triCoord[3 * i + 1];
//     double x1 = triCoord[3 * ((i + 1) % 3) + 0];
//     double y1 = triCoord[3 * ((i + 1) % 3) + 1];
//     dist[i] = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
//     k += 0.5 * dist[i];
//   }
//   double const area = std::sqrt(k * (k - dist[0]) * (k - dist[1]) * (k - dist[2]));
//   return area / k;
// }

// static double getOuterRadius(const std::vector<double> &triCoord)
// {
//   // radius of circle circumscribing a triangle
//   double dist[3], k = 0.0;
//   for(int i = 0; i < 3; i++) {
//     double x0 = triCoord[3 * i + 0];
//     double y0 = triCoord[3 * i + 1];
//     double x1 = triCoord[3 * ((i + 1) % 3) + 0];
//     double y1 = triCoord[3 * ((i + 1) % 3) + 1];
//     dist[i] = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
//     k += 0.5 * dist[i];
//   }
//   double const area = std::sqrt(k * (k - dist[0]) * (k - dist[1]) * (k - dist[2]));
//   return dist[0] * dist[1] * dist[2] / (4 * area);
// }

// static double hCircumInscr(const std::vector<double> &triCoord)
// {
//   double hInscribed = 2. * getInnerRadius(triCoord);
//   double hCircumscribed = 2. * getOuterRadius(triCoord);
//   return (hInscribed + hCircumscribed) / 2.;
// }

// Diameter of circle of equivalent area
inline double hEquivalentCircle(const double jac) { return 2. * sqrt(jac / (2. * M_PI)); }

static double compute_tau(const std::vector<double> &triCoord, int dim, double velocity[3],
                          double viscosity, int nFunctions, std::vector<double> &gradphi,
                          double jac)
{
  // Unused parameters unless using hFortin above
  UNUSED(triCoord, nFunctions, gradphi);

  double normV = 0.;
  for(int i = 0; i < dim; ++i) {
    normV += velocity[i] * velocity[i];
  }
  normV = sqrt(normV);

  // double h1 = hFortin(dim, velocity, normV, nFunctions, gradphi);
  if(normV < 1e-12) {
    return 0.;
  }
  // double h2 = hCircumInscr(triCoord);
  double h3 = hEquivalentCircle(jac);

  double h = h3;

  double r1 = 2. * normV / h;
  double r2 = 4. * viscosity / (h * h);

  return 1. / sqrt(r1 * r1 + 9. * r2 * r2);
}

void feSysElm_Stokes_SUPG_PSPG::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idP = 1;
  _fieldsLayoutI = {_idU, _idP};
  _fieldsLayoutJ = {_idU, _idP};
  _nComponents = space[_idU]->getNumComponents();
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsP = space[_idP]->getNumFunctions();
  _f.resize(_nComponents);
  _residual.resize(_nComponents);
  _u.resize(_nComponents);
  _gradu.resize(_nComponents * _nComponents);
  _gradp.resize(_nComponents);
  _gradPhiU.resize(_nComponents * _nFunctionsU);
  _gradPhiP.resize(_nComponents * _nFunctionsP);

  int numUniqueMixedDerivatives[3] = {1, 3, 6};
  _hessu.resize(numUniqueMixedDerivatives[_nComponents - 1] * _nComponents);
  _hessPhiU.resize(numUniqueMixedDerivatives[_nComponents - 1] * _nFunctionsU);
}

void feSysElm_Stokes_SUPG_PSPG::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  feErrorMsg(FE_STATUS_ERROR, "Implement FE matrix for feSysElm_Stokes_SUPG_PSPG weak form.");
  // Computed with finite differences
}

void feSysElm_Stokes_SUPG_PSPG::computeBe(feBilinearForm *form)
{
  double jac, coeff, rho, mu, tau, dotProdU, dotProdP, residualQ[2];

  for(int k = 0, cnt; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficients
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);
    rho = (*_density)(form->_tn, _pos);
    mu = (*_viscosity)(form->_tn, _pos);
    (*_volumeForce)(form->_tn, _pos, _f);

    // Get u, grad_u, grad_p, gradphiU, gradphiP and _hessu (2nd order derivatives)
    // u
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    // grad_u
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());
    // grad_p
    form->_intSpaces[_idP]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idP], k, form->_transformation, _gradp.data());
    // grad_phiU
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiU.data());
    // grad_phiP
    form->_intSpaces[_idP]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiP.data());
    // hess_u
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalHessian(
      form->_sol[_idU], _nComponents, k, form->_transformation, _hessu.data());
    // hessPhiU
    form->_intSpaces[_idU]->getFunctionsPhysicalHessianAtQuadNode(k, form->_transformation,
                                                                  _hessPhiU.data());

    // Stabilization parameter
    tau = compute_tau(form->_geoCoord, _nComponents, _u.data(), mu, _nFunctionsU, _gradPhiU, jac);

    // Residual for Navier-Stokes (laplacian form):
    double lapu[2] = {_hessu[0] + _hessu[2], _hessu[3] + _hessu[5]};
    _residual[0] = -_gradp[0] + mu * lapu[0] + rho * _f[0];
    _residual[1] = -_gradp[1] + mu * lapu[1] + rho * _f[1];

    cnt = 0;
    for(int iU = 0; iU < _nFunctionsU; ++iU) {
      double d2phidx2 = _hessPhiU[iU * 3 + 0];
      double d2phidy2 = _hessPhiU[iU * 3 + 2];
      double lap_phiU = d2phidx2 + d2phidy2;

      double residualV = mu * lap_phiU;

      // Dot product only has one non-zero term
      dotProdU = residualV * _residual[iU % _nComponents];

      form->_Be[cnt++] += coeff * tau * dotProdU * jac * _wQuad[k];
    }

    for(int iP = 0; iP < _nFunctionsP; ++iP) {
      residualQ[0] = -_gradPhiP[iP * _nComponents + 0];
      residualQ[1] = -_gradPhiP[iP * _nComponents + 1];

      dotProdP = _residual[0] * residualQ[0] + _residual[1] * residualQ[1];

      form->_Be[cnt++] -= coeff * tau * dotProdP * jac * _wQuad[k];
    }
  }
}

void feSysElm_NS_SUPG_PSPG::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idP = 1;
  _fieldsLayoutI = {_idU, _idP};
  _fieldsLayoutJ = {_idU, _idP};
  _nComponents = space[_idU]->getNumComponents();
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsP = space[_idP]->getNumFunctions();
  _f.resize(_nComponents);
  _residual.resize(_nComponents);
  _u.resize(_nComponents);
  _gradu.resize(_nComponents * _nComponents);
  _gradp.resize(_nComponents);
  _gradPhiU.resize(_nComponents * _nFunctionsU);
  _gradPhiP.resize(_nComponents * _nFunctionsP);
  _uDotGradu.resize(_nComponents);
  _uDotGradPhiu.resize(_nComponents);

  int numUniqueMixedDerivatives[3] = {1, 3, 6};
  _hessu.resize(numUniqueMixedDerivatives[_nComponents - 1] * _nComponents);
  _hessPhiU.resize(numUniqueMixedDerivatives[_nComponents - 1] * _nFunctionsU);
}

void feSysElm_NS_SUPG_PSPG::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  // Computed with finite differences
}

void feSysElm_NS_SUPG_PSPG::computeBe(feBilinearForm *form)
{
  double jac, coeff, rho, mu, tau, dotProdU, dotProdP, residualQ[2];

  for(int k = 0, cnt; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficients
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);
    rho = (*_density)(form->_tn, _pos);
    mu = (*_viscosity)(form->_tn, _pos);
    (*_volumeForce)(form->_tn, _pos, _f);

    // Get u, grad_u, grad_p, gradphiU, gradphiP and _hessu (2nd order derivatives)
    // u
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    // grad_u
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], _nComponents, k, form->_transformation, _gradu.data());
    // grad_p
    form->_intSpaces[_idP]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idP], k, form->_transformation, _gradp.data());
    // grad_phiU
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiU.data());
    // grad_phiP
    form->_intSpaces[_idP]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiP.data());
    // hess_u
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalHessian(
      form->_sol[_idU], _nComponents, k, form->_transformation, _hessu.data());
    // hessPhiU
    form->_intSpaces[_idU]->getFunctionsPhysicalHessianAtQuadNode(k, form->_transformation,
                                                                  _hessPhiU.data());

    // Stabilization parameter
    tau = compute_tau(form->_geoCoord, _nComponents, _u.data(), mu, _nFunctionsU, _gradPhiU, jac);

    // Compute (u dot gradu)
    std::fill(_uDotGradu.begin(), _uDotGradu.end(), 0.);
    for(int m = 0; m < _nComponents; ++m) {
      for(int n = 0; n < _nComponents; ++n) {
        _uDotGradu[m] += _u[n] * _gradu[m * _nComponents + n];
      }
    }

    // Residual for Navier-Stokes (laplacian form):
    double lapu[2] = {_hessu[0] + _hessu[2], _hessu[3] + _hessu[5]};
    _residual[0] = -_uDotGradu[0] - _gradp[0] + mu * lapu[0] + rho * _f[0];
    _residual[1] = -_uDotGradu[1] - _gradp[1] + mu * lapu[1] + rho * _f[1];

    cnt = 0;
    for(int iU = 0; iU < _nFunctionsU; ++iU) {
      // Compute (u dot grad_v)
      double uDotGradPhiU = 0.;
      for(int n = 0; n < _nComponents; ++n) {
        uDotGradPhiU += _u[n] * _gradPhiU[iU * _nComponents + n];
      }

      double d2phidx2 = _hessPhiU[iU * 3 + 0];
      double d2phidy2 = _hessPhiU[iU * 3 + 2];
      double lap_phiU = d2phidx2 + d2phidy2;

      double residualV = -uDotGradPhiU + mu * lap_phiU;

      // Dot product only has one non-zero term
      dotProdU = residualV * _residual[iU % _nComponents];

      form->_Be[cnt++] += coeff * tau * dotProdU * jac * _wQuad[k];
    }

    for(int iP = 0; iP < _nFunctionsP; ++iP) {
      residualQ[0] = -_gradPhiP[iP * _nComponents + 0];
      residualQ[1] = -_gradPhiP[iP * _nComponents + 1];

      dotProdP = _residual[0] * residualQ[0] + _residual[1] * residualQ[1];

      form->_Be[cnt++] -= coeff * tau * dotProdP * jac * _wQuad[k];
    }
  }
}

void feSysElm_GLS_Stokes::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idP = 1;
  _fieldsLayoutI = {_idU, _idP};
  _fieldsLayoutJ = {_idU, _idP};
  _nComponents = space[_idU]->getNumComponents();
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsP = space[_idP]->getNumFunctions();
  _f.resize(_nComponents);
  _residual.resize(_nComponents);
  _u.resize(_nComponents);
  _gradu.resize(_nComponents * _nComponents);
  _gradp.resize(_nComponents);
  _gradPhiU.resize(_nComponents * _nFunctionsU);
  _gradPhiP.resize(_nComponents * _nFunctionsP);

  int numUniqueMixedDerivatives[3] = {1, 3, 6};
  _hessu.resize(numUniqueMixedDerivatives[_nComponents - 1] * _nComponents);
  _hessPhiU.resize(numUniqueMixedDerivatives[_nComponents - 1] * _nFunctionsU);
}

void feSysElm_GLS_Stokes::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  // ...
}

void feSysElm_GLS_Stokes::computeBe(feBilinearForm *form)
{
  UNUSED(form);
  // ...
}

void feSysElm_GLS_NavierStokes::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _idP = 1;
  _fieldsLayoutI = {_idU, _idP};
  _fieldsLayoutJ = {_idU, _idP};
  _nComponents = space[_idU]->getNumComponents();
  _nFunctionsU = space[_idU]->getNumFunctions();
  _nFunctionsP = space[_idP]->getNumFunctions();
  _f.resize(_nComponents);
  _residual.resize(_nComponents);
  _u.resize(_nComponents);
  _gradu.resize(_nComponents * _nComponents);
  _gradp.resize(_nComponents);
  _gradPhiU.resize(_nComponents * _nFunctionsU);
  _gradPhiP.resize(_nComponents * _nFunctionsP);
  _uDotGradu.resize(_nComponents);
  _uDotGradPhiu.resize(_nComponents);

  int numUniqueMixedDerivatives[3] = {1, 3, 6};
  _hessu.resize(numUniqueMixedDerivatives[_nComponents - 1] * _nComponents);
  _hessPhiU.resize(numUniqueMixedDerivatives[_nComponents - 1] * _nFunctionsU);
}

void feSysElm_GLS_NavierStokes::computeAe(feBilinearForm *form)
{
  UNUSED(form);
  // ...
}

void feSysElm_GLS_NavierStokes::computeBe(feBilinearForm *form)
{
  UNUSED(form);
  // ...
}