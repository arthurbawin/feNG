#include "feSysElm.h"
#include "feBilinearForm.h"

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
      form->_Be[i] -= _phiU[i] * _S[ i % _nComponents ] * jac * _wQuad[k];
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
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhi.data());

    for(int i = 0; i < _nFunctions; ++i) {
    	// Divergence of test function:
    	div_v = _gradPhi[i*_nComponents + (i % _nComponents)];
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
  double jac, doubleContraction, diffusivity;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar diffusivity k(t,x)
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    diffusivity = (*_diffusivity)(form->_tn, _pos);

    // Compute compacted grad(phi) without the trivial zeros
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhi.data());

    // Increment A with -diff * grad(u) : grad(v)
    for(int i = 0; i < _nFunctions; ++i) {
    	for(int j = 0; j < _nFunctions; ++j) {
	    	doubleContraction = 0.;
	    	// Only contractions of the form
	    	// (phi_i, 0) : (phi_j, 0) or (0, phi_i) : (0, phi_j)
	    	// are nonzero (mixed products are zero)
	    	if( (i % _nComponents) == (j % _nComponents) ) {
			    int iOffset = i*_nComponents;
			    int jOffset = j*_nComponents;
			    // Only "nComponents" terms in the sum are nonzero, since the vector test
			    // functions are v_i = (phi_i, 0) or (0, phi_i)
			    for(int m = 0; m < _nComponents; ++m) {
			    	doubleContraction += _gradPhi[iOffset + m] * _gradPhi[jOffset + m];
				  }
		      form->_Ae[i][j] -= diffusivity * doubleContraction * jac * _wQuad[k];
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

    // Compute gradient of current solution:
    // gradu = [grad(u_1) ... grad(u_nComponents)]
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(form->_sol[_idU], _nComponents, 
      k, form->_transformation, _gradu.data());

    // Compute compacted grad(phi) without the trivial zeros
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhi.data());

    for(int i = 0; i < _nFunctions; ++i) {
    	doubleContraction = 0.;
	    int offset = i*_nComponents;
	    // Only "nComponents" terms in the sum are nonzero, since the vector test
	    // functions are v_i = (phi_i, 0) or (0, phi_i)
	    for(int m = 0; m < _nComponents; ++m) {
	    	doubleContraction += _gradu[(i % _nComponents) * _nComponents + m] * _gradPhi[offset + m];
		  }
      form->_Be[i] += diffusivity * doubleContraction * jac * _wQuad[k];
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
    form->_intSpaces[_idV]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiV.data());

    for(int i = 0; i < _nFunctionsV; ++i) {
      for(int j = 0; j < _nFunctionsU; ++j) {
        div_v = _gradPhiV[i*_nComponents + (i % _nComponents)];
        form->_Ae[i][j] -= coeff * _phiU[j] * div_v * jac * _wQuad[k];
        // feInfo("Adding %f - %f - %f - %f", 
        //   coeff * _phiU[j] * div_v * jac * _wQuad[k],
        //   coeff, _phiU[j], div_v);
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
    form->_intSpaces[_idV]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiV.data());

    for(int i = 0; i < _nFunctionsV; ++i) {
      div_v = _gradPhiV[i*_nComponents + (i % _nComponents)];
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
  // double jac, coeff, div_v;
  // for(int k = 0; k < _nQuad; ++k) {
  //   jac = form->_J[_nQuad * form->_numElem + k];

  //   // Evaluate scalar coefficient
  //   form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
  //   coeff = (*_coeff)(form->_tn, _pos);

  //   // Get phiU
  //   form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);

  //   // Compute compacted grad(phiV) without the trivial zeros
  //   form->_intSpaces[_idV]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, 
  //     form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, _gradPhiV.data());

  //   for(int i = 0; i < _nFunctionsV; ++i) {
  //     for(int j = 0; j < _nFunctionsU; ++j) {
  //       div_v = _gradPhiV[i*_nComponents + (i % _nComponents)];
  //       form->_Ae[i][j] -= coeff * _phiU[j] * div_v * jac * _wQuad[k];
  //       // feInfo("Adding %f - %f - %f - %f", 
  //       //   coeff * _phiU[j] * div_v * jac * _wQuad[k],
  //       //   coeff, _phiU[j], div_v);
  //     }
  //   }
  // }
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
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(form->_sol[_idU], _nComponents, 
      k,form->_transformation, _gradu.data());

    div_u = 0.;
    for(int iComp = 0; iComp < _nComponents; ++iComp){
      div_u += _gradu[(iComp + 1)*(iComp + 1) - 1];
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
  // TODO
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
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(form->_sol[_idU], _nComponents, 
      k, form->_transformation, _gradu.data());

    // Only for 2D fields
    curl_u = _gradu[2] - _gradu[1];

    for(int i = 0; i < _nFunctionsV; ++i) {
      form->_Be[i] -= coeff * curl_u * _phiV[i] * jac * _wQuad[k];
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
  // std::vector<double> test(_nComponents, 0.);
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    // Get u, grad_u, phiU and gradPhiU
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(form->_sol[_idU], _nComponents, 
      k, form->_transformation, _gradu.data());
    form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiU.data());

    for(int i = 0; i < _nFunctions; ++i) {
      for(int j = 0; j < _nFunctions; ++j) {

        // Compute (u0 dot gradu) + (u dot gradu0)
        std::fill(_uDotGradu.begin(), _uDotGradu.end(), 0.);
        // std::fill(test.begin(), test.end(), 0.);

        for(int n = 0; n < _nComponents; ++n) {
          _uDotGradu[j % _nComponents] += _u[n] * _gradPhiU[ j*_nComponents + n];
        }
        for(int n = 0; n < _nComponents; ++n) {
          _uDotGradu[n] += _phiU[j] * _gradu[ n*_nComponents + (j % _nComponents)];
        }

        // for(int n = 0; n < _nComponents; ++n) {
        //   test[j % _nComponents] += _u[n] * _gradPhiU[ j*_nComponents + n];
        //   test[n] += _phiU[j] * _gradu[ n*_nComponents + (j % _nComponents)];
        // }

        // for(int n = 0; n < _nComponents; ++n) {
        //   if(fabs(_uDotGradu[n] - test[n]) > 1e-14){
        //     feInfo("_uDotGradu[%d] = %1.18e", n, _uDotGradu[n]);
        //     feInfo("      test[%d] = %1.18e", n, test[n]);
        //     exit(-1);
        //   }
        // }

        // Only one nonzero component to both v_i and udotgradu, hence dot product
        // ((u0 dot gradu) + (u dot gradu0)) dot v is only one term
        form->_Ae[i][j] += coeff * _phiU[i] * _uDotGradu[ i % _nComponents ] * jac * _wQuad[k];
      }
    }
  }
}

void feSysElm_VectorConvectiveAcceleration::computeBe(feBilinearForm *form)
{
  double jac, coeff;
  double dotprod;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate scalar coefficient
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    coeff = (*_coeff)(form->_tn, _pos);

    // Get u, grad_u and phiU
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode(form->_sol[_idU], k, _u, _nComponents);
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(form->_sol[_idU], _nComponents, 
      k, form->_transformation, _gradu.data());
    form->_intSpaces[_idU]->getFunctionsAtQuadNode(k, _phiU);

    // Compute (u dot gradu)
    std::fill(_uDotGradu.begin(), _uDotGradu.end(), 0.);
    for(int m = 0; m < _nComponents; ++m) {
      for(int n = 0; n < _nComponents; ++n) {
        _uDotGradu[m] += _u[n] * _gradu[ m*_nComponents + n];
      }
    }

    // Increment with (u dot gradu) dot v
    for(int i = 0; i < _nFunctions; ++i) {
      form->_Be[i] -= coeff * _phiU[i] * _uDotGradu[ i % _nComponents ] * jac * _wQuad[k];
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
  _nFunctions = space[_idU]->getNumFunctions();
  _nComponents = space[_idU]->getNumComponents();
  _gradu.resize(_nComponents * _nComponents);
  _symmetricGradu.resize(_nComponents * _nComponents);
  _gradPhiU.resize(_nComponents * _nFunctions);
}

void feSysElm_DivergenceNewtonianStress::computeAe(feBilinearForm *form)
{
  // TO DO
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
    form->_intSpaces[_idU]->interpolateVectorFieldAtQuadNode_physicalGradient(form->_sol[_idU], _nComponents, 
      k, form->_transformation, _gradu.data());
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation, _gradPhiU.data());

    for(int m = 0; m < _nComponents; ++m) {
      for(int n = 0; n < _nComponents; ++n) {
        _symmetricGradu[m*_nComponents+n] = _gradu[m*_nComponents + n] + _gradu[n*_nComponents + m];
      }
    }

    for(int i = 0; i < _nFunctions; ++i) {

      // Divergence of test function
      div_v = _gradPhiU[i*_nComponents + (i % _nComponents)];
      // Double contraction d(u) : grad(v)
      doubleContraction = 0.;
      for(int m = 0; m < _nComponents; ++m) {
        doubleContraction += _symmetricGradu[(i % _nComponents)*_nComponents + m] * _gradPhiU[i*_nComponents + m];
      }

      form->_Be[i] -= coeff * (-p * div_v + mu * doubleContraction) * jac * _wQuad[k];
    }
  }
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////

void feSysElm_ZeroBlock::createElementarySystem(std::vector<feSpace *> &space)
{
  _idV = 0;
  _idU = 1;
  _fieldsLayoutI = {_idV};
  _fieldsLayoutJ = {_idU};
  _nFunctions = space[_idU]->getNumFunctions();
}