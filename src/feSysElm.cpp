#include "feSysElm.h"
#include "feBilinearForm.h"

// -----------------------------------------------------------------------------
// Source weak form
// -----------------------------------------------------------------------------
void feSysElm_Source::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _feU.resize(_nFunctions);
}

void feSysElm_Source::computeBe(feBilinearForm *form)
{
  for(int k = 0; k < _nQuad; ++k) {
    double jac = form->_J[_nQuad * form->_numElem + k];
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    double S = (*_source)(form->_tn, _pos);
    for(int i = 0; i < _nFunctions; ++i) {
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= _feU[i] * S * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// Transient weak form (mass matrix)
// -----------------------------------------------------------------------------
void feSysElm_Mass::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _feU.resize(_nFunctions);
}

void feSysElm_Mass::computeAe(feBilinearForm *form)
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

void feSysElm_Mass::computeBe(feBilinearForm *form)
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
// Diffusion (stiffness matrix)
// -----------------------------------------------------------------------------
template<int dim>
void feSysElm_Diffusion<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _gradPhi.resize(dim * _nFunctions);
}

template<int dim>
void feSysElm_Diffusion<dim>::computeAe(feBilinearForm *form)
{
  double jac;
  for(int k = 0; k < _nQuad; ++k) {

    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    double kD = (*_diffusivity)(form->_tn, _pos);

    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, 
    	form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, _gradPhi.data());

    for(int i = 0; i < _nFunctions; ++i) {
      for(int j = 0; j < _nFunctions; ++j) {
        for(int iDim = 0; iDim < dim; ++iDim){
          form->_Ae[i][j] += _gradPhi[i*dim + iDim] * _gradPhi[j*dim + iDim] * kD * jac * _wQuad[k];
          // feInfo("nQuad in A = %d - Adding %f", _nQuad, _gradPhi[i*dim + iDim] * _gradPhi[j*dim + iDim] * kD * jac * _wQuad[k]);
        }
      }
    }
  }
}

template<int dim>
void feSysElm_Diffusion<dim>::computeBe(feBilinearForm *form)
{
  double jac, grad_u[3] = {0., 0., 0.};
  for(int k = 0; k < _nQuad; ++k) {

    jac = form->_J[_nQuad * form->_numElem + k];

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    double kD = (*_diffusivity)(form->_tn, _pos);

    form->_intSpaces[_idU]->interpolateFieldAtQuadNode_physicalGradient(form->_sol[_idU], k, jac, form->_geoCoord, 
    	form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, grad_u);

    // double res = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k)/jac;

    // feInfo("Error =  %1.6e", fabs(grad_u[0] - res));

    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, 
    	form->_dxdr, form->_dxds, form->_dxdt, form->_drdx, form->_drdy, form->_drdz, _gradPhi.data());
    
    for(int i = 0; i < _nFunctions; ++i) {
    	for(int iDim = 0; iDim < dim; ++iDim){
      	form->_Be[i] -= _gradPhi[i*dim + iDim] * grad_u[iDim] * kD * jac * _wQuad[k];
    	}
    }
  }
}

template class feSysElm_Diffusion<0>;
template class feSysElm_Diffusion<1>;
template class feSysElm_Diffusion<2>;