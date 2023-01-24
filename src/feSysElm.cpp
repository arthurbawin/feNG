#include "feSysElm.h"
#include "feBilinearForm.h"

// -----------------------------------------------------------------------------
// Source weak form
// -----------------------------------------------------------------------------
template<int dim>
void feSysElm_Source<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _feU.resize(_nFunctions);
}

template<int dim>
void feSysElm_Source<dim>::computeBe(feBilinearForm *form)
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
// Transient term (mass matrix)
// -----------------------------------------------------------------------------
template<int dim>
void feSysElm_Mass<dim>::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _feU.resize(_nFunctions);
}

template<int dim>
void feSysElm_Mass<dim>::computeAe(feBilinearForm *form)
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

template<int dim>
void feSysElm_Mass<dim>::computeBe(feBilinearForm *form)
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
  _feU.resize(_nFunctions);
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

    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, _gradPhi.data());

    for(int i = 0; i < _nFunctions; ++i) {
      for(int j = 0; j < _nFunctions; ++j) {
        for(int iDim = 0; iDim < dim; ++iDim){
          form->_Ae[i][j] += _gradPhi[iDim * _nFunctions + i] * _gradPhi[iDim * _nFunctions + j] * kD * jac * _wQuad[k];
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
    form->_intSpaces[_idU]->interpolateFieldAtQuadNode_physicalGradient(form->_sol[_idU], k, jac, form->_geoCoord, grad_u);
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, jac, form->_geoCoord, _gradPhi.data());
    for(int i = 0; i < _nFunctions; ++i) {
    	for(int iDim = 0; iDim < dim; ++iDim){
      	form->_Be[i] -= _gradPhi[iDim * _nFunctions + i] * grad_u[iDim] * kD * jac * _wQuad[k];
    	}
    }
  }
}

template class feSysElm_Source<0>;
template class feSysElm_Source<1>;
template class feSysElm_Source<2>;
template class feSysElm_Mass<0>;
template class feSysElm_Mass<1>;
template class feSysElm_Mass<2>;
template class feSysElm_Diffusion<0>;
template class feSysElm_Diffusion<1>;
template class feSysElm_Diffusion<2>;