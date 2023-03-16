#include "feSysElm.h"
#include "feBilinearForm.h"

// -----------------------------------------------------------------------------
// 1D Neumann boundary condition (= source)
// -----------------------------------------------------------------------------
void feSysElm_1D_NeumannBC::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _feU.resize(_nFunctions);
}

void feSysElm_1D_NeumannBC::computeBe(feBilinearForm *form)
{
  for(int k = 0; k < _nQuad; ++k) {
    double jac = form->_J[_nQuad * form->_numElem + k];
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    double h = _neumannBC->eval(form->_tn, _pos);
    for(int i = 0; i < _nFunctions; ++i) {
      _feU[i] = form->_intSpaces[_idU]->getFunctionAtQuadNode(i, k);
      form->_Be[i] -= _feU[i] * h * jac * _wQuad[k];
    }
  }
}

// -----------------------------------------------------------------------------
// 1D SUPG Stabilization
// -----------------------------------------------------------------------------
void feSysElm_1D_SUPGStab::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _feU.resize(_nFunctions);
  _feUdx.resize(_nFunctions);
}

void feSysElm_1D_SUPGStab::computeAe(feBilinearForm *form)
{
  // Finite differences
}

void feSysElm_1D_SUPGStab::computeBe(feBilinearForm *form)
{
  // ...
}