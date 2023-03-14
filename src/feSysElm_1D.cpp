#include "feSysElm.h"
#include "feSysElm_Ext.h"
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
  double u, dudx, jac, cVelocity;
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];

    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);
    dudx = form->_intSpaces[_idU]->interpolateFieldAtQuadNode_rDerivative(form->_sol[_idU], k);
    dudx /= jac;

    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    cVelocity = (*_velocity)(form->_tn, _pos);

    // Compute element Peclet and optimal numerical diffusivity coth(Peh) - 1/Peh
    double he = 2.0 * jac;
    double Peh = cVelocity * he / (2. * _diffusivity);
    double beta = cosh(Peh) / sinh(Peh) - (1. / Peh);
    double tau = beta * he / (2. * cVelocity);
    // FIXME : completer avec d2dx2
    // FIXME : choisir le signe du terme source
    double residual = cVelocity * dudx - (*_source)(form->_tn, _pos);

    for(int i = 0; i < _nFunctions; ++i) {
      _feUdx[i] = form->_intSpaces[_idU]->getdFunctiondrAtQuadNode(i, k);
      _feUdx[i] /= jac;

      form->_Be[i] -= residual * tau * cVelocity * _feUdx[i] * jac * _wQuad[k];
    }
  }
}