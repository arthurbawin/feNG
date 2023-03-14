#include "feSysElm.h"
#include "feBilinearForm.h"

// -----------------------------------------------------------------------------
// 2D SUPG Stabilization
// -----------------------------------------------------------------------------
void feSysElm_2D_SUPGStab::createElementarySystem(std::vector<feSpace *> &space)
{
  _idU = 0;
  _fieldsLayoutI[0] = _idU;
  _fieldsLayoutJ[0] = _idU;
  _nFunctions = space[_idU]->getNumFunctions();
  _phiU.resize(_nFunctions);
  _gradPhiU.resize(2 * _nFunctions);
}

void feSysElm_2D_SUPGStab::computeAe(feBilinearForm *form)
{
  // Finite differences
}

static double tau2DTri(int dim, double velocity[3], double diffusivity, int nFunctions,
                       std::vector<double> &gradphi)
{
  // Compute hTau from (13.17) in Fortin & Garon
  double normV =
    sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]);

  double res = 0., dotprod;
  for(int i = 0; i < nFunctions; ++i) {
    dotprod = 0.;
    for(int iDim = 0; iDim < dim; ++iDim) {
      dotprod += velocity[iDim] / normV * gradphi[i * dim + iDim];
    }
    res += dotprod * dotprod;
  }

  double hTau = sqrt(2.) / sqrt(res);

  double rho = 1.;
  double cp = 1.;
  return 1.0 / sqrt(4. * normV * normV / (hTau * hTau) +
                    144. * diffusivity * diffusivity /
                      ((rho * cp * hTau * hTau) * (rho * cp * hTau * hTau)));
}

void feSysElm_2D_SUPGStab::computeBe(feBilinearForm *form)
{
  double jac, u, grad_u[3] = {0., 0., 0.}, tau, diffusivity, source, reaction;
  std::vector<double> c(2, 0.0);
  for(int k = 0; k < _nQuad; ++k) {
    jac = form->_J[_nQuad * form->_numElem + k];
    form->_cnc->computeElementTransformation(form->_geoCoord, k, jac, form->_transformation);

    // Evaluate the imposed velocity field, diffusivity, reaction coefficient and source term
    form->_geoSpace->interpolateVectorFieldAtQuadNode(form->_geoCoord, k, _pos);
    (*_velocity)(form->_tn, _pos, c);
    diffusivity = (*_diffusivity)(form->_tn, _pos);
    reaction = (*_reactionCoeff)(form->_tn, _pos);
    source = (*_source)(form->_tn, _pos);

    // Compute u
    u = form->_intSpaces[_idU]->interpolateFieldAtQuadNode(form->_sol[_idU], k);

    // Compute grad(u)
    form->_intSpaces[_idU]->interpolateFieldAtQuadNode_physicalGradient(
      form->_sol[_idU], k, form->_transformation, grad_u);

    // Compute grad(phi)
    form->_intSpaces[_idU]->getFunctionsPhysicalGradientAtQuadNode(k, form->_transformation,
                                                                   _gradPhiU.data());

    // Compute stabilization parameter tau
    tau = tau2DTri(2, c.data(), diffusivity, _nFunctions, _gradPhiU);

    // Compute PDE residual:
    // Reaction term
    double residual = reaction * u;
    // Convection term
    for(int iDim = 0; iDim < 2; ++iDim) {
      residual += c[iDim] * grad_u[iDim];
    }
    // Diffusion term (0 for linear elements)
    // ...
    /// Source term
    residual += source;

    for(int i = 0; i < _nFunctions; ++i) {
      for(int iDim = 0; iDim < 2; ++iDim) {
        form->_Be[i] -= residual * tau * c[iDim] * _gradPhiU[i * 2 + iDim] * jac * _wQuad[k];
      }
    }
  }
}