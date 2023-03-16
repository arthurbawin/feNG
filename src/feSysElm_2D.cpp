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
  // ...
}