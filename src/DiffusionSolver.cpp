
#include "DiffusionSolver.h"

#include "TransientAdapter.h"
#include "feConstants.h"
#include "feTimeIntegration.h"
#include "feNorm.h"

feStatus createDiffusionSolver(
  DiffusionSolver                     *&solver,
  const Parameters::Diffusion          &diff_parameters,
  const Parameters::NonLinearSolver    &NLSolver_parameters,
  const Parameters::TimeIntegration    &TimeIntegration_parameters,
  const std::vector<FEDescriptor>      &fieldDescriptors,
  const std::vector<BoundaryCondition> &boundaryConditions,
  feFunction                           *diffusivity,
  PostProcCallback                     *postProcCallback,
  feFunction                           *source)
{
  if (fieldDescriptors.size() != 1)
  {
    return feErrorMsg(FE_STATUS_ERROR,
                      "Diffusion solver expects 1 unknown finite element "
                      "space, but %d field descriptor(s) were provided.",
                      fieldDescriptors.size());
  }

  solver = new DiffusionSolver(diff_parameters,
                               NLSolver_parameters,
                               TimeIntegration_parameters,
                               fieldDescriptors,
                               boundaryConditions,
                               diffusivity,
                               postProcCallback,
                               source);

  return FE_STATUS_OK;
}

DiffusionSolver::DiffusionSolver(
  const Parameters::Diffusion          &diff_parameters,
  const Parameters::NonLinearSolver    &NLSolver_parameters,
  const Parameters::TimeIntegration    &TimeIntegration_parameters,
  const std::vector<FEDescriptor>      &fieldDescriptors,
  const std::vector<BoundaryCondition> &boundaryConditions,
  feFunction                           *diffusivity,
  PostProcCallback                     *postProcCallback,
  feFunction                           *source)
  : SolverBase(NLSolver_parameters,
               TimeIntegration_parameters,
               fieldDescriptors,
               boundaryConditions,
               postProcCallback)
  , _diff_parameters(diff_parameters)
  , _diffusivity(diffusivity)
  , _source(source)
{
  _numFields = 1;

  if (source == nullptr)
  {
    _source = new feConstantFunction(0.);
  }
}

feStatus
DiffusionSolver::createBilinearForms(const std::vector<feSpace *>  &spaces,
                                     std::vector<feBilinearForm *> &forms) const
{
  forms.clear();

  feSpace *u = spaces[0];

  feBilinearForm *dudt, *diff, *src;
  feCheckReturn(createBilinearForm(
    dudt, {u}, new feSysElm_TransientMass(&scalarConstant::one)));
  feCheckReturn(
    createBilinearForm(diff, {u}, new feSysElm_Diffusion<2>(_diffusivity)));
  feCheckReturn(createBilinearForm(src, {u}, new feSysElm_Source(_source)));
  forms = {dudt, diff, src};

  //
  // TODO: Loop over boundaries and impose flux BC if needed
  //

  return FE_STATUS_OK;
}

//
// Compute error on scalar fields
//
feStatus DiffusionSolver::computeError(const std::vector<feSpace *> &spaces,
                                       feMesh                       */*mesh*/,
                                       feSolution                   *sol,
                                       const int /*iIntervalReferenceSolution*/,
                                       TransientAdapter    &/*adapter*/,
                                       std::vector<double> &errors) const
{
  const feFunction* exactSolution = _scalarExactSolutions[0].second;

  if (exactSolution)
  {
    feSpace *u = spaces[0];
    feNorm *errorU;
    feCheck(createNorm(errorU, L2_ERROR, {u}, sol, exactSolution));
    errors[0] = errorU->compute();
    delete errorU;
  }

  return FE_STATUS_OK;
}