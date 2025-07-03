
#include "DiffusionSolver.h"

#include "TransientAdapter.h"
#include "feConstants.h"
#include "feNorm.h"
#include "feTimeIntegration.h"

feStatus createDiffusionSolver(
  DiffusionSolver                  *&solver,
  const Parameters::Diffusion       &diff_parameters,
  const Parameters::NonLinearSolver &NLSolver_parameters,
  const Parameters::TimeIntegration &TimeIntegration_parameters,
  const std::vector<FEDescriptor>   &fieldDescriptors,
  const std::vector<const BoundaryConditions::BoundaryCondition *>
                   &boundaryConditions,
  feFunction       *diffusivity,
  PostProcCallback *postProcCallback,
  feFunction       *source,
  feVectorFunction *vectorSource)
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
                               source,
                               vectorSource);

  return FE_STATUS_OK;
}

DiffusionSolver::DiffusionSolver(
  const Parameters::Diffusion       &diff_parameters,
  const Parameters::NonLinearSolver &NLSolver_parameters,
  const Parameters::TimeIntegration &TimeIntegration_parameters,
  const std::vector<FEDescriptor>   &fieldDescriptors,
  const std::vector<const BoundaryConditions::BoundaryCondition *>
                   &boundaryConditions,
  feFunction       *diffusivity,
  PostProcCallback *postProcCallback,
  feFunction       *source,
  feVectorFunction *vectorSource)
  : SolverBase(NLSolver_parameters,
               TimeIntegration_parameters,
               fieldDescriptors,
               boundaryConditions,
               postProcCallback)
  , _diff_parameters(diff_parameters)
  , _diffusivity(diffusivity)
  , _source(source)
  , _vectorSource(vectorSource)
{
  _numFields = 1;

  if (source == nullptr)
  {
    _source = new feConstantFunction(0.);
  }
  if (vectorSource == nullptr)
  {
    vectorSource = new feConstantVectorFunction({0., 0., 0.});
  }
}

feStatus
DiffusionSolver::createBilinearForms(const std::vector<feSpace *>  &spaces,
                                     std::vector<feBilinearForm *> &forms) const
{
  forms.clear();

  feSpace *u = spaces[0];

  feBilinearForm *dudt, *diff, *src;

  if(u->isVectorValued())
  {
    feSysElm *f1 = new feSysElm_TransientVectorMass<2>(&scalarConstant::one);
    feSysElm *f2 = new feSysElm_VectorDiffusion<2>(&scalarConstant::one, _diffusivity);
    feSysElm *f3 = new feSysElm_VectorSource<2>(_vectorSource);

    feCheckReturn(createBilinearForm(dudt, {u}, f1));
    feCheckReturn(createBilinearForm(diff, {u}, f2));
    feCheckReturn(createBilinearForm( src, {u}, f3));
  } else {
    feSysElm *f1 = new feSysElm_TransientMass(&scalarConstant::one);
    feSysElm *f2 = new feSysElm_Diffusion<2>(_diffusivity);
    feSysElm *f3 = new feSysElm_Source(_source);

    feCheckReturn(createBilinearForm(dudt, {u}, f1));
    feCheckReturn(createBilinearForm(diff, {u}, f2));
    feCheckReturn(createBilinearForm( src, {u}, f3));
  }
  
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
                                       feMesh                       *mesh,
                                       feSolution                   *sol,
                                       const int iIntervalReferenceSolution,
                                       TransientAdapter    &adapter,
                                       std::vector<double> &errors) const
{
  const feFunction *exactSolution = _scalarExactSolutions[0].second;

  feSpace *u      = spaces[0];
  feNorm  *errorU = nullptr;

  if (exactSolution)
  {
    feInfo("Computing error from analytical solution");
    feCheckReturn(createNorm(errorU, L2_ERROR, {u}, sol, exactSolution));
    errors[0] = errorU->compute();
  }

  if (adapter._adapt_parameters.readReferenceSolution)
  {
    feInfo("Computing error from finer reference solution");
    feCheckReturn(createNorm(errorU, NONE, {u}, sol));
    feCheckReturn(errorU->computeErrorNormFromExternalSolution(
      adapter.referenceTestCases[iIntervalReferenceSolution].mesh,
      adapter.referenceTestCases[iIntervalReferenceSolution].sol,
      adapter.referenceTestCases[iIntervalReferenceSolution].spaces,
      mesh,
      sol,
      errors));
  }

  delete errorU;

  return FE_STATUS_OK;
}