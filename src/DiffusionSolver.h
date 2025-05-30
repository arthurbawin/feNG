#ifndef _DIFFUSION_SOLVER_
#define _DIFFUSION_SOLVER_

#include "GenericSolver.h"
#include "feBoundaryConditions.h"
#include "feExporter.h"
#include "feFunction.h"
#include "feMessage.h"
#include "feNonLinearSolver.h"
#include "feParameters.h"

/*

  Generic solver for the unsteady heat equation:

  du
  -- - ∇ . (k ∇ u) + f = 0
  dt

*/

class DiffusionSolver;

feStatus createDiffusionSolver(
  DiffusionSolver                  *&solver,
  const Parameters::Diffusion       &diff_parameters,
  const Parameters::NonLinearSolver &NLSolver_parameters,
  const Parameters::TimeIntegration &TimeIntegration_parameters,
  const std::vector<FEDescriptor>   &fieldDescriptors,
  const std::vector<const BoundaryConditions::BoundaryCondition *>
                   &boundaryConditions,
  feFunction       *diffusivity,
  PostProcCallback *postProcCallback = nullptr,
  feFunction       *source           = nullptr);

class DiffusionSolver : public SolverBase
{
public:
  Parameters::Diffusion _diff_parameters;

  feFunction *_diffusivity;
  feFunction *_source;

public:
  DiffusionSolver(
    const Parameters::Diffusion       &diff_parameters,
    const Parameters::NonLinearSolver &NLSolver_parameters,
    const Parameters::TimeIntegration &TimeIntegration_parameters,
    const std::vector<FEDescriptor>   &fieldDescriptors,
    const std::vector<const BoundaryConditions::BoundaryCondition *>
                     &boundaryConditions,
    feFunction       *diffusivity,
    PostProcCallback *postProcCallback,
    feFunction       *source);
  ~DiffusionSolver() {}

  feStatus createBilinearForms(const std::vector<feSpace *>  &spaces,
                               std::vector<feBilinearForm *> &forms) const;

  //
  // Compute error on scalar field u
  //
  feStatus computeError(const std::vector<feSpace *> &spaces,
                        feMesh                       *mesh,
                        feSolution                   *sol,
                        const int            iIntervalReferenceSolution,
                        TransientAdapter    &adapter,
                        std::vector<double> &errors) const;
};

#endif