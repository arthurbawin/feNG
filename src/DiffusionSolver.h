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
  DiffusionSolver                     *&solver,
  const Parameters::Diffusion          &diff_parameters,
  const Parameters::NonLinearSolver    &NLSolver_parameters,
  const Parameters::TimeIntegration    &TimeIntegration_parameters,
  const std::vector<FEDescriptor>      &fieldDescriptors,
  const std::vector<BoundaryCondition> &boundaryConditions,
  feFunction                           *diffusivity,
  PostProcCallback                     *postProcCallback = nullptr,
  feFunction                           *source           = nullptr);

class DiffusionSolver : public SolverBase
{
public:
  Parameters::Diffusion _diff_parameters;

  feFunction *_diffusivity;
  feFunction *_source;

public:
  DiffusionSolver(const Parameters::Diffusion       &diff_parameters,
                  const Parameters::NonLinearSolver &NLSolver_parameters,
                  const Parameters::TimeIntegration &TimeIntegration_parameters,
                  const std::vector<FEDescriptor>   &fieldDescriptors,
                  const std::vector<BoundaryCondition> &boundaryConditions,
                  feFunction                           *diffusivity,
                  PostProcCallback                     *postProcCallback,
                  feFunction                           *source);
  ~DiffusionSolver() {}

  feStatus solve(feMesh                 *mesh,
                 feSolution             *sol,
                 feMetaNumber           *numbering,
                 std::vector<feSpace *> &spaces,
                 feExportData           &exportData) const override;

  feStatus solve(feMesh                 *mesh,
                 feSolution             *sol,
                 feMetaNumber           *numbering,
                 std::vector<feSpace *> &spaces,
                 feExportData           &exportData,
                 TransientAdapter       &adapter,
                 const int               iInterval) const override;

  //
  // Compute error on scalar fields p, phi, mu w.r.t. reference finer solution
  //
  feStatus computeError(const std::vector<feSpace *> &spaces,
                        feMesh                       *mesh,
                        feSolution                   *sol,
                        const int            iIntervalReferenceSolution,
                        TransientAdapter    &adapter,
                        std::vector<double> &errors) const;

  feStatus updateBeforeFixedPointIteration(const int iter);
};

#endif