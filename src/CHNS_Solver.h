#ifndef _CHNS_SOLVER_
#define _CHNS_SOLVER_

#include "GenericSolver.h"
#include "feBoundaryConditions.h"
#include "feExporter.h"
#include "feFunction.h"
#include "feMessage.h"
#include "feNonLinearSolver.h"
#include "feParameters.h"

/*

  Generic Cahn-Hilliard Navier-Stokes solver.

*/

class CHNS_Solver;

feStatus
createCHNS_Solver(CHNS_Solver                      *&solver,
                  const Parameters::CHNS            &CHNS_parameters,
                  const Parameters::NonLinearSolver &NLSolver_parameters,
                  const Parameters::TimeIntegration &TimeIntegration_parameters,
                  const std::vector<FEDescriptor>   &fieldDescriptors,
                  const std::vector<BoundaryCondition> &boundaryConditions,
                  PostProcCallback *postProcCallback = nullptr,
                  feFunction       *pSource          = nullptr,
                  feVectorFunction *uSource          = nullptr,
                  feFunction       *phiSource        = nullptr,
                  feFunction       *muSource         = nullptr);

class CHNS_Solver : public SolverBase
{
public:
  Parameters::CHNS _CHNS_parameters;

  feFunction       *_density;
  feFunction       *_drhodphi;
  feFunction       *_viscosity;
  feFunction       *_dviscdphi;
  feFunction       *_mobility;
  feVectorFunction *_volumeForce;

  feFunction       *_pSource;
  feVectorFunction *_uSource;
  feFunction       *_phiSource;
  feFunction       *_muSource;

public:
  CHNS_Solver(const Parameters::CHNS               &CHNS_parameters,
              const Parameters::NonLinearSolver    &NLSolver_parameters,
              const Parameters::TimeIntegration    &TimeIntegration_parameters,
              const std::vector<FEDescriptor>      &fieldDescriptors,
              const std::vector<BoundaryCondition> &boundaryConditions,
              PostProcCallback                     *postProcCallback,
              feFunction                           *pSource,
              feVectorFunction                     *uSource,
              feFunction                           *phiSource,
              feFunction                           *muSource);
  ~CHNS_Solver()
  {
    delete _density;
    delete _drhodphi;
    delete _viscosity;
    delete _dviscdphi;
    delete _mobility;
    delete _volumeForce;

    delete _pSource;
    delete _uSource;
    delete _phiSource;
    delete _muSource;
  }

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

  // /*
  //  * This function computes various postprocessing quantities. Since the
  //  * postprocessing varies a lot from one benchmark to another, the function
  //  * must be implemented individually in each benchmark executable.
  //  * This is not ideal, should think of a better way.
  //  */
  // feStatus computePostProcessing(const std::vector<feSpace *> &spaces,
  //                                feMesh                       *mesh,
  //                                feSolution                   *sol,
  //                                PostProcData &data) const override;

  feStatus updateBeforeFixedPointIteration(const int iter);
};

#endif