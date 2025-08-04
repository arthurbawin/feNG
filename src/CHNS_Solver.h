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

feStatus createCHNS_Solver(
  CHNS_Solver                      *&solver,
  const Parameters::CHNS            &CHNS_parameters,
  const Parameters::NonLinearSolver &NLSolver_parameters,
  const Parameters::TimeIntegration &TimeIntegration_parameters,
  const std::vector<FEDescriptor>   &fieldDescriptors,
  const std::vector<const BoundaryConditions::BoundaryCondition *>
                   &boundaryConditions,
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

  feConstantFunction *_coeffFluxContactAngle;

  feFunction       *_pSource;
  feVectorFunction *_uSource;
  feFunction       *_phiSource;
  feFunction       *_muSource;

  // The indicator function for each parameter,
  // returning rho(phi) or mu(phi) if in fluid A (resp. B)
  // and 0 otherwise.
  feFunction *_densityIndicatorFluidA;
  feFunction *_densityIndicatorFluidB;
  feFunction *_viscosityIndicatorFluidA;
  feFunction *_viscosityIndicatorFluidB;

public:
  CHNS_Solver(const Parameters::CHNS            &CHNS_parameters,
              const Parameters::NonLinearSolver &NLSolver_parameters,
              const Parameters::TimeIntegration &TimeIntegration_parameters,
              const std::vector<FEDescriptor>   &fieldDescriptors,
              const std::vector<const BoundaryConditions::BoundaryCondition *>
                               &boundaryConditions,
              PostProcCallback *postProcCallback,
              feFunction       *pSource,
              feVectorFunction *uSource,
              feFunction       *phiSource,
              feFunction       *muSource);
  ~CHNS_Solver()
  {
    delete _density;
    delete _drhodphi;
    delete _viscosity;
    delete _dviscdphi;
    delete _mobility;
    delete _volumeForce;
    
    delete _densityIndicatorFluidA;
    delete _densityIndicatorFluidB;
    delete _viscosityIndicatorFluidA;
    delete _viscosityIndicatorFluidB;
  }

  double getInterfaceThickness() const { return _CHNS_parameters.epsilon; }

  // Get the const density or viscosity callback or their
  // indicator function.
  // This is used in post-processing callbacks, to access
  // the rho(phi) and mu(phi) relations.
  const feFunction *densityCallback() const { return _density; };
  const feFunction *dynamicViscosityCallback() const { return _viscosity; };
  const feFunction *densityIndicatorCallback(const int whichFluid) const
  {
    return (whichFluid == 0) ? _densityIndicatorFluidA :
                               _densityIndicatorFluidB;
  };
  const feFunction *
  dynamicViscosityIndicatorCallback(const int whichFluid) const
  {
    return (whichFluid == 0) ? _viscosityIndicatorFluidA :
                               _viscosityIndicatorFluidB;
  };

  feStatus createBilinearForms(const std::vector<feSpace *>  &spaces,
                               std::vector<feBilinearForm *> &forms) const;

  void updateBeforeFixedPointIteration(const int iter) override;

  //
  // Compute error on scalar fields p, phi, mu w.r.t. reference finer solution
  //
  feStatus computeError(const std::vector<feSpace *> &spaces,
                        feMesh                       *mesh,
                        feSolution                   *sol,
                        const int            iIntervalReferenceSolution,
                        TransientAdapter    &adapter,
                        std::vector<double> &errors) const;
};

#endif