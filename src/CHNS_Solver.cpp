
#include "CHNS_Solver.h"

#include "TransientAdapter.h"
#include "feTimeIntegration.h"

feStatus
createCHNS_Solver(CHNS_Solver                      *&solver,
                  const Parameters::CHNS            &CHNS_parameters,
                  const Parameters::NonLinearSolver &NLSolver_parameters,
                  const Parameters::TimeIntegration &TimeIntegration_parameters,
                  const std::vector<FEDescriptor>   &fieldDescriptors,
                  const std::vector<BoundaryCondition> &boundaryConditions,
                  PostProcCallback                     *postProcCallback,
                  feFunction                           *pSource,
                  feVectorFunction                     *uSource,
                  feFunction                           *phiSource,
                  feFunction                           *muSource)
{
  if (fieldDescriptors.size() != 4)
  {
    return feErrorMsg(FE_STATUS_ERROR,
                      "CHNS Solver expects 4 unknown finite element spaces: u, "
                      "p, phi, mu. Only %d field descriptor(s) were provided.",
                      fieldDescriptors.size());
  }

  std::string expectedName[4] = {"U", "P", "Phi", "Mu"};

  for (size_t i = 0; i < fieldDescriptors.size(); ++i)
  {
    if (fieldDescriptors[i]._fieldName != expectedName[i])
    {
      return feErrorMsg(FE_STATUS_ERROR,
                        "FE space %d: Expected field name: %s - Received: %s",
                        i,
                        expectedName[i].data(),
                        fieldDescriptors[i]._fieldName.data());
    }
  }

  solver = new CHNS_Solver(CHNS_parameters,
                           NLSolver_parameters,
                           TimeIntegration_parameters,
                           fieldDescriptors,
                           boundaryConditions,
                           postProcCallback,
                           pSource,
                           uSource,
                           phiSource,
                           muSource);

  return FE_STATUS_OK;
}

double density_f(const feFunctionArguments &args,
                 const std::vector<double> &par)
{
  // Assumes phi is in [-1,1]
  // Property of A when phi = 1, of B when phi = -1
  double phi  = args.u;
  double rhoA = par[0];
  double rhoB = par[1];
  return (rhoA - rhoB) / 2. * phi + (rhoA + rhoB) / 2.;
}

double densityLimiter_f(const feFunctionArguments &args,
                        const std::vector<double> &par)
{
  double phi  = fmax(-1., fmin(1., args.u));
  double rhoA = par[0];
  double rhoB = par[1];
  return (rhoA - rhoB) / 2. * phi + (rhoA + rhoB) / 2.;
}

double drhodphi_f(const feFunctionArguments & /*args*/,
                  const std::vector<double> &par)
{
  double rhoA = par[0];
  double rhoB = par[1];
  return (rhoA - rhoB) / 2.;
}

double viscosity_f(const feFunctionArguments &args,
                   const std::vector<double> &par)
{
  double phi   = args.u;
  double viscA = par[0];
  double viscB = par[1];
  return (viscA - viscB) / 2. * phi + (viscA + viscB) / 2.;
}

double viscosityLimiter_f(const feFunctionArguments &args,
                          const std::vector<double> &par)
{
  double phi   = fmax(-1., fmin(1., args.u));
  double viscA = par[0];
  double viscB = par[1];
  return (viscA - viscB) / 2. * phi + (viscA + viscB) / 2.;
}

double dviscdphi_f(const feFunctionArguments & /*args*/,
                   const std::vector<double> &par)
{
  double viscA = par[0];
  double viscB = par[1];
  return (viscA - viscB) / 2.;
}

double degenerateMobility_f(const feFunctionArguments &args,
                            const std::vector<double> &par)
{
  double       phi         = args.u;
  const double mobilityVal = par[0];
  return mobilityVal * fabs(1. - phi * phi);
}

CHNS_Solver::CHNS_Solver(
  const Parameters::CHNS               &CHNS_parameters,
  const Parameters::NonLinearSolver    &NLSolver_parameters,
  const Parameters::TimeIntegration    &TimeIntegration_parameters,
  const std::vector<FEDescriptor>      &fieldDescriptors,
  const std::vector<BoundaryCondition> &boundaryConditions,
  PostProcCallback                     *postProcCallback,
  feFunction                           *pSource,
  feVectorFunction                     *uSource,
  feFunction                           *phiSource,
  feFunction                           *muSource)
  : SolverBase(NLSolver_parameters,
               TimeIntegration_parameters,
               fieldDescriptors,
               boundaryConditions,
               postProcCallback)
  , _CHNS_parameters(CHNS_parameters)
  , _pSource(pSource)
  , _uSource(uSource)
  , _phiSource(phiSource)
  , _muSource(muSource)
{
  const double rhoA = CHNS_parameters.fluids[0].density;
  const double rhoB = CHNS_parameters.fluids[1].density;

  const double viscA = CHNS_parameters.fluids[0].dynamic_viscosity;
  const double viscB = CHNS_parameters.fluids[1].dynamic_viscosity;

  const double M = CHNS_parameters.mobility;

  //
  // Initialize CHNS callbacks
  //
  if (CHNS_parameters.limitPhaseMarker)
  {
    _density   = new feFunction(densityLimiter_f, {rhoA, rhoB});
    _viscosity = new feFunction(viscosityLimiter_f, {viscA, viscB});
  }
  else
  {
    _density   = new feFunction(density_f, {rhoA, rhoB});
    _viscosity = new feFunction(viscosity_f, {viscA, viscB});
  }

  _drhodphi  = new feFunction(drhodphi_f, {rhoA, rhoB});
  _dviscdphi = new feFunction(dviscdphi_f, {viscA, viscB});

  switch (CHNS_parameters.mobilityType)
  {
    case Parameters::CHNS::MobilityType::constant:
      _mobility = new feConstantFunction({M});
      break;
    case Parameters::CHNS::MobilityType::degenerate:
      _mobility = new feFunction(degenerateMobility_f, {M});
      break;
  }

  _volumeForce =
    new feConstantVectorFunction({0., CHNS_parameters.gravity, 0.});

  if (pSource == nullptr)
  {
    _pSource   = new feConstantFunction(0.);
    _uSource   = new feConstantVectorFunction({0., 0., 0.});
    _phiSource = new feConstantFunction(0.);
    _muSource  = new feConstantFunction(0.);
  }

  //
  // Initialize derived coefficients
  //
  _CHNS_parameters.mass_alpha = (rhoB - rhoA) / (rhoA + rhoB);
}

feStatus CHNS_Solver::solve(feMesh                 *mesh,
                            feSolution             *sol,
                            feMetaNumber           *numbering,
                            std::vector<feSpace *> &spaces,
                            feExportData           &exportData) const
{
  //
  // Time stepping is controlled by the TimeIntegration parameters
  //
  const int    nT = _TimeIntegration_parameters.nTimeSteps;
  const double t0 = _TimeIntegration_parameters.t_initial;
  const double tf = _TimeIntegration_parameters.t_final;

  feSpace *u   = spaces[0];
  feSpace *p   = spaces[1];
  feSpace *phi = spaces[2];
  feSpace *mu  = spaces[3];

  std::vector<feBilinearForm *> forms = {};

  std::vector<double> scalarParameters;

  feBilinearForm *CHNS;

  switch (_CHNS_parameters.formulation)
  {
    case Parameters::CHNS::Formulation::volumeAveraged:

      scalarParameters.push_back(_CHNS_parameters.surfaceTension);
      scalarParameters.push_back(_CHNS_parameters.epsilon);

      feCheckReturn(
        createBilinearForm(CHNS,
                           {u, p, phi, mu},
                           new CHNS_VolumeAveraged<2>(_density,
                                                      _drhodphi,
                                                      _viscosity,
                                                      _dviscdphi,
                                                      _mobility,
                                                      _volumeForce,
                                                      _pSource,
                                                      _uSource,
                                                      _phiSource,
                                                      _muSource,
                                                      scalarParameters)));
      break;

    case Parameters::CHNS::Formulation::massAveraged:

      scalarParameters.push_back(_CHNS_parameters.mass_alpha);
      scalarParameters.push_back(_CHNS_parameters.surfaceTension);
      scalarParameters.push_back(_CHNS_parameters.epsilon);

      feCheckReturn(
        createBilinearForm(CHNS,
                           {u, p, phi, mu},
                           new feSysElm_CHNS_Alternative<2>(_density,
                                                            _drhodphi,
                                                            _viscosity,
                                                            _dviscdphi,
                                                            _mobility,
                                                            _volumeForce,
                                                            _pSource,
                                                            _uSource,
                                                            _phiSource,
                                                            _muSource,
                                                            scalarParameters)));
      break;
  }

  forms.push_back(CHNS);

  //
  // TODO: Loop over boundaries and impose flux BC if needed
  //

  feLinearSystem *system = nullptr;
#if defined(HAVE_MKL)
  feCheckReturn(
    createLinearSystem(system, MKLPARDISO, forms, numbering->getNbUnknowns()));
#elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
  feCheckReturn(
    createLinearSystem(system, PETSC_MUMPS, forms, numbering->getNbUnknowns()));
#else
  feCheckReturn(
    createLinearSystem(system, PETSC, forms, numbering->getNbUnknowns()));
#endif

  std::vector<feNorm *> norms = {};

  feNLSolverOptions NLoptions{
    _NLSolver_parameters.tolResidual,
    _NLSolver_parameters.tolCorrection,
    _NLSolver_parameters.tolDivergence,
    _NLSolver_parameters.maxIter,
    _NLSolver_parameters.recomputeJacobianEveryNsteps,
    _NLSolver_parameters.residualDecrease,
  };

  TimeIntegrator *timeIntegrator = nullptr;
  feCheckReturn(createTimeIntegrator(timeIntegrator,
                                     timeIntegratorScheme::BDF2,
                                     NLoptions,
                                     system,
                                     sol,
                                     mesh,
                                     norms,
                                     exportData,
                                     t0,
                                     tf,
                                     nT));

  static_cast<BDF2Integrator*>(timeIntegrator)->setStartingMethod(BDF2Starter::InitialCondition);

  // Solve
  feCheckReturn(timeIntegrator->makeSteps(nT));

  delete timeIntegrator;
  delete system;
  for (feBilinearForm *f : forms)
    delete f;

  return FE_STATUS_OK;
}

feStatus CHNS_Solver::solve(feMesh                 *mesh,
                            feSolution             *sol,
                            feMetaNumber           *numbering,
                            std::vector<feSpace *> &spaces,
                            feExportData           &exportData,
                            TransientAdapter       &adapter,
                            const int               iInterval) const
{
  //
  // Time stepping is controlled by the transient adapter
  //
  const int    nT = adapter._adapt_parameters.nTimeStepsPerInterval;
  const double t0 = adapter.currentTime;
  const double tf = t0 + adapter._time_parameters.dt * nT;
  const double dt = adapter._time_parameters.dt;

  feSpace *u   = spaces[0];
  feSpace *p   = spaces[1];
  feSpace *phi = spaces[2];
  feSpace *mu  = spaces[3];

  std::vector<feBilinearForm *> forms = {};

  std::vector<double> scalarParameters;

  feBilinearForm *CHNS;

  switch (_CHNS_parameters.formulation)
  {
    case Parameters::CHNS::Formulation::volumeAveraged:

      scalarParameters.push_back(_CHNS_parameters.surfaceTension);
      scalarParameters.push_back(_CHNS_parameters.epsilon);

      feCheckReturn(
        createBilinearForm(CHNS,
                           {u, p, phi, mu},
                           new CHNS_VolumeAveraged<2>(_density,
                                                      _drhodphi,
                                                      _viscosity,
                                                      _dviscdphi,
                                                      _mobility,
                                                      _volumeForce,
                                                      _pSource,
                                                      _uSource,
                                                      _phiSource,
                                                      _muSource,
                                                      scalarParameters)));
      break;

    case Parameters::CHNS::Formulation::massAveraged:

      scalarParameters.push_back(_CHNS_parameters.mass_alpha);
      scalarParameters.push_back(_CHNS_parameters.surfaceTension);
      scalarParameters.push_back(_CHNS_parameters.epsilon);

      feCheckReturn(
        createBilinearForm(CHNS,
                           {u, p, phi, mu},
                           new feSysElm_CHNS_Alternative<2>(_density,
                                                            _drhodphi,
                                                            _viscosity,
                                                            _dviscdphi,
                                                            _mobility,
                                                            _volumeForce,
                                                            _pSource,
                                                            _uSource,
                                                            _phiSource,
                                                            _muSource,
                                                            scalarParameters)));
      break;
  }

  forms.push_back(CHNS);

  //
  // TODO: Loop over boundaries and impose flux BC if needed
  //

  feLinearSystem *system = nullptr;
#if defined(HAVE_MKL)
  feCheckReturn(
    createLinearSystem(system, MKLPARDISO, forms, numbering->getNbUnknowns()));
#elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
  feCheckReturn(
    createLinearSystem(system, PETSC_MUMPS, forms, numbering->getNbUnknowns()));
#else
  feCheckReturn(
    createLinearSystem(system, PETSC, forms, numbering->getNbUnknowns()));
#endif

  std::vector<feNorm *> norms = {};

  feNLSolverOptions NLoptions{
    _NLSolver_parameters.tolResidual,
    _NLSolver_parameters.tolCorrection,
    _NLSolver_parameters.tolDivergence,
    _NLSolver_parameters.maxIter,
    _NLSolver_parameters.recomputeJacobianEveryNsteps,
    _NLSolver_parameters.residualDecrease,
  };

  TimeIntegrator *timeIntegrator = nullptr;
  feCheckReturn(createTimeIntegrator(timeIntegrator,
                                     timeIntegratorScheme::BDF2,
                                     NLoptions,
                                     system,
                                     sol,
                                     mesh,
                                     norms,
                                     exportData,
                                     t0,
                                     tf,
                                     nT));

  // Restart from solution container if not the first sub-interval
  // (restart from last solutions of previous sub-interval)
  if (iInterval > 0)
  {
    feCheckReturn(
      timeIntegrator->restartFromContainer(*adapter.allContainers[iInterval]));
    timeIntegrator->setCurrentStep(adapter.currentTimeStep);

    feInfo("Time           of timeIntegrator : %f",
           timeIntegrator->getCurrentTime());
    feInfo("Previous times of timeIntegrator : %f",
           timeIntegrator->getTime()[0]);
    feInfo("Previous times of timeIntegrator : %f",
           timeIntegrator->getTime()[1]);
    feInfo("Previous times of timeIntegrator : %f",
           timeIntegrator->getTime()[2]);
    feInfo("Step           of timeIntegrator : %d",
           timeIntegrator->getCurrentStep());
  }

  // Solve and compute metrics at each time step
  for (int iStep = 0; iStep < nT; ++iStep)
  {
    tic();
    feCheckReturn(timeIntegrator->makeSteps(1));
    feInfo("Computed time step in %f s", toc());

    // Reconstruct phi
    std::string   meshFile = adapter.allOptions[iInterval].backgroundMeshfile;
    std::string   recoveryFile = adapter._io_parameters.writeDir + "phi.msh";
    feNewRecovery recPhi(phi,
                         0,
                         mesh,
                         sol,
                         meshFile,
                         recoveryFile,
                         false,
                         false,
                         false,
                         nullptr,
                         numbering);

    // Compute metric field and increment metric integral Hi
    feMetric *metricField = adapter.allMetrics[iInterval];
    metricField->setRecoveredFields({&recPhi});
    metricField->setCurrentTime(timeIntegrator->getCurrentTime());
    metricField->setMetricScaling(false);
    feCheckReturn(metricField->computeMetrics());

    // Trapeze rule
    if (iStep == 0 || iStep == nT - 1)
      metricField->addMetricsToOther(dt / 2., adapter.allHi[iInterval]);
    else
      metricField->addMetricsToOther(dt, adapter.allHi[iInterval]);

    adapter.currentTimeStep++;
  }

  adapter.currentTime = timeIntegrator->getCurrentTime();

  delete adapter.allContainers[iInterval];
  adapter.allContainers[iInterval] = new feSolutionContainer();
  adapter.allContainers[iInterval]->NaNify();
  *(adapter.allContainers[iInterval]) = timeIntegrator->getSolutionContainer();

  delete timeIntegrator;
  delete system;
  for (feBilinearForm *f : forms)
    delete f;

  return FE_STATUS_OK;
}

feStatus CHNS_Solver::updateBeforeFixedPointIteration(const int iter)
{
  const Parameters::CHNS::MeshAdaptation &p = _CHNS_parameters.adaptation;

  // Reduce interface thickness
  if (iter > 0 && (iter % p.epsilonDecreaseFrequency) == 0)
  {
    const double prev = _CHNS_parameters.epsilon;
    _CHNS_parameters.epsilon /= p.epsilonDecrease;
    feInfo("CHNS epsilon was updated from %+-1.6e to %+-1.6e",
      prev, _CHNS_parameters.epsilon);
  }

  // TODO: Update mobility
  // // Update functions depending on epsilon or M
  // _mobility->resetParameters({_M});

  return FE_STATUS_OK;
}

//
// Compute error on scalar fields
//
feStatus CHNS_Solver::computeError(const std::vector<feSpace *> &spaces,
                                   feMesh                       *mesh,
                                   feSolution                   *sol,
                                   const int         iIntervalReferenceSolution,
                                   TransientAdapter &adapter,
                                   std::vector<double> &errors) const
{
  if (adapter._adapt_parameters.readReferenceSolution)
  {
    feInfo("Computing error for CHNS based on reference finer solution");

    errors.resize(4, 0.);

    // feSpace *u   = spaces[0];
    feSpace *p   = spaces[1];
    feSpace *phi = spaces[2];
    feSpace *mu  = spaces[3];

    feNorm *errorOnScalarSpaces = nullptr;
    feCheckReturn(createNorm(errorOnScalarSpaces, NONE, {p, phi, mu}, sol));
    feCheckReturn(errorOnScalarSpaces->computeErrorNormFromExternalSolution(
      adapter.referenceTestCases[iIntervalReferenceSolution].mesh,
      adapter.referenceTestCases[iIntervalReferenceSolution].sol,
      adapter.referenceTestCases[iIntervalReferenceSolution].spaces,
      mesh,
      sol,
      errors));
    feInfo("Computed error on P  : %+-1.6e", errors[0]);
    feInfo("Computed error on Phi: %+-1.6e", errors[1]);
    feInfo("Computed error on Mu : %+-1.6e", errors[2]);

    delete errorOnScalarSpaces;
  }
  else
  {
    feInfo("Not computing CHNS error as no reference solution was provided.");
  }
  return FE_STATUS_OK;
}