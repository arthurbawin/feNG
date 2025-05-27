
#include "GenericSolver.h"

#include "TransientAdapter.h"
#include "feNorm.h"
#include "feSolution.h"
#include "feSolutionContainer.h"
#include "feTimeIntegration.h"

feStatus SolverBase::createSpaces(feMesh                 *mesh,
                                  std::vector<feSpace *> &spaces,
                                  std::vector<feSpace *> &essentialSpaces) const
{
  spaces.clear();
  essentialSpaces.clear();

  //
  // Unknown fields
  //
  for (const auto &d : _fieldDescriptors)
  {
    void *initFct =
      d._isScalar ? (void *)d._scalarField : (void *)d._vectorField;

    feSpace *s;
    feCheckReturn(createFiniteElementSpace(s,
                                           mesh,
                                           d._feType,
                                           d._degree,
                                           d._fieldName,
                                           d._physicalEntityName,
                                           d._degreeQuadrature,
                                           initFct));
    spaces.push_back(s);
  }

  //
  // Boundary conditions
  //
  for (const auto &b : _boundaryConditions)
  {
    const FEDescriptor &d = b._descriptor;

    void *initFct =
      d._isScalar ? (void *)d._scalarField : (void *)d._vectorField;

    feSpace *s;
    feCheckReturn(createFiniteElementSpace(s,
                                           mesh,
                                           d._feType,
                                           d._degree,
                                           d._fieldName,
                                           d._physicalEntityName,
                                           d._degreeQuadrature,
                                           initFct));
    spaces.push_back(s);

    // Print a summary of the boundary condition
    feInfo("");

    //
    // Assign essential components
    //
    if (d._isEssential)
    {
      essentialSpaces.push_back(s);
    }
    else
    {
      for (int i = 0; i < 3; ++i)
      {
        if (d._essentialComponents[i] == 1)
        {
          s->setEssentialComponent(i, true);
        }
      }
    }
  }

  return FE_STATUS_OK;
}

feStatus SolverBase::solve(feMesh                 *mesh,
                           feSolution             *sol,
                           feMetaNumber           *numbering,
                           std::vector<feSpace *> &spaces,
                           feExportData           &exportData) const
{
  //
  // Time stepping is controlled by the TimeIntegration parameters
  //
  const Parameters::TimeIntegration &tp = _TimeIntegration_parameters;
  const int                          nT = tp.nTimeSteps;
  const double                       t0 = tp.t_initial;
  const double                       tf = tp.t_final;

  //
  // Linear and bilinear forms
  //
  std::vector<feBilinearForm *> forms;
  this->createBilinearForms(spaces, forms);

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

  //
  // Until NLoptions and timeIntegratorScheme are fully replaced by feParameters
  // structures:
  //
  feNLSolverOptions NLoptions{
    _NLSolver_parameters.tolResidual,
    _NLSolver_parameters.tolCorrection,
    _NLSolver_parameters.tolDivergence,
    _NLSolver_parameters.maxIter,
    _NLSolver_parameters.recomputeJacobianEveryNsteps,
    _NLSolver_parameters.residualDecrease,
  };

  timeIntegratorScheme scheme;

  switch (tp.method)
  {
    case Parameters::TimeIntegration::TimeIntegrationMethod::stationary:
      scheme = timeIntegratorScheme::STATIONARY;
      break;
    case Parameters::TimeIntegration::TimeIntegrationMethod::bdf1:
      scheme = timeIntegratorScheme::BDF1;
      break;
    case Parameters::TimeIntegration::TimeIntegrationMethod::bdf2:
      scheme = timeIntegratorScheme::BDF2;
      break;
    default:
      return feErrorMsg(FE_STATUS_ERROR,
                        "Time integration method not handled.");
  }

  TimeIntegrator *timeIntegrator = nullptr;
  feCheckReturn(createTimeIntegrator(timeIntegrator,
                                     scheme,
                                     NLoptions,
                                     system,
                                     sol,
                                     mesh,
                                     norms,
                                     exportData,
                                     t0,
                                     tf,
                                     nT));

  // Set BDF2 starter if needed
  if (tp.method == Parameters::TimeIntegration::TimeIntegrationMethod::bdf2)
  {
    if (tp.bdf2starter == Parameters::TimeIntegration::BDF2Starter::bdf1)
    {
      static_cast<BDF2Integrator *>(timeIntegrator)
        ->setStartingMethod(BDF2Starter::BDF1);
    }
    if (tp.bdf2starter ==
        Parameters::TimeIntegration::BDF2Starter::exactSolution)
    {
      static_cast<BDF2Integrator *>(timeIntegrator)
        ->setStartingMethod(BDF2Starter::InitialCondition);
    }
  }

  // Solve
  feCheckReturn(timeIntegrator->makeSteps(nT));

  delete timeIntegrator;
  delete system;
  for (feBilinearForm *f : forms)
    delete f;

  return FE_STATUS_OK;
}

feStatus SolverBase::solve(feMesh                 *mesh,
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
  const Parameters::TimeIntegration &tp = _TimeIntegration_parameters;
  const int    nT = adapter._adapt_parameters.nTimeStepsPerInterval;
  const double t0 = adapter.currentTime;
  const double tf = t0 + adapter._time_parameters.dt * nT;
  const double dt = adapter._time_parameters.dt;

  //
  // Linear and bilinear forms
  //
  std::vector<feBilinearForm *> forms;
  this->createBilinearForms(spaces, forms);

  //
  // Linear system, nonlinear solver and time integration
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

  //
  // Until NLoptions and timeIntegratorScheme are fully replaced by feParameters
  // structures:
  //
  feNLSolverOptions NLoptions{
    _NLSolver_parameters.tolResidual,
    _NLSolver_parameters.tolCorrection,
    _NLSolver_parameters.tolDivergence,
    _NLSolver_parameters.maxIter,
    _NLSolver_parameters.recomputeJacobianEveryNsteps,
    _NLSolver_parameters.residualDecrease,
  };

  timeIntegratorScheme scheme;

  switch (tp.method)
  {
    case Parameters::TimeIntegration::TimeIntegrationMethod::stationary:
      scheme = timeIntegratorScheme::STATIONARY;
      break;
    case Parameters::TimeIntegration::TimeIntegrationMethod::bdf1:
      scheme = timeIntegratorScheme::BDF1;
      break;
    case Parameters::TimeIntegration::TimeIntegrationMethod::bdf2:
      scheme = timeIntegratorScheme::BDF2;
      break;
    default:
      return feErrorMsg(FE_STATUS_ERROR,
                        "Time integration method not handled.");
  }

  TimeIntegrator *timeIntegrator = nullptr;
  feCheckReturn(createTimeIntegrator(timeIntegrator,
                                     scheme,
                                     NLoptions,
                                     system,
                                     sol,
                                     mesh,
                                     norms,
                                     exportData,
                                     t0,
                                     tf,
                                     nT));

  // Set BDF2 starter if needed
  if (tp.method == Parameters::TimeIntegration::TimeIntegrationMethod::bdf2)
  {
    if (tp.bdf2starter == Parameters::TimeIntegration::BDF2Starter::bdf1)
    {
      static_cast<BDF2Integrator *>(timeIntegrator)
        ->setStartingMethod(BDF2Starter::BDF1);
    }
    if (tp.bdf2starter ==
        Parameters::TimeIntegration::BDF2Starter::exactSolution)
    {
      static_cast<BDF2Integrator *>(timeIntegrator)
        ->setStartingMethod(BDF2Starter::InitialCondition);
    }
  }

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

    //
    // Compute metric field
    //
    feMetric *metricField = adapter.allMetrics[iInterval];

    if (adapter._adapt_parameters.useExactDerivatives)
    {
      // Use analytical field derivatives
      // Check that the derivatives were given
      const int order = adapter._adapt_parameters.orderForAdaptation;
      if ((order == 1 && metricField->_options.secondDerivatives == nullptr) ||
          (order == 2 && metricField->_options.thirdDerivatives == nullptr) ||
          (order == 3 && metricField->_options.fourthDerivatives == nullptr) ||
          (order == 4 && metricField->_options.fifthDerivatives == nullptr))
      {
        std::cout << metricField->_options.secondDerivatives << std::endl;
        return feErrorMsg(
          FE_STATUS_ERROR,
          "The metric field is set to be computed with exact field "
          "derivatives,\n but no callback for derivatives of order k+1 was "
          "provided for field of order k = %d driving the adaptation.\n",
          order);
      }
      metricField->setCurrentTime(timeIntegrator->getCurrentTime());
      metricField->setMetricScaling(false);
      feCheckReturn(metricField->computeMetrics());
    }
    else
    {
      // Reconstruct prescribed field for adaptation
      const int   id = adapter._adapt_parameters.spaceIDForAdaptation;
      feSpace    *spaceToReconstruct = spaces[id];
      std::string meshFile = adapter.allOptions[iInterval].backgroundMeshfile;
      std::string recoveryFile =
        adapter._io_parameters.writeDir + "reconstruction.msh";
      feNewRecovery recoveredField(spaceToReconstruct,
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

      metricField->setRecoveredFields({&recoveredField});
      metricField->setCurrentTime(timeIntegrator->getCurrentTime());
      metricField->setMetricScaling(false);
      feCheckReturn(metricField->computeMetrics());
    }

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

feStatus
SolverBase::projectSolution(feMesh2DP1                  *currentMesh,
                            feMesh2DP1                  *nextMesh,
                            feSolutionContainer         *container,
                            feSolutionContainer         &bufferContainer,
                            feSolutionContainer        *&nextContainer,
                            feMetaNumber                *numbering,
                            const std::vector<feSpace *> spaces,
                            const std::vector<feSpace *> essentialSpaces) const
{
  // Copy current container into the buffer
  bufferContainer = *container;

  // For checks : last BDF solutions
  std::vector<feSolution> currentSolutions(bufferContainer.getNbSol());
  for (int i = 0; i < bufferContainer.getNbSol(); ++i)
  {
    currentSolutions[i] =
      feSolution(bufferContainer, i, spaces, essentialSpaces);
  }

  // Projection
  std::vector<feSpace *> nextSpaces, nextEssentialSpaces;
  feCheckReturn(this->createSpaces(nextMesh, nextSpaces, nextEssentialSpaces));
  feMetaNumber nextNumbering(nextMesh, nextSpaces, nextEssentialSpaces);
  feCheckReturn(currentMesh->transfer(nextMesh,
                                      numbering,
                                      &nextNumbering,
                                      &bufferContainer,
                                      spaces,
                                      essentialSpaces,
                                      nextSpaces));

  delete nextContainer;
  nextContainer = new feSolutionContainer(bufferContainer);

  // Check projection (visualization) and integral of each projected solution
  // in container
  for (int i = 0; i < bufferContainer.getNbSol(); ++i)
  {
    feInfo("Checking solution %d in bufferContainer :", i);
    feSolution projectedSolution(bufferContainer, i, spaces, essentialSpaces);
    // feExporter *nextExporter;
    // feCheck(createVisualizationExporter(nextExporter, VTK, &nextNumbering,
    // &projectedSolution, nextMesh, nextSpaces));
    // feCheck(nextExporter->writeStep("afterTransfer" + std::to_string(i) +
    // ".vtk")); delete nextExporter;

    // Integral of highest dimensionsal fields before and after transfer
    for (size_t j = 0; j < spaces.size(); ++j)
    {
      if (spaces[j]->getDim() == 2)
      {
        feNorm *integral, *integralNextMesh, *intDot, *intDotNextMesh;
        feCheck(
          createNorm(integral, INTEGRAL, {spaces[j]}, &currentSolutions[i]));
        feCheck(createNorm(
          integralNextMesh, INTEGRAL, {nextSpaces[j]}, &projectedSolution));
        feInfo("Integral of field %s on current mesh = %+-1.6e",
               spaces[j]->getFieldID().data(),
               integral->compute());
        feInfo("Integral of field %s on next    mesh = %+-1.6e",
               nextSpaces[j]->getFieldID().data(),
               integralNextMesh->compute());
        feCheck(
          createNorm(intDot, INTEGRAL_DT, {spaces[j]}, &currentSolutions[i]));
        feCheck(createNorm(
          intDotNextMesh, INTEGRAL_DT, {nextSpaces[j]}, &projectedSolution));
        feInfo("Integral of field %s dot on current mesh = %+-1.6e",
               spaces[j]->getFieldID().data(),
               intDot->compute());
        feInfo("Integral of field %s dot on next    mesh = %+-1.6e",
               nextSpaces[j]->getFieldID().data(),
               intDotNextMesh->compute());
        delete integral;
        delete integralNextMesh;
        delete intDot;
        delete intDotNextMesh;
      }
    }
  }

  for (feSpace *s : nextSpaces)
    delete s;

  return FE_STATUS_OK;
}

feStatus SolverBase::readReferenceTestCase(const std::string meshName,
                                           const std::string solutionFileName,
                                           feMesh          *&mesh,
                                           std::vector<feSpace *> &spaces,
                                           feMetaNumber          *&numbering,
                                           feSolution            *&sol) const
{
  mesh = new feMesh2DP1(meshName);

  // Essential spaces are not needed for error computation, and are not
  // returned
  std::vector<feSpace *> refEssentialSpaces;
  this->createSpaces(mesh, spaces, refEssentialSpaces);
  numbering = new feMetaNumber(mesh, spaces, refEssentialSpaces);

  std::filebuf fb;
  if (fb.open(solutionFileName, std::ios::in))
  {
    fb.close();

    // Read solution file if it exists
    sol = new feSolution(numbering->getNbDOFs(), solutionFileName);
  }
  else
  {
    return feErrorMsg(FE_STATUS_ERROR,
                      "Could not create reference test case because file "
                      "\"%s\" does not exist.",
                      solutionFileName.data());
  }
  return FE_STATUS_OK;
}