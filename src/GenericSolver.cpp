
#include "GenericSolver.h"

#include "TransientAdapter.h"
#include "feNorm.h"
#include "feSolution.h"
#include "feSolutionContainer.h"
#include "feTimeIntegration.h"

feStatus
SolverBase::addBoundarySpace(const BoundaryConditions::BoundaryCondition *bc,
                             const FEDescriptor     &spaceDescriptor,
                             feMesh                 *mesh,
                             std::vector<feSpace *> &spaces,
                             std::vector<feSpace *> &essentialSpaces,
                             size_t                 &spaceIndex) const
{
  const FEDescriptor &d = spaceDescriptor;

  // Check if space was already added
  for(size_t i = 0; i < spaces.size(); ++i)
  {
    const auto &s = spaces[i];
    if(s->representsSameFieldAs(d))
    {
      spaceIndex = i;
      feInfo("Skipping existing space for %s - %s",
        d._fieldName.data(), d._physicalEntityName.data());
      return FE_STATUS_OK;
    }
  }

  void *initFct = d._isScalar ? (void *)d._scalarField : (void *)d._vectorField;

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
  spaceIndex = spaces.size() - 1;

  feInfo("\tAdded boundary condition:");
  feInfo("\t\tType : %s",
         BoundaryConditions::typeToString(bc->_boundaryType).data());
  feInfo("\t\tField: %s", d._fieldName.data());
  feInfo("\t\tPhysical entity: %s", d._physicalEntityName.data());

  // if(bc->_boundaryType == BoundaryConditions::Type::periodic)
  // {
  //   feInfo("\t\tMatching periodic boundary: ");
  // }

  feInfo("\t\tIs essential ?: %s", d._isEssential ? "yes" : "no");
  feInfo("\t\tEssential components: %d - %d - %d",
         d._essentialComponents[0],
         d._essentialComponents[1],
         d._essentialComponents[2]);

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

  return FE_STATUS_OK;
}

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
    // For most BC, the test and trial spaces are identical.
    // They differ for e.g. a contact angle condition
    // in the Cahn-Hilliard model, where the flux of the phase
    // marker is tested in the equation of the potential.
    const FEDescriptor &dtrial = b->_descriptor_trialSpace;
    const FEDescriptor &dtest  = b->_descriptor_testSpace;

    size_t bSpaceIndex;
    feCheckReturn(this->addBoundarySpace(b, dtrial, mesh, spaces, essentialSpaces, bSpaceIndex));

    if (!dtrial.representsSameFieldAs(dtest))
    {
      size_t unusedIndex;
      feCheckReturn(this->addBoundarySpace(b, dtest, mesh, spaces, essentialSpaces, unusedIndex));
    }

    //
    // Periodic conditions
    //
    if (b->_boundaryType == BoundaryConditions::Type::periodic)
    {
      const BoundaryConditions::ScalarPeriodic *bP =
        static_cast<const BoundaryConditions::ScalarPeriodic *>(b);
      const FEDescriptor &d0 = bP->_descriptor_trialSpace;
      const FEDescriptor &m0 = bP->_matching_descriptor;

      if (d0._fieldName != m0._fieldName)
      {
        return feErrorMsg(
          FE_STATUS_ERROR,
          "Could not create periodic BC"
          " because field descriptors have different field names: %s vs %s",
          d0._fieldName.data(),
          m0._fieldName.data());
      }

      if (d0._degree != m0._degree)
      {
        return feErrorMsg(FE_STATUS_ERROR,
                          "Could not create periodic BC"
                          " because field descriptors have different "
                          "polynomial degree: %d vs %d",
                          d0._degree,
                          m0._degree);
      }

      // Compare number of mesh vertices on matching boundaries
      int n0 = mesh->getCncGeoByName(d0._physicalEntityName)->getNumVertices();
      int n1 = mesh->getCncGeoByName(m0._physicalEntityName)->getNumVertices();

      if (n0 != n1)
      {
        return feErrorMsg(
          FE_STATUS_ERROR,
          "Could not create periodic BC because matching boundaries do "
          "not have the same number of mesh vertices: %d vs %d",
          n0, n1);
      }

      // Add matching feSpace and set both spaces as matching
      size_t b1SpaceIndex;
      feCheckReturn(this->addBoundarySpace(b, m0, mesh, spaces, essentialSpaces, b1SpaceIndex));

      feSpace *b0Space = spaces[bSpaceIndex];
      feSpace *b1Space = spaces[b1SpaceIndex];
      b0Space->setPeriodic(true);
      b1Space->setPeriodic(true);
      b0Space->setPeriodicMaster(true);
      b1Space->setPeriodicSlave(true);
      b0Space->setMatchingPeriodicSpace(b1Space);
      b1Space->setMatchingPeriodicSpace(b0Space);
      b0Space->setPeriodicOffset(bP->_offset);
      b1Space->setPeriodicOffset(bP->_offset);

      if(b0Space == b1Space || b0Space->representsSameFieldAs(*b1Space))
      {
        return feErrorMsg(FE_STATUS_ERROR, "FIXME: Space was set as matching itself...");
      }

      feInfo("Set space %s - %s matching for %s - %s",
        b0Space->getFieldID().data(), b0Space->getCncGeoID().data(),
        b1Space->getFieldID().data(), b1Space->getCncGeoID().data());

      // // Get matching boundary
      // bool found = false;
      // for (const auto &b1 : _boundaryConditions)
      // {
      //   if (b1->_boundaryType == BoundaryConditions::Type::periodic)
      //   {
      //     const BoundaryConditions::ScalarPeriodic *b1P =
      //       static_cast<const BoundaryConditions::ScalarPeriodic *>(b1);
      //     const FEDescriptor &d1 = b1P->_descriptor_trialSpace;
      //     const FEDescriptor &m1 = b1P->_matching_descriptor;

      //     if (d0._physicalEntityName == m1._physicalEntityName &&
      //         m0._physicalEntityName == d1._physicalEntityName &&
      //         d0._fieldName == m1._fieldName && m0._fieldName == d1._fieldName)
      //     {
      //       if(found)
      //       {
      //         return feErrorMsg(
      //           FE_STATUS_ERROR,
      //           "Periodic BC for field %s on entity %s has more than one "
      //           "matching boundary!", d0._fieldName.data(), d0._physicalEntityName.data());
      //       }

      //       found = true;
      //       feInfo("Matching periodic BC found");

      //       // Compare number of mesh vertices on matching boundaries
      //       int n0 =
      //         mesh->getCncGeoByName(d0._physicalEntityName)->getNumVertices();
      //       int n1 =
      //         mesh->getCncGeoByName(d1._physicalEntityName)->getNumVertices();

      //       if (n0 != n1)
      //       {
      //         return feErrorMsg(
      //           FE_STATUS_ERROR,
      //           "Could not create periodic BC because matching boundaries do "
      //           "not have the same number of mesh vertices: %d vs %d",
      //           n0, n1);
      //       }

      //       // Add matching feSpace and set both spaces as matching
      //       const FEDescriptor &d1trial = b1->_descriptor_trialSpace;

      //       size_t b1SpaceIndex;
      //       feCheckReturn(this->addBoundarySpace(b1, d1trial, mesh, spaces, essentialSpaces, b1SpaceIndex));

      //       feSpace *b0Space = spaces[bSpaceIndex];
      //       feSpace *b1Space = spaces[b1SpaceIndex];
      //       b0Space->setPeriodic(true);
      //       b1Space->setPeriodic(true);
      //       b0Space->setMatchingPeriodicSpace(b1Space);
      //       b1Space->setMatchingPeriodicSpace(b0Space);

      //       if(b0Space == b1Space || b0Space->representsSameFieldAs(*b1Space))
      //       {
      //         return feErrorMsg(FE_STATUS_ERROR, "FIXME: Space was set as matching itself...");
      //       }

      //       feInfo("Set space %s - %s matching for %s - %s",
      //         b0Space->getFieldID().data(), b0Space->getCncGeoID().data(),
      //         b1Space->getFieldID().data(), b1Space->getCncGeoID().data());

      //       break;
      //     }
      //   }
      // }

      // if (!found)
      // {
      //   return feErrorMsg(FE_STATUS_ERROR, "Did not find matching periodic BC");
      // }
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
  feCheckReturn(this->createBilinearForms(spaces, forms));

  feLinearSystem *system = nullptr;
#if defined(HAVE_MKL)
  feCheckReturn(
    createLinearSystem(system, MKLPARDISO, forms, numbering));
#elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
  feCheckReturn(
    createLinearSystem(system, PETSC_MUMPS, forms, numbering));
#else
  feCheckReturn(
    createLinearSystem(system, PETSC, forms, numbering));
#endif

  std::vector<feNorm *> norms = {};

  if (_normTypes.size() > 0)
  {
    for (size_t i = 0; i < _normTypes.size(); ++i)
    {
      const int      index = _normTypes[i].first;
      const normType nt    = _normTypes[i].second;
      if (nt == L2_ERROR)
      {
        bool found = false;
        for (const auto &p : _scalarExactSolutions)
        {
          if (p.first == index)
          {
            found                           = true;
            const feFunction *exactSolution = p.second;
            if (exactSolution)
            {
              feNorm *n;
              feCheckReturn(
                createNorm(n, L2_ERROR, {spaces[index]}, sol, exactSolution));
              norms.push_back(n);
            }
          }
        }
        if (!found)
        {
          return feErrorMsg(
            FE_STATUS_ERROR,
            "An L2_ERROR norm is set to be computed for FE space %d,"
            " but no scalar reference field was provided for this space.",
            index);
        }
      }

      else if (nt == H1_ERROR || nt == SEMI_H1_ERROR || nt == VECTOR_L2_ERROR ||
               nt == VECTOR_SEMI_H1_ERROR)
      {
        bool found = false;
        for (const auto &p : _vectorExactSolutions)
        {
          if (p.first == index)
          {
            found                                 = true;
            const feVectorFunction *exactSolution = p.second;
            if (exactSolution)
            {
              feNorm *n;
              feCheckReturn(createNorm(
                n, nt, {spaces[index]}, sol, nullptr, exactSolution));
              norms.push_back(n);
            }
          }
        }
        if (!found)
        {
          return feErrorMsg(
            FE_STATUS_ERROR,
            "An H1_ERROR/VECTOR_L2_ERROR/VECTOR_SEMI_H1_ERROR "
            " norm is set to be computed for FE space %d,"
            " but no vector reference field was provided for this space.",
            index);
        }
      }

      else
      {
        return feErrorMsg(FE_STATUS_ERROR,
                          "Prescribed norm in GenericSolver is not handled.");
      }
    }
  }

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

  _postProcessingData = timeIntegrator->getPostProcessingData();

  delete timeIntegrator;
  delete system;
  for (feBilinearForm *f : forms)
    delete f;
  for (feNorm *norm : norms)
    delete norm;

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
  feCheckReturn(this->createBilinearForms(spaces, forms));

  //
  // Linear system, nonlinear solver and time integration
  //
  feLinearSystem *system = nullptr;
#if defined(HAVE_MKL)
  feCheckReturn(
    createLinearSystem(system, MKLPARDISO, forms, numbering));
#elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
  feCheckReturn(
    createLinearSystem(system, PETSC_MUMPS, forms, numbering));
#else
  feCheckReturn(
    createLinearSystem(system, PETSC, forms, numbering));
#endif

  std::vector<feNorm *> norms = {};

  ///////////////////////////////////////
  feNorm *norm1;
  if (_scalarExactSolutions.size() > 0)
  {
    const feFunction *exactSolution = _scalarExactSolutions[0].second;
    if (exactSolution)
    {
      feCheckReturn(
        createNorm(norm1, L2_ERROR, {spaces[0]}, sol, exactSolution));
      norms.push_back(norm1);
    }
  }
  ///////////////////////////////////////

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

  if (adapter._adapt_parameters.adapt)
  {
    // Need to compute metrics at t = t_i (beginning of time interval) as well
    if (adapter._adapt_parameters.useExactDerivatives)
    {
      feMetric *metricField = adapter.allMetrics[iInterval];

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
      // metricField->setMetricsToIdentity();

      // Increment integral for t = ti
      metricField->addMetricsToOther(dt / 2., adapter.allHi[iInterval]);
    }
    else
    {
      feMetric *metricField = adapter.allMetrics[iInterval];

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

      // Increment integral for t = ti
      metricField->addMetricsToOther(dt / 2., adapter.allHi[iInterval]);
    }
  }

  // Solve and compute metrics at each time step
  for (int iStep = 0; iStep < nT; ++iStep)
  {
    tic();
    feCheckReturn(timeIntegrator->makeSteps(1));
    feInfo("Computed time step in %f s", toc());

    if (adapter._adapt_parameters.adapt)
    {
      //
      // Compute metric field
      //
      feMetric *metricField = adapter.allMetrics[iInterval];

      if (adapter._adapt_parameters.useExactDerivatives)
      {
        // Use analytical field derivatives
        // Check that the derivatives were given
        const int order = adapter._adapt_parameters.orderForAdaptation;
        if ((order == 1 &&
             metricField->_options.secondDerivatives == nullptr) ||
            (order == 2 && metricField->_options.thirdDerivatives == nullptr) ||
            (order == 3 &&
             metricField->_options.fourthDerivatives == nullptr) ||
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
        // metricField->setMetricsToIdentity();
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
      if (iStep == nT - 1)
      {
        feInfo("With dt/2 at t = %f", timeIntegrator->getCurrentTime());
        metricField->addMetricsToOther(dt / 2., adapter.allHi[iInterval]);
      }
      else
      {
        feInfo("With dt   at t = %f", timeIntegrator->getCurrentTime());
        metricField->addMetricsToOther(dt, adapter.allHi[iInterval]);
      }
    }

    adapter.currentTimeStep++;
  }

  adapter.currentTime = timeIntegrator->getCurrentTime();

  feInfo("Time integrator computed %d norms", norms.size());

  std::vector<std::vector<double>> postProc =
    timeIntegrator->getPostProcessingData();
  if (postProc.size() >= 2)
  {
    for (int j = 0; j < nT + 1; ++j)
    {
      feInfo("Error at step %2d (t = %1.3e) = %1.6e",
             j,
             postProc[0][j],
             postProc[1][j]);
      adapter.allErrorsAllTimeSteps[adapter.iFixedPoint][0].push_back(
        postProc[0][j]);
      adapter.allErrorsAllTimeSteps[adapter.iFixedPoint][1].push_back(
        postProc[1][j]);
    }
  }

  delete adapter.allContainers[iInterval];
  adapter.allContainers[iInterval] = new feSolutionContainer();
  adapter.allContainers[iInterval]->NaNify();
  *(adapter.allContainers[iInterval]) = timeIntegrator->getSolutionContainer();

  delete timeIntegrator;
  delete system;
  for (feBilinearForm *f : forms)
    delete f;
  for (feNorm *norm : norms)
    delete norm;

  feInfo("Done solving on this interval");

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
  feInfo("Starting with mesh %p", mesh);

  mesh = new feMesh2DP1(meshName);

  feInfo("Created mesh %p", mesh);

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