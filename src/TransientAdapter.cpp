
#include "TransientAdapter.h"

#include <numeric>

#include "feOptionsParser.h"

#if defined(HAVE_GMSH)
#  include "gmsh.h"
#endif

TransientAdapter::TransientAdapter(
  const Parameters::MeshAdaptation  &adaptParameters,
  const Parameters::TimeIntegration &timeParameters,
  const Parameters::IO              &IOParameters)
  : _adapt_parameters(adaptParameters)
  , _time_parameters(timeParameters)
  , _io_parameters(IOParameters)
{
  const int    nI     = adaptParameters.nIntervals;
  const int    nT     = adaptParameters.nTimeStepsPerInterval;
  const double t0     = timeParameters.t_initial;
  const double tf     = timeParameters.t_final;
  currentTime         = t0;
  _time_parameters.dt = (tf - t0) / (nI * nT);

  feInfo("TRANSIENT ADAPTER:");
  feInfo("\tInitial  time               : %1.6e", t0);
  feInfo("\tFinal    time               : %1.6e", tf);
  feInfo("\t# of subintervals           : %d", nI);
  feInfo("\t# time steps per subinterval: %d", nT);
  feInfo("\tTotal number of time steps  : %d", nI * nT);
  feInfo("\tTime step dt                : %1.6e", _time_parameters.dt);
}

feStatus TransientAdapter::allocate(const SolverBase &solver)
{
  const Parameters::MeshAdaptation &ap = _adapt_parameters;
  const Parameters::IO             &io = _io_parameters;
  const int                         nI = ap.nIntervals;

  // Create problem structures for each time interval
  allMeshes.resize(nI, nullptr);
  allMetrics.resize(nI, nullptr);
  allHi.resize(nI);
  allContainers.resize(nI, nullptr);
  allNumberings.resize(nI, nullptr);
  allSpaces.resize(nI);
  allEssentialSpaces.resize(nI);

  allErrors.resize(ap.nFixedPoint, std::vector<std::vector<double>>(nI));
  const int nUnknownFields = solver.getNumFields();
  for (int i = 0; i < ap.nFixedPoint; ++i)
    for (int j = 0; j < nI; ++j)
      allErrors[i][j].resize(nUnknownFields, 0.);

  spaceTimeComplexity.resize(ap.nFixedPoint);

  referenceTestCases.resize(ap.nIntervalsReferenceSolution);

  // Metric options common to all meshes
  feMetricOptions metricOptions(io.initialMesh);
  metricOptions.polynomialDegree = ap.orderForAdaptation;

  // Check that the provided polynomial order for mesh adaptation
  // matches the one of the computed field in the solver.
  // This would still work with mismatched degrees, but would
  // yield non-optimal adapted meshes.
  const FEDescriptor &fieldForAdaptation =
    solver._fieldDescriptors[ap.spaceIDForAdaptation];
  const int solvedDegree = fieldForAdaptation._degree;
  if (solvedDegree != ap.orderForAdaptation)
  {
    return feErrorMsg(
      FE_STATUS_ERROR,
      "The prescribed field degree for mesh adaptation (%d) does not match the "
      "degree of this field in the solver (%d).\n"
      "Field ID given to drive the mesh adaptation procedure: %d\n"
      "Associated field name   in solver: %s\n"
      "Associated field entity in solver: %s\n"
      "Associated field degree in solver: %d\n",
      ap.orderForAdaptation,
      solvedDegree,
      ap.spaceIDForAdaptation,
      fieldForAdaptation._fieldName.data(),
      fieldForAdaptation._physicalEntityName.data(),
      fieldForAdaptation._degree);
  }

  switch (ap.orderForAdaptation)
  {
    case 1:
      metricOptions.method = adaptationMethod::ANISO_P1;
      break;
    case 2:
      metricOptions.method = adaptationMethod::ANISO_P2;
      break;
    default:
      metricOptions.method = adaptationMethod::ANISO_PN;
  }

  metricOptions.nTargetVertices   = ap.targetVertices;
  metricOptions.LpNorm            = ap.Lp_norm;
  metricOptions.hMin              = 1e-10;
  metricOptions.hMax              = 10.;
  metricOptions.enableGradation   = ap.enableGradation;
  metricOptions.gradation         = ap.gradation;
  metricOptions.gradationSpace    = GradationSpace::Mixed;
  metricOptions.mmgInputMeshfile  = io.writeDir + "inputMMG.msh";
  metricOptions.mmgOutputMeshfile = io.writeDir + "outputMMG.msh";
  metricOptions.debug             = false;

  if (ap.useExactDerivatives)
  {
    metricOptions.useAnalyticDerivatives = true;
  }

  allOptions.resize(nI, metricOptions);

  //
  // Read initial meshes and create metric fields
  //
  for (int iI = 0; iI < nI; ++iI)
  {
    // Different adapted mesh name for each sub-interval
    allOptions[iI].adaptedMeshName =
      io.writeDir + "adapted_ifp_0_interval_" + std::to_string(iI) + ".msh";

    if (ap.restartFromAdaptedMeshes)
    {
      std::string meshFile = io.writeDir + ap.restartMeshesRoot + "_interval" +
                             std::to_string(iI) + ".msh";
      allMeshes[iI] = new feMesh2DP1(meshFile);
    }
    else
    {
      allMeshes[iI] = new feMesh2DP1(io.initialMesh);
    }

    feCheckReturn(
      createMetricField(allMetrics[iI], allOptions[iI], allMeshes[iI]));
    allMetrics[iI]->setMetricsToZero();
    allMetrics[iI]->getMetrics(allHi[iI]);
    allMetrics[iI]->setMetricsToIdentity();
  }

  // Read reference solutions for error computations
  feInfo("TRANSIENT ADAPTER:");
  feInfo("\tReading reference solution on each interval");
  if (ap.readReferenceSolution)
  {
    // Number of subintervals should be a divider
    // of the number of reference subintervals
    if (ap.nIntervalsReferenceSolution % nI != 0)
    {
      return feErrorMsg(
        FE_STATUS_ERROR,
        "The number of sub-intervals (%d) for this simulation is not a divider "
        "of the number of sub-intervals for the reference solution (%d).",
        ap.nIntervalsReferenceSolution,
        nI);
    }

    for (int iIref = 0; iIref < ap.nIntervalsReferenceSolution; ++iIref)
    {
      std::string refMeshFile = io.writeDir + ap.referenceMeshRoot +
                                "_interval" + std::to_string(iIref) + ".msh";
      std::string refSolutionFile = io.writeDir + ap.referenceSolutionRoot +
                                    "_interval" + std::to_string(iIref) +
                                    ".txt";

      feMesh *mesh = referenceTestCases[iIref].mesh;

      solver.readReferenceTestCase(refMeshFile,
                                   refSolutionFile,
                                   mesh,
                                   referenceTestCases[iIref].spaces,
                                   referenceTestCases[iIref].numbering,
                                   referenceTestCases[iIref].sol);
    }
  }

  // Open post-processing file
  postProcFile.open(io.writeDir + "postProcessing_iConv" + std::to_string(0) +
                    ".txt");

  wasAllocated = true;

  return FE_STATUS_OK;
}

feStatus
TransientAdapter::setExactDerivatives(const int               order,
                                      const feVectorFunction *derivatives)
{
  if (!wasAllocated)
  {
    return feErrorMsg(
      FE_STATUS_ERROR,
      "Cannot set exact derivatives for metric computation before "
      "TransientAdapter is allocated. Call allocate(...) function first.");
  }
  for (int iI = 0; iI < _adapt_parameters.nIntervals; ++iI)
  {
    // Set callback for the current metric fields (already created)
    // and for the options, which are used to create the subsequent
    // metric fields, at the next fixed point iteration.
    switch (order)
    {
      case 1:
        allOptions[iI].firstDerivatives           = derivatives;
        allMetrics[iI]->_options.firstDerivatives = derivatives;
        break;
      case 2:
        allOptions[iI].secondDerivatives           = derivatives;
        allMetrics[iI]->_options.secondDerivatives = derivatives;
        break;
      case 3:
        allOptions[iI].thirdDerivatives           = derivatives;
        allMetrics[iI]->_options.thirdDerivatives = derivatives;
        break;
      case 4:
        allOptions[iI].fourthDerivatives           = derivatives;
        allMetrics[iI]->_options.fourthDerivatives = derivatives;
        break;
      case 5:
        allOptions[iI].fifthDerivatives           = derivatives;
        allMetrics[iI]->_options.fifthDerivatives = derivatives;
        break;
      default:
        return feErrorMsg(FE_STATUS_ERROR,
                          "Derivatives of order %d are not handled.",
                          order);
    }
  }
  return FE_STATUS_OK;
}

feStatus TransientAdapter::readAdaptedMeshes()
{
  for (int iI = 0; iI < _adapt_parameters.nIntervals; ++iI)
  {
    // Read adapted mesh
    feInfo("Deleting mesh %d", iI);
    delete allMeshes[iI];
    std::string meshFile = _io_parameters.writeDir + "adapted_ifp_" +
                           std::to_string(iFixedPoint - 1) + "_interval_" +
                           std::to_string(iI) + ".msh";
    allMeshes[iI] = new feMesh2DP1(meshFile);

    // Create new metric field
    allOptions[iI].backgroundMeshfile = meshFile;
    allOptions[iI].adaptedMeshName = _io_parameters.writeDir + "adapted_ifp_" +
                                     std::to_string(iFixedPoint) +
                                     "_interval_" + std::to_string(iI) + ".msh";
    delete allMetrics[iI];
    allMetrics[iI] = nullptr;
    feCheckReturn(
      createMetricField(allMetrics[iI], allOptions[iI], allMeshes[iI]));
    allHi[iI].clear();
    allMetrics[iI]->setMetricsToZero();
    allMetrics[iI]->getMetrics(allHi[iI]);
    allMetrics[iI]->setMetricsToIdentity();
  }

  return FE_STATUS_OK;
}

void TransientAdapter::updateBeforeFixedPointIteration()
{
  Parameters::MeshAdaptation &param = _adapt_parameters;
  if (iFixedPoint > 0)
  {
    //
    // Increase Lp norm after first iteration
    //
    param.Lp_norm = 100.;
    for (int iI = 0; iI < param.nIntervals; ++iI)
    {
      allOptions[iI].LpNorm = param.Lp_norm;
    }

    //
    // Modify target number of vertices
    //
    if (iFixedPoint % param.verticesIncreaseFrequency == 0)
    {
      const double prev = param.targetVertices;
      param.targetVertices *= param.targetVerticesIncrease;
      for (int iI = 0; iI < param.nIntervals; ++iI)
      {
        allOptions[iI].nTargetVertices = param.targetVertices;
      }
      feInfo("Transient adapter: target number of vertices was updated from %f "
             "to %f",
             prev,
             param.targetVertices);
    }
  }
}

feStatus TransientAdapter::computeSpaceTimeMetrics()
{
#if defined(HAVE_GMSH)
  const Parameters::MeshAdaptation &ap = _adapt_parameters;
  const int                         nI = ap.nIntervals;
  //
  // Compute metric scaling
  //
  // Space-time complexity from Nv :
  // const double targetNst = (double) adapter.nVertices * adapter.nIntervals *
  // adapter.nTimeStepsPerInterval;
  const double targetNst            = (double)ap.targetVertices * nI;
  const double dim                  = ap.dim;
  const double p                    = ap.Lp_norm;
  const double nT                   = (double)ap.nTimeStepsPerInterval;
  const double expDeterminantGlobal = p / (2. * p + dim);
  const double expDeterminantLocal  = -1. / (2. * p + dim);

  double sumKj = 0.;
  for (int iI = 0; iI < nI; ++iI)
  {
    feMetric *metricField = allMetrics[iI];
    metricField->setMetrics(allHi[iI]);
    const double Kj =
      metricField->computeIntegralOfDeterminant(expDeterminantGlobal);
    // Use simplification for constant time step. Otherwise, multiply Kj
    // by integral of (time step)^-1 on sub-interval before summing.
    sumKj += Kj;
  }

  // Integral of (time step)^-1 is nT, exponents of nT reduce to 2/d
  const double globalScalingFactor = pow(targetNst / (nT * sumKj), 2. / dim);

  std::vector<int> numVerticesInAdaptedMeshes(nI, 0);

  //
  // Apply scaling and gradation
  //
  // for(int iI = 0; iI < nI; ++iI)
  // {
  //   feMetric *metricField = allMetrics[iI];
  //   feInfo("Scaling metrics");
  //   metricField->scaleMetricsByDeterminant(globalScalingFactor,
  //   expDeterminantLocal); feInfo("Applying gradation");
  //   metricField->applyGradation();
  // }
  for (int iI = 0; iI < nI; ++iI)
  {
    feMetric *metricField = allMetrics[iI];
    metricField->scaleMetricsByDeterminant(globalScalingFactor,
                                           expDeterminantLocal);

    initializeGmsh();
    feMetricOptions my_options = metricField->getOptions();
    gmsh::open(my_options.backgroundMeshfile);
    // openGmshModel(my_options.backgroundMeshfile);
    metricField->setDebugPrefix(std::to_string(iI));

    if (ap.enableSpaceTimeGradation)
    {
      metricField->applyGradation();
    }

    metricField->writeMetricField();
    setVerbose(0);
    allMeshes[iI]->adapt(metricField);
    setVerbose(1);

    // Read the adapted mesh to check the number of vertices
    feMesh2DP1 mesh(my_options.adaptedMeshName);
    numVerticesInAdaptedMeshes[iI] = mesh.getNumVertices();
  }

  feInfo("Number of vertices in adapted meshes :");
  for (int iI = 0; iI < nI; ++iI)
    feInfo("Mesh %d : %d", iI, numVerticesInAdaptedMeshes[iI]);
  int totalVertices = std::accumulate(numVerticesInAdaptedMeshes.begin(),
                                      numVerticesInAdaptedMeshes.end(),
                                      0);

  const int effectiveNst = totalVertices * ap.nTimeStepsPerInterval;

  spaceTimeComplexity[iFixedPoint] = effectiveNst;

  feInfo("Total  number of space      vertices : %d", totalVertices);
  feInfo("Total  number of space-time vertices : %d", effectiveNst);
  feInfo("Target number of space-time vertices : %f", targetNst);
#else
  feWarning("Gmsh is required to compute and write metric fields.");
#endif
  return FE_STATUS_OK;
}

feStatus TransientAdapter::fixedPointIteration(SolverBase &solver)
{
  const Parameters::MeshAdaptation &ap = _adapt_parameters;
  const int                         nI = ap.nIntervals;
  //
  // Updates in adapter and/or solver
  //
  this->updateBeforeFixedPointIteration();
  solver.updateBeforeFixedPointIteration(iFixedPoint);

  // Reset time and number of time steps
  currentTime     = _time_parameters.t_initial;
  currentTimeStep = 0;

  // Reinitialize integrals of Hessians
  for (int iI = 0; iI < nI; ++iI)
    for (auto &integralHessian : this->allHi[iI])
      integralHessian.second.setToZero();

  // Read all adapted meshes if not at the first fixed-point iteration
  if (iFixedPoint > 0)
    feCheckReturn(this->readAdaptedMeshes());

  PostProcData postProcData;

  //
  // Solve on each subinterval
  // Compute metric every prescribed number of time steps
  //
  for (int iI = 0; iI < nI; ++iI)
  {
    feInfo("TRANSIENT ADAPTER:");
    feInfo("\tGlobal fixed-point iteration %d/%d",
           iFixedPoint + 1,
           ap.nFixedPoint);
    feInfo("\tSub-interval %d/%d", iI + 1, nI);

    if (iI > 0)
    {
      //
      // Transfer solution to the mesh of the next sub-interval
      //
      feCheck(solver.projectSolution(allMeshes[iI - 1],
                                     allMeshes[iI],
                                     allContainers[iI - 1],
                                     bufferContainer,
                                     allContainers[iI],
                                     allNumberings[iI - 1],
                                     allSpaces[iI - 1],
                                     allEssentialSpaces[iI - 1]));

      feInfo("");
      feInfo("DOF in source container : %d",
             allContainers[iI - 1]->getNbDOFs());
      feInfo("DOF in buffer container : %d", bufferContainer.getNbDOFs());
      if (allContainers[iI])
        feInfo("DOF in target container : %d", allContainers[iI]->getNbDOFs());
      else
        feInfo("DOF in target container is null");
      feInfo("");
    }

    //
    // Solve on sub-interval and compute metrics
    //
    feMesh                 *mesh            = allMeshes[iI];
    std::vector<feSpace *> &spaces          = allSpaces[iI];
    std::vector<feSpace *> &essentialSpaces = allEssentialSpaces[iI];
    feCheck(solver.createSpaces(mesh, spaces, essentialSpaces));
    feMetaNumber *numbering = new feMetaNumber(mesh, spaces, essentialSpaces);
    if (iFixedPoint > 0)
      delete allNumberings[iI];
    allNumberings[iI] = numbering;
    feSolution sol(numbering->getNbDOFs(), spaces, essentialSpaces);

    // Start from projected solution if not first sub-interval
    if (iI > 0)
    {
      sol.setSolFromContainer(allContainers[iI], 0);
      for (auto *s : spaces)
        s->setDOFInitialization(dofInitialization::PREVIOUS_SOL);
      for (auto *s : essentialSpaces)
        s->setDOFInitialization(dofInitialization::NODEWISE);
    }

    feExporter *exporter = nullptr;
    std::string vtkFileRoot =
      _io_parameters.writeDir + "sol_iConv" + std::to_string(0);
    if (_io_parameters.exportVisualizationFile)
    {
      feCheck(createVisualizationExporter(
        exporter, VTK, numbering, &sol, mesh, spaces));
    }
    feExportData exportData = {exporter,
                               _io_parameters.exportFrequency,
                               vtkFileRoot};

    //
    // Solve. The time stepping information is provided by *this.
    //
    feCheck(solver.solve(mesh, &sol, numbering, spaces, exportData, *this, iI));

    //
    // Compute error on this sub-interval
    //
    const int refToCurrentRatio = ap.nIntervalsReferenceSolution / nI;
    const int iIref             = ((iI + 1) * refToCurrentRatio) - 1;

    feInfo("Computing error:");
    if (ap.readReferenceSolution)
    {
      feInfo("If reference solution is used, comparing:");
      feInfo("Reference solution at end of interval: %d (end time = %1.4e)",
             iIref,
             this->referenceTestCases[iIref].sol->getCurrentTime());
      feInfo("Current   solution at end of interval: %d (end time = %1.4e)",
             iI,
             sol.getCurrentTime());
    }

    feCheck(solver.computeError(
      spaces, mesh, &sol, iIref, *this, allErrors[iFixedPoint][iI]));

    //
    // Post-process
    //
    feCheck(solver.computePostProcessing(spaces, &sol, postProcData));
    for (auto val : postProcData.res)
    {
      postProcFile << val << "\t";
    }
    postProcFile << "\n";
    postProcFile.flush();

    // // Write the isoline 0 if it is a CHNS problem
    // if(solver->isCHNS()) {
    //   std::ofstream isolineFile;
    //   std::string name = solver->getWriteDir() + "isoline_iConv" +
    //   std::to_string(iConv) + "_ifp" + std::to_string(ifp) + "_endOfInterval"
    //   + std::to_string(iI) + ".txt"; isolineFile.open(name);

    //   // Print time
    //   isolineFile << sol.getCurrentTime() << "\n";

    //   // Print points on the isoline
    //   for(size_t iPt = 0; iPt < postProcData.isoline.size() / 3; ++iPt) {
    //     isolineFile << postProcData.isoline[3 * iPt + 0] << "\t" <<
    //     postProcData.isoline[3 * iPt + 1] << "\n";
    //   }
    //   isolineFile.close();
    // }

    delete exporter;

    for (int jI = 0; jI < nI; ++jI)
    {
      if (allContainers[jI] != nullptr)
        feInfo("Computed interval %d - Container %d has %d DOF",
               iI,
               jI,
               allContainers[jI]->getNbDOFs());
    }
  } // for iI (sub-intervals)

  if (ap.adapt)
  {
    this->computeSpaceTimeMetrics();
    // this->adapt();
  }

  ++iFixedPoint;

  return FE_STATUS_OK;
}

feStatus TransientAdapter::computeSpaceTimeErrors(
  std::vector<std::vector<double>> &spaceTimeErrors)
{
  const Parameters::MeshAdaptation  &ap        = _adapt_parameters;
  const Parameters::TimeIntegration &tp        = _time_parameters;
  const int                          nI        = ap.nIntervals;
  const size_t                       numFields = allErrors[0][0].size();

  spaceTimeErrors.clear();
  spaceTimeErrors.resize(_adapt_parameters.nFixedPoint);

  for (int ifp = 0; ifp < _adapt_parameters.nFixedPoint; ++ifp)
  {
    // The space-time (integral) error for each field and this fixed-point
    // iteration
    spaceTimeErrors[ifp].resize(numFields, 0.);

    for (int iI = 0; iI < nI; ++iI)
    {
      for (size_t j = 0; j < numFields; ++j)
      {
        switch (ap.spaceTimeNorm)
        {
          case Parameters::MeshAdaptation::SpaceTimeErrorNorm::L1:
            // L1 error from Alauzet et al.'s paper
            spaceTimeErrors[ifp][j] +=
              tp.dt * ap.nTimeStepsPerInterval * allErrors[ifp][iI][j];
            break;
          case Parameters::MeshAdaptation::SpaceTimeErrorNorm::Linf:
            // Linf error in time
            spaceTimeErrors[ifp][j] =
              fmax(spaceTimeErrors[ifp][j], allErrors[ifp][iI][j]);
            break;
        }
      }
    }

    // // Print error of the last fixed-point iteration to file
    // if (ifp == adapter.nFixedPoint - 1)
    // {
    //   const int targetNst =
    //     nI * adapter.nVertices * adapter.nTimeStepsPerInterval;

    //   errorFile << adapter.nVertices << "\t";
    //   errorFile << (effectiveNst / adapter.nTimeStepsPerInterval)
    //             << "\t"; // sum of all effective vertices
    //   errorFile << nI << "\t";
    //   errorFile << adapter.nTimeStepsPerInterval << "\t";

    //   errorFile << targetNst << "\t";
    //   errorFile << effectiveNst << "\t";

    //   errorFile << adapter.dt;
    //   for (int j = 0; j < solver->getNumFields(); ++j)
    //   {
    //     errorFile << "\t" << std::setprecision(8) << e_st[j];
    //   }
    //   errorFile << "\n";
    //   errorFile.flush();
    // }
  }

  return FE_STATUS_OK;
}

feStatus TransientAdapter::finalize()
{
  postProcFile.close();

  return FE_STATUS_OK;
}