
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

  referenceTestCases.resize(ap.nIntervalsReferenceSolution);

  // Metric options common to all meshes
  feMetricOptions metricOptions(io.initialMesh);
  metricOptions.polynomialDegree = ap.orderForAdaptation;

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

  allOptions.resize(nI, metricOptions);

  for (int iI = 0; iI < nI; ++iI)
  {
    // Different adapted mesh name for each sub-interval
    allOptions[iI].adaptedMeshName =
      io.writeDir + "adapted_ifp_0_interval_" + std::to_string(iI) + ".msh";

    allMeshes[iI] = new feMesh2DP1(io.initialMesh);

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

feStatus TransientAdapter::updateBeforeFixedPointIteration()
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

  return FE_STATUS_OK;
}

feStatus TransientAdapter::computeSpaceTimeMetrics()
{
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
    metricField->applyGradation();
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
    feInfo("%d : %d", iI, numVerticesInAdaptedMeshes[iI]);
  int totalVertices = std::accumulate(numVerticesInAdaptedMeshes.begin(),
                                      numVerticesInAdaptedMeshes.end(),
                                      0);

  const int effectiveNst = totalVertices * ap.nTimeStepsPerInterval;

  feInfo("Total  number of space      vertices : %d", totalVertices);
  feInfo("Total  number of space-time vertices : %d", effectiveNst);
  feInfo("Target number of space-time vertices : %f", targetNst);

  return FE_STATUS_OK;
}

feStatus TransientAdapter::adaptAllMeshes()
{
  // #if defined(HAVE_GMSH)
  //   const Parameters::MeshAdaptation &ap = _adapt_parameters;
  //   const int    nI     = ap.nIntervals;

  //   std::vector<int> numVerticesInAdaptedMeshes(nI, 0);

  //   //
  //   // Adapt on each sub-interval
  //   //
  //   for(int iI = 0; iI < nI; ++iI)
  //   {
  //     feMetric *metricField = allMetrics[iI];
  //     // metricField->scaleMetricsByDeterminant(globalScalingFactor,
  //     expDeterminantLocal);

  //     initializeGmsh();
  //     feMetricOptions my_options = metricField->getOptions();
  //     gmsh::open(my_options.backgroundMeshfile);
  //     // openGmshModel(my_options.backgroundMeshfile);
  //     metricField->setDebugPrefix(std::to_string(iI));
  //     // metricField->applyGradation();
  //     metricField->writeMetricField();
  //     setVerbose(0);
  //     allMeshes[iI]->adapt(metricField);
  //     setVerbose(1);

  //     // Read the adapted mesh to check the number of vertices
  //     feMesh2DP1 mesh(my_options.adaptedMeshName);
  //     numVerticesInAdaptedMeshes[iI] = mesh.getNumVertices();
  //   }


  //   feInfo("Number of vertices in adapted meshes :");
  //   for(int iI = 0; iI < nI; ++iI)
  //     feInfo("%d : %d", iI, numVerticesInAdaptedMeshes[iI]);
  //   int totalVertices = std::accumulate(numVerticesInAdaptedMeshes.begin(),
  //   numVerticesInAdaptedMeshes.end(), 0);

  //   const double targetNst = (double) ap.targetVertices * nI;
  //   const int effectiveNst = totalVertices * ap.nTimeStepsPerInterval;

  //   feInfo("Total  number of space      vertices : %d", totalVertices);
  //   feInfo("Total  number of space-time vertices : %d", effectiveNst);
  //   feInfo("Target number of space-time vertices : %f", targetNst);
  // #else
  //   feWarning("Gmsh is required to write the metric field");
  // #endif
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
    if (exportSolution)
    {
      feCheck(createVisualizationExporter(
        exporter, VTK, numbering, &sol, mesh, spaces));
    }
    feExportData exportData = {exporter, exportFrequency, vtkFileRoot};

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

  if (adapt)
  {
    this->computeSpaceTimeMetrics();
    // this->adapt();
  }

  ++iFixedPoint;

  return FE_STATUS_OK;
}

feStatus TransientAdapter::finalize()
{
  postProcFile.close();

  return FE_STATUS_OK;
}