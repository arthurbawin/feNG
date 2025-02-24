
#include "feTimeIntegration.h"
#include "feMessage.h"

extern int FE_VERBOSE;

feStatus createTimeIntegrator(TimeIntegrator *&solver, timeIntegratorScheme scheme,
                              feNLSolverOptions NLoptions, feLinearSystem *system,
                              feSolution *solution, feMesh *mesh,
                              std::vector<feNorm*> &norms, feExportData exportData, double tBegin,
                              double tEnd, int nTimeSteps)
{
  UNUSED(mesh, tBegin, tEnd, nTimeSteps);
  switch(scheme) {
    case timeIntegratorScheme::STATIONARY:
      solver = new StationaryIntegrator(system, solution, mesh, NLoptions, norms, exportData);
      break;
    case timeIntegratorScheme::BDF1:
      solver = new BDF1Integrator(system, solution, mesh, NLoptions, norms, exportData, tBegin, tEnd, nTimeSteps);
      break;
    case timeIntegratorScheme::BDF2:
      solver = new BDF2Integrator(system, solution, mesh, NLoptions, norms, exportData, tBegin, tEnd, nTimeSteps);
      break;
    // case DC2F:
    //   solver = new DC2FSolver(NLoptions, metaNumber, system, solution, norms, mesh, exportData,
    //                           tBegin, tEnd, nTimeSteps);
    //   break;
    // case DC3:
    //   solver = new DC3Solver(NLoptions, metaNumber, system, solution, norms, mesh, exportData,
    //                          tBegin, tEnd, nTimeSteps, initializationCode);
    //   break;
    // case DC3F:
    //   solver = new DC3FSolver(NLoptions, metaNumber, system, solution, norms, mesh, exportData,
    //                           tBegin, tEnd, nTimeSteps);
    //   break;
    default:
      return feErrorMsg(FE_STATUS_ERROR, "Unsupported time integration scheme.");
  }
  return FE_STATUS_OK;
}

TimeIntegrator::TimeIntegrator(feLinearSystem *linearSystem,
                 feSolution *sol,
                 feMesh *mesh,
                 feNLSolverOptions &NLoptions,
                 std::vector<feNorm*> &postProcessing,
                 feExportData &exportData,
                 int order,
                 double t0, double tEnd, int nTimeSteps)
  : _linearSystem(linearSystem), _currentSolution(sol), _mesh(mesh), _NLoptions(NLoptions), _postProcessing(postProcessing),
  _exportData(exportData), _order(order), _t0(t0), _tEnd(tEnd), _nTimeSteps(nTimeSteps), _currentTime(t0), _currentStep(0)
{
  _dt = (tEnd - t0) / (double) nTimeSteps;
  _currentSolution->initializeTemporalSolution(_t0, _tEnd, _nTimeSteps);
  _timeHistory.resize(_nTimeSteps, 0.);

  // Allocate vectors for post-processing operators
  _postProcessingData.resize(postProcessing.size());
  for(auto &vec : _postProcessingData)
    vec.resize(nTimeSteps, 0.);
}

//
// Initialize the time and solution array
//
void TimeIntegrator::initialize()
{
  _t[0] = _currentTime;
  _currentSolution->setCurrentTime(_currentTime);
  _currentSolution->initializeUnknowns(_mesh);
  _currentSolution->initializeEssentialBC(_mesh);
}

void TimeIntegrator::updateTime(const int iStep, const double dt)
{
  // Rotate time and time difference of previous steps
  for(int i = _nSol - 1; i > 0; i--) {
    _t[i] = _t[i-1];
    _deltaT[i] = _deltaT[i-1];
  }

  _currentTime += dt;
  _t[0] =  _currentTime;
  _deltaT[0] = _t[0] - _t[1];
  _timeHistory[iStep] = _currentTime;

  _currentSolution->setCurrentTime(_currentTime);
}

feStatus TimeIntegrator::writeVisualizationFile(const int iStep, const feExportData &data)
{
  if(data.exporter == nullptr)
    return FE_STATUS_OK;

  if((iStep % data.exportEveryNSteps) == 0) {
    std::string fileName = data.fileNameRoot + "_" + _integratorName + "_" + std::to_string(iStep) + ".vtk";
    if(_integratorName == "Stationary")
      fileName = data.fileNameRoot + "_" + _integratorName + ".vtk";
    feStatus s = data.exporter->writeStep(fileName);
    return s;
  }

  return FE_STATUS_OK;
}

//
// Stationary solver
//
StationaryIntegrator::StationaryIntegrator(feLinearSystem *linearSystem,
          					                       feSolution *sol,
                                           feMesh *mesh,
          					                       feNLSolverOptions &NLoptions,
          					                       std::vector<feNorm*> &postProcessing,
          					                       feExportData &exportData)
  : TimeIntegrator(linearSystem, sol, mesh, NLoptions, postProcessing, exportData, 0)
{
  _integratorName = "Stationary";
  _nSol = 1;
  _t.resize(_nSol, 0.);
  _deltaT.resize(_nSol, 0.);
  _sC = new BDFContainer(_nSol, _currentTime, _currentSolution->getNumDOFs());
  this->initialize();
  // Copy the initialized solution into the container
  _sC->setCurrentSolution(_currentSolution);
  _sC->setCurrentSolutionDot(_currentSolution);
}

feStatus StationaryIntegrator::restartFromContainer(const feSolutionContainer &/* container */)
{
  // Nothing to do for a stationary solver
  return FE_STATUS_OK;
}

feStatus StationaryIntegrator::makeStep(double /* dt */)
{
  return this->makeSteps(0);
}

feStatus StationaryIntegrator::makeSteps(int /* nSteps */)
{
  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "STATIONARY SOLVER:");
  feInfoCond(FE_VERBOSE > 0, "\t\tSolving for steady state solution");

  _sC->computeBDFCoefficients(_order, _deltaT);

  // Solve nonlinear problem
  feStatus s = solveNewtonRaphson(_linearSystem, _currentSolution, _sC, _NLoptions);
  if(s == FE_STATUS_ERROR) { return s; }

  feInfoCond(FE_VERBOSE > 0, "\t\tSteady-state solution found");

  // Compute post-processing indicators (L2 norms, etc.)
  for(size_t i = 0; i < _postProcessing.size(); ++i) {
    _postProcessingData[i][0] = _postProcessing[i]->compute();
  }

  // Export solution for visualization
  s = this->writeVisualizationFile(0, _exportData);
  if(s != FE_STATUS_OK) return s;

  return FE_STATUS_OK;
}

//
// BDF1 time integrator
//
BDF1Integrator::BDF1Integrator(feLinearSystem *linearSystem,
                               feSolution *sol,
                               feMesh *mesh,
                               feNLSolverOptions &NLoptions,
                               std::vector<feNorm*> &postProcessing,
                               feExportData &exportData,
                               double t0, double tEnd, int nTimeSteps)
  : TimeIntegrator(linearSystem, sol, mesh, NLoptions,
    postProcessing, exportData, 1, t0, tEnd, nTimeSteps)
{
  _integratorName = "BDF1";
  // Initialize the solution container with 2 solutions
  _nSol = 2;
  _t.resize(_nSol, 0.);
  _deltaT.resize(_nSol, 0.);
  _sC = new BDFContainer(_nSol, _currentTime, _currentSolution->getNumDOFs());
  this->initialize();
  // Copy the initialized solution into the container
  _sC->setCurrentSolution(_currentSolution);
  _sC->setCurrentSolutionDot(_currentSolution);
}

feStatus BDF1Integrator::restartFromContainer(const feSolutionContainer &container)
{
  if(_nSol != container.getNbSol()) {
    return feErrorMsg(FE_STATUS_ERROR, "Could not restart %s time integrator from solution container: "
      "integrator needs %d solution arrays but container stores %d",
      _integratorName.data(),
      _nSol,
      container.getNbSol());
  }

  int nDOF = _currentSolution->getNumDOFs();
  if(nDOF != container.getNbDOFs()) {
    return feErrorMsg(FE_STATUS_ERROR, "Could not restart %s time integrator from solution container: "
      "current solution has %d DOF but container stores %d DOF per solution",
      _integratorName.data(),
      nDOF,
      container.getNbDOFs());
  }

  // Can't use copy assignment because a BDFContainer is specialized,
  // but only through the computeTimeDerivative function...
  // This should be reworked entirely.
  // In the meantime:
  // Copy data from container : 
  const std::vector<double> &tC = container.getTime();
  const std::vector<double> &deltaTC = container.getTimeDifferences();

  for(int i = 0; i < _nSol; ++i) {
    _t[i] = tC[i];
    _deltaT[i] = deltaTC[i];
    _sC->setSolution(container.getSolution(i), i);
    _sC->setSolutionDot(container.getSolutionDot(i), i);
    _sC->setResidual(container.getFResidual(i), i);
  }

  return FE_STATUS_OK;
}

feStatus BDF1Integrator::makeStep(double /* dt */)
{
  return feErrorMsg(FE_STATUS_ERROR, "Reimplement single makeStep");
}

feStatus BDF1Integrator::makeSteps(int nSteps)
{
  if(_quietStart) {
    feInfoCond(FE_VERBOSE > 0, "\t\t\tBDF1 : Advancing from t = %1.4e to t = %1.4e in %d step(s)", 
      _currentTime, _currentTime + nSteps * _dt, nSteps);
  } else {
    feInfoCond(FE_VERBOSE > 0, "");
    feInfoCond(FE_VERBOSE > 0, "BDF1 SOLVER:");
    feInfoCond(FE_VERBOSE > 0, "\t\tAdvancing from t = %1.4e to t = %1.4e in %d step(s)", 
      _currentTime, _currentTime + nSteps * _dt, nSteps);
  }

  if(_currentStep == 0)
  {
    if(_exportInitialCondition) {
      // Export initial solution
      feStatus s = this->writeVisualizationFile(0, _exportData);
      if(s != FE_STATUS_OK) return s;
    }
  }

  for(int iStep = 0; iStep < nSteps; ++iStep, ++_currentStep)
  {
    if(!_quietStart) {
      feInfoCond(FE_VERBOSE > 0, "\t\tBDF1 - Overall step = %d - Step %d/%d from t = %1.4e to t = %1.4e", 
        _currentStep+1, iStep+1, nSteps, _currentTime, _currentTime + _dt);
    }

    this->updateTime(iStep, _dt);
    _sC->rotate(_dt);
    _sC->computeBDFCoefficients(_order, _deltaT);
    _currentSolution->setC0(_sC->getC0());
    _currentSolution->initializeEssentialBC(_mesh, _sC);

    // Solve nonlinear problem
    feStatus s = solveNewtonRaphson(_linearSystem, _currentSolution, _sC, _NLoptions);
    if(s == FE_STATUS_ERROR) { return s; }

    // Compute post-processing indicators (L2 norms, etc.)
    for(size_t i = 0; i < _postProcessing.size(); ++i) {
      _postProcessingData[i][0] = _postProcessing[i]->compute();
    }

    // Export solution for visualization
    s = this->writeVisualizationFile(_currentStep+1, _exportData);
    if(s != FE_STATUS_OK) return s;
  }

  return FE_STATUS_OK;
}

//
// BDF2 time integrator
//
BDF2Integrator::BDF2Integrator(feLinearSystem *linearSystem,
                               feSolution *sol,
                               feMesh *mesh,
                               feNLSolverOptions &NLoptions,
                               std::vector<feNorm*> &postProcessing,
                               feExportData &exportData,
                               double t0, double tEnd, int nTimeSteps)
  : TimeIntegrator(linearSystem, sol, mesh, NLoptions,
    postProcessing, exportData, 2, t0, tEnd, nTimeSteps)
{
  _integratorName = "BDF2";
  // Initialize the solution container with 2 solutions
  _nSol = 3;
  _t.resize(_nSol, 0.);
  _deltaT.resize(_nSol, 0.);
  _sC = new BDFContainer(_nSol, _currentTime, _currentSolution->getNumDOFs());
  this->initialize();
  // Copy the initialized solution into the container
  _sC->setCurrentSolution(_currentSolution);
  _sC->setCurrentSolutionDot(_currentSolution);
}

feStatus BDF2Integrator::restartFromContainer(const feSolutionContainer &container)
{
  if(_nSol != container.getNbSol()) {
    return feErrorMsg(FE_STATUS_ERROR, "Could not restart %s time integrator from solution container: "
      "integrator needs %d solution arrays but container stores %d",
      _integratorName.data(),
      _nSol,
      container.getNbSol());
  }

  int nDOF = _currentSolution->getNumDOFs();
  if(nDOF != container.getNbDOFs()) {
    return feErrorMsg(FE_STATUS_ERROR, "Could not restart %s time integrator from solution container: "
      "current solution has %d DOF but container stores %d DOF per solution",
      _integratorName.data(),
      nDOF,
      container.getNbDOFs());
  }

  // Copy data from container :
  const std::vector<double> &tC = container.getTime();
  const std::vector<double> &deltaTC = container.getTimeDifferences();

  for(int i = 0; i < _nSol; ++i) {
    _t[i] = tC[i];
    _deltaT[i] = deltaTC[i];
    _sC->setSolution(container.getSolution(i), i);
    _sC->setSolutionDot(container.getSolutionDot(i), i);
    _sC->setResidual(container.getFResidual(i), i);
  }

  return FE_STATUS_OK;
}

feStatus BDF2Integrator::makeStep(double /* dt */)
{
  return feErrorMsg(FE_STATUS_ERROR, "Reimplement single makeStep");
}

feStatus BDF2Integrator::startWithBDF1(const double timeStepRatio)
{
  double dtBDF1 = timeStepRatio * _dt;
  double dt2 = _dt - dtBDF1;

  // Assign initial condition to solution at position n-1
  _sC->setSolution(_currentSolution->getSolution(), 2);
  _sC->setSolutionDot(_currentSolution->getSolutionDot(), 2);

  BDF1Integrator solver(_linearSystem, _currentSolution, 
    _mesh, _NLoptions, _postProcessing, _exportData, 0, dtBDF1, 1);
  solver.setQuietStart(true);
  solver.setExportInitialCondition(false);
  feStatus s = solver.makeSteps(1);
  if(s == FE_STATUS_ERROR) { return s; }

  feInfoCond(FE_VERBOSE > 0, "\t\t\tContinuing with BDF2 from t = %1.4e to t = %1.4e", 
    _currentTime, _currentTime + _dt);

  // Assign BDF1 solution at position n and n+1 ("current solution")
  _sC->setSolution(_currentSolution->getSolution(), 1);
  _sC->setSolutionDot(_currentSolution->getSolutionDot(), 1);
  _sC->setSolution(_currentSolution->getSolution(), 0);
  _sC->setSolutionDot(_currentSolution->getSolutionDot(), 0);

  // Manually update time in BDF2 solver and container
  // Time and deltaT at position 2 are already 0
  _currentTime = _dt;
  _t[0] = _currentTime;
  _deltaT[0] = _currentTime - dtBDF1;
  _sC->setTime(0, _dt);
  _sC->setTimeDifference(1, _currentTime - dtBDF1);
  _t[1] = dtBDF1;
  _deltaT[1] = dtBDF1;
  _sC->setTime(1, dtBDF1);
  _sC->setTimeDifference(1, dtBDF1);

  std::vector<double> bdfCoeff(3);
  bdfCoeff[0] = (dtBDF1 + 2.*dt2)/(dtBDF1 + dt2) * 1. / dt2;
  bdfCoeff[1] = - (1. / dt2 + 1. / dtBDF1);
  bdfCoeff[2] = dt2/(dtBDF1 + dt2) * 1. / dtBDF1;
  _sC->setBDFCoefficients(bdfCoeff);

  _currentSolution->setC0(_sC->getC0());
  _currentSolution->setCurrentTime(_currentTime);
  _currentSolution->initializeEssentialBC(_mesh, _sC);

  // Solve nonlinear problem
  s = solveNewtonRaphson(_linearSystem, _currentSolution, _sC, _NLoptions);
  if(s == FE_STATUS_ERROR) { return s; }

  return FE_STATUS_OK;
}

//
// Integrate nSteps.
//
// If the integrator starts from t0 = 0, we need to start the BDF2
// then perform n-1 steps.
//
// If the integrator restarts from an arbitrary t0, the solution container
// is supposed to contain the solutions for the last 2 time steps, and we
// can simply integrate n steps.
//
feStatus BDF2Integrator::makeSteps(int nSteps)
{
  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "BDF2 SOLVER:");
  feInfoCond(FE_VERBOSE > 0, "\t\tAdvancing from t = %1.4e to t = %1.4e in %d step(s)", 
    _currentTime, _currentTime + nSteps * _dt, nSteps);

  int start = 0;

  if(_currentStep == 0)
  {
    start = 1;

    // Export initial solution
    feStatus s = this->writeVisualizationFile(0, _exportData);
    if(s != FE_STATUS_OK) return s;

    switch(_startMethod) {
      case BDF2Starter::BDF1:
      {
        // Use BDF1 to get intermediary solution at time t0 + timeStepRatio * dt
        // Then use variable time step BDF2 to get first solution at t0 + dt
        double timeStepRatio = 0.1;
        feInfoCond(FE_VERBOSE > 0, "\t\tBDF2 - Overall step = 1 - Step 1/%d from t = %1.4e to t = %1.4e", 
          nSteps, _currentTime, _currentTime + _dt);
        feInfoCond(FE_VERBOSE > 0, "\t\t\tStarting with BDF1 and intermediary time step %1.4e", timeStepRatio * _dt);

        s = startWithBDF1(timeStepRatio);
        if(s != FE_STATUS_OK) return s;
        break;
      }

      case BDF2Starter::InitialCondition:
        // Simply set the prescribed solution at the first time step
        // Start the BDF2 by computing initial conditions at time n and n-1
        // Only for verification purposes with exact solution, and
        // not suitable for actual simulation where solution at first step is not known.
        this->updateTime(0, _dt);
        _sC->rotate(_dt);
        _currentSolution->initializeUnknowns(_mesh);
        _currentSolution->initializeEssentialBC(_mesh);
        _sC->setCurrentSolution(_currentSolution);
        break;
    }

    // Export solution for visualization
    _currentStep++;
    s = this->writeVisualizationFile(_currentStep, _exportData);
    if(s != FE_STATUS_OK) return s;
  }

  for(int iStep = start; iStep < nSteps; ++iStep, ++_currentStep)
  {
    feInfoCond(FE_VERBOSE > 0, "\t\tBDF2 - Overall step = %d - Step %d/%d from t = %1.4e to t = %1.4e", 
      _currentStep+1, iStep+1, nSteps, _currentTime, _currentTime + _dt);

    this->updateTime(iStep, _dt);
    _sC->rotate(_dt);
    _sC->computeBDFCoefficients(_order, _deltaT);
    _currentSolution->setC0(_sC->getC0());
    // _currentSolution->setCurrentTime(_currentTime); // Done in updateTime
    _currentSolution->initializeEssentialBC(_mesh, _sC);

    // Solve nonlinear problem
    feStatus s = solveNewtonRaphson(_linearSystem, _currentSolution, _sC, _NLoptions);
    if(s == FE_STATUS_ERROR) { return s; }

    // Compute post-processing indicators (L2 norms, etc.)
    for(size_t i = 0; i < _postProcessing.size(); ++i) {
      _postProcessingData[i][0] = _postProcessing[i]->compute();
    }

    // Export solution for visualization
    s = this->writeVisualizationFile(_currentStep+1, _exportData);
    if(s != FE_STATUS_OK) return s;
  }

  return FE_STATUS_OK;
}