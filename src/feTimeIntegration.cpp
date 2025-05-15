
#include "feTimeIntegration.h"
#include "feMessage.h"

extern int FE_VERBOSE;

feStatus createTimeIntegrator(TimeIntegrator *&solver, timeIntegratorScheme scheme,
                              feNLSolverOptions NLoptions, feLinearSystem *system,
                              feSolution *solution, feMesh *mesh,
                              std::vector<feNorm*> &norms, feExportData exportData, double tBegin,
                              double tEnd, int nTimeSteps)
{
  switch(scheme) {
    case timeIntegratorScheme::STATIONARY:
      solver = new StationaryIntegrator(system, solution, mesh, NLoptions, norms, exportData);
      break;

    // BDF schemes
    case timeIntegratorScheme::BDF1:
      solver = new BDF1Integrator(system, solution, mesh, NLoptions, norms, exportData, tBegin, tEnd, nTimeSteps);
      break;
    case timeIntegratorScheme::BDF2:
      solver = new BDF2Integrator(system, solution, mesh, NLoptions, norms, exportData, tBegin, tEnd, nTimeSteps);
      break;

    // DC/BDF1 schemes
    case timeIntegratorScheme::DC2BDF1:
      solver = new DCBDF1Integrator(2, system, solution, mesh, NLoptions, norms, exportData, tBegin, tEnd, nTimeSteps);
      break;
    case timeIntegratorScheme::DC3BDF1:
      solver = new DCBDF1Integrator(3, system, solution, mesh, NLoptions, norms, exportData, tBegin, tEnd, nTimeSteps);
      break;

    // DC/BDF2 schemes
    case timeIntegratorScheme::DC3BDF2:
      solver = new DC3BDF2Integrator(system, solution, mesh, NLoptions, norms, exportData, tBegin, tEnd, nTimeSteps);
      break;

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

  // Allocate vectors for post-processing operators.
  // Norms are computed for the initial condition as well,
  // thus nTimeSteps + 1 times.
  _postProcessingData.resize(postProcessing.size() + 1);
  for(auto &vec : _postProcessingData)
    vec.resize(nTimeSteps+1, 0.);
}

//
// Initialize the time and solution array
//
void TimeIntegrator::initialize()
{
  _t[0] = _currentTime;
  _currentSolution->setCurrentTime(_currentTime);
  _currentSolution->initialize(_mesh);
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

void TimeIntegrator::resetTime()
{
  _currentTime = _t0;
  for(int i = _nSol - 1; i > 0; i--) {
    _t[i] = 0.;
    _deltaT[i] = 0.;
  }
  _t[0] = _currentTime;
}

void TimeIntegrator::computePostProcessing(const int iStep, const double t)
{
  // First vector is the time array
  _postProcessingData[0][iStep] = t;

  // Remaining are the results of the norms and indicators
  for(size_t i = 0; i < _postProcessing.size(); ++i) {
    _postProcessingData[i+1][iStep] = _postProcessing[i]->compute();
  }
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

feStatus TimeIntegrator::postProcessingAndVisualization(const int localStep, const int globalStep,
                                                        const double t, const feExportData &data)
{
  this->computePostProcessing(localStep, t);
  feCheckReturn(this->writeVisualizationFile(globalStep, data));
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
  computePostProcessing(0, 0.);

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

  *_sC = container;

  for(int i = 0; i < _nSol; ++i) {
    _t[i] = tC[i];
  }
  for(int i = 0; i < _nSol-1; ++i) {
    _deltaT[i] = _t[i] - _t[i+1];
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

    // Compute post-processing indicators (L2 norms, etc.)
    // on initial condition
    computePostProcessing(0, _t[0]);
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
    feCheckReturn(solveNewtonRaphson(_linearSystem, _currentSolution, _sC, _NLoptions));

    // Compute post-processing indicators (L2 norms, etc.)
    computePostProcessing(iStep+1, _t[0]);

    // Export solution for visualization
    feCheckReturn(this->writeVisualizationFile(_currentStep+1, _exportData));
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
  // Initialize the solution container with 3 solutions
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
  // const std::vector<double> &deltaTC = container.getTimeDifferences();

  // Call custom copy constructor for the internal solution container
  *_sC = container;

  for(int i = 0; i < _nSol; ++i) {
    _t[i] = tC[i];
  }
  for(int i = 0; i < _nSol-1; ++i) {
    _deltaT[i] = _t[i] - _t[i+1];
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

int currentBDFStep;

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

    // Compute post-processing indicators (L2 norms, etc.)
    // on initial condition
    computePostProcessing(0, _t[0]);

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
        feInfoCond(FE_VERBOSE > 0, "\t\tBDF2 - Overall step = 1 - Step 1/%d from t = %1.4e to t = %1.4e", 
          nSteps, _currentTime, _currentTime + _dt);
        feInfoCond(FE_VERBOSE > 0, "\t\t\tStarting with exact solution");
        this->updateTime(0, _dt);
        _sC->rotate(_dt);
        _currentSolution->initializeUnknowns(_mesh);
        _currentSolution->initializeEssentialBC(_mesh);
        _sC->setCurrentSolution(_currentSolution);
        break;
    }

    computePostProcessing(1, _t[0]);

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

    // feInfo("Integrator info:");
    // feInfo("_nSol            : %d", _nSol);
    // feInfo("_order           : %d", _order);
    // feInfo("_t0              : %+-1.20e", _t0);
    // feInfo("_tEnd            : %+-1.20e", _tEnd);
    // feInfo("_dt              : %+-1.20e", _dt);
    // feInfo("_nTimeSteps :    : %d", _nTimeSteps);
    // feInfo("_currentTime     : %+-1.20e", _currentTime);
    // feInfo("_currentStep     : %d", _currentStep);
    // feInfo("time             : %+-1.20e - %+-1.20e - %+-1.20e", _t[0], _t[1], _t[2]);
    // feInfo("deltaT           : %+-1.20e - %+-1.20e - %+-1.20e", _deltaT[0], _deltaT[1], _deltaT[2]);
    // feInfo("time history     : %+-1.20e - %+-1.20e - %+-1.20e", _timeHistory[0], _timeHistory[1], _timeHistory[2]);
    // feInfo("BDF coefficients : %+-1.20e - %+-1.20e - %+-1.20e",
    //   _sC->getBDFCoefficients()[0],
    //   _sC->getBDFCoefficients()[1],
    //   _sC->getBDFCoefficients()[2]);

    _currentSolution->setC0(_sC->getC0());
    // _currentSolution->setCurrentTime(_currentTime); // Done in updateTime

    // for(const auto &val : _currentSolution->getSolution()) {
    //   if(std::isnan(val)) {
    //     feInfo("Nan alert before initialize BC");
    //     exit(-1);
    //   }
    // }

    _currentSolution->initializeEssentialBC(_mesh, _sC);

    // Solve nonlinear problem
    currentBDFStep = _currentStep+1;
    feStatus s = solveNewtonRaphson(_linearSystem, _currentSolution, _sC, _NLoptions);
    if(s == FE_STATUS_ERROR) { return s; }

    // Compute post-processing indicators (L2 norms, etc.)
    computePostProcessing(iStep+1, _t[0]);

    // Export solution for visualization
    s = this->writeVisualizationFile(_currentStep+1, _exportData);
    if(s != FE_STATUS_OK) return s;
  }

  return FE_STATUS_OK;
}

//
// Deferred correction DC/BDF1 time integrators
//
DCBDF1Integrator::DCBDF1Integrator(const int order,
                                   feLinearSystem *linearSystem,
                                   feSolution *sol,
                                   feMesh *mesh,
                                   feNLSolverOptions &NLoptions,
                                   std::vector<feNorm*> &postProcessing,
                                   feExportData &exportData,
                                   double t0, double tEnd, int nTimeSteps)
  : TimeIntegrator(linearSystem, sol, mesh, NLoptions,
    postProcessing, exportData, order, t0, tEnd, nTimeSteps)
{
  _integratorName = "DC" + std::to_string(_order) + "BDF1";
  _nSol = 2;
  _t.resize(_nSol, 0.);
  _deltaT.resize(_nSol, 0.);

  // Initialize time and solution wrapper
  this->initialize();

  // Copy the initialized solution into the containers
  _scBDF1 = new BDFContainer(_nSol, _currentTime, _currentSolution->getNumDOFs());
  _scBDF1->setCurrentSolution(_currentSolution);
  // _scBDF1->setCurrentSolutionDot(_currentSolution);

  // Initialize DC containers of decreasing order.
  // First solution is the DC-_order solution
  // Second is the DC-_order-1 solution, etc.
  //
  // The highest order DC container stores 2 solutions, as for the BDF1 scheme,
  // but intermediary DC containers must store _order solutions to improve the previous one.
  int nSol_intermediaryDC = _order;

  _scDC.resize(_order-1);
  for(size_t i = 0; i < _scDC.size(); ++i)
  {
    if(i == 0)
      _scDC[i] = new DCBDFContainer(_order - i, 1, _nSol, _currentTime, _currentSolution->getNumDOFs());
    else
      _scDC[i] = new DCBDFContainer(_order - i, 1, nSol_intermediaryDC, _currentTime, _currentSolution->getNumDOFs());

    _scDC[i]->setCurrentSolution(_currentSolution);
  }
}

feStatus DCBDF1Integrator::restartFromContainer(const feSolutionContainer &container)
{
  UNUSED(container);
  return feErrorMsg(FE_STATUS_ERROR, "Reimplement restart for DCBDF1Integrator");
}

feStatus DCBDF1Integrator::makeStep(double /* dt */)
{
  return feErrorMsg(FE_STATUS_ERROR, "Reimplement single makeStep");
}

feStatus DCBDF1Integrator::makeSteps(int nSteps)
{
  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "DC%d/BDF1 SOLVER:", _order);
  feInfoCond(FE_VERBOSE > 0, "\t\tAdvancing from t = %1.4e to t = %1.4e in %d step(s)", 
    _currentTime, _currentTime + nSteps * _dt, nSteps);

  int start = 0;

  //
  // First time step
  // 
  if(_currentStep == 0)
  {
    start = 1;

    feInfoCond(FE_VERBOSE > 0, "\t\tDC%d/BDF1 - Overall step = 1 - Step 1/%d from t = %1.4e to t = %1.4e", 
      _order, nSteps, _currentTime, _currentTime + _dt);
    feInfoCond(FE_VERBOSE > 0, "");
    feInfoCond(FE_VERBOSE > 0, "\t\tBeginning of DC%d/BDF1 initialization procedure", _order);

    //
    // Compute time derivative at t = t_0, required for the first corrections.
    // Set this first uDot in each container.
    //
    // feInfo("Entering first solve for udot with containers:");
    // for(size_t i = 0; i < _scBDF1->_sol[0].size(); ++i) {
    //   feInfo("BDF1: s[0] = %+-1.2e s[1] = %+-1.2e [[d[0] = %+-1.2e]] d[1] = %+-1.2e  "
    //            "DC2: s[0] = %+-1.2e s[1] = %+-1.2e d[0] = %+-1.2e d[1] = %+-1.2e c = %+-1.3e ",
    //       _scBDF1->_sol[0][i], _scBDF1->_sol[1][i], _scBDF1->_solDot[0][i], _scBDF1->_solDot[1][i],
    //       _scDC[0]->_sol[0][i], _scDC[0]->_sol[1][i], _scDC[0]->_solDot[0][i], _scDC[0]->_solDot[1][i], _scDC[0]->_d[i]);
    //   if((int) i+1 == _linearSystem->getSystemSize())
    //     feInfo("=====================================================================================");
    // }

    // _scBDF1->setSolutionDot(_currentSolution->getSolutionDot(), 0);
    feInfoCond(FE_VERBOSE > 0, "\t\t\tSolving linear problem for initial time derivative");
    bool solveForTimeDerivative = true;
    feCheckReturn(solveNewtonRaphson(_linearSystem, _currentSolution, _scBDF1, _NLoptions, solveForTimeDerivative));
    // for(auto &val : _currentSolution->getSolutionDot()) {
    //   val = 0.;
    //   const double t = _t[0];
    //   val = -2./((1.+t)*(1.+t));
    // }

    // for(auto &val : _currentSolution->getSolutionDot())
    //   feInfo("computed soldot %f", val);
    //     // exit(-1);
    // for(auto &val : _scBDF1->_solDot[0])
    //   feInfo("bdf1     soldot %f", val);

    // uDot(t = t0) must be known by all but the most precise container
    // Here we just set it in all containers
    _scBDF1->setSolutionDot(_currentSolution->getSolutionDot(), 0);
    for(auto &container : _scDC) {
      container->setSolutionDot(_currentSolution->getSolutionDot(), 0);
    }

    // feInfo("Leaving first solve for udot with containers:");
    // for(size_t i = 0; i < _scBDF1->_sol[0].size(); ++i) {
    //   feInfo("BDF1: s[0] = %+-1.2e s[1] = %+-1.2e [[d[0] = %+-1.2e]] d[1] = %+-1.2e  "
    //            "DC2: s[0] = %+-1.2e s[1] = %+-1.2e d[0] = %+-1.2e d[1] = %+-1.2e c = %+-1.3e ",
    //       _scBDF1->_sol[0][i], _scBDF1->_sol[1][i], _scBDF1->_solDot[0][i], _scBDF1->_solDot[1][i],
    //       _scDC[0]->_sol[0][i], _scDC[0]->_sol[1][i], _scDC[0]->_solDot[0][i], _scDC[0]->_solDot[1][i], _scDC[0]->_d[i]);
    //   if((int) i+1 == _linearSystem->getSystemSize())
    //     feInfo("=====================================================================================");
    // }

    // Compute error after the initial uDot has been set
    feCheckReturn(postProcessingAndVisualization(0, 0, _t[0], _exportData));

    if(_order == 2) {
      // Nothing to do to start the DC2/BDF1
      start = 0;
    } else {
      // feInfo("Testing DC2/BDF1 first");
      // exit(-1);
      if(_order == 3) // DC3/BDF1
      {
        // feInfo("Testing DC2/BDF1 first");
        // exit(-1);
        // Advance BDF1 and DC2/BDF1 until t = t_2
        for(int iStep = 0; iStep < 2; ++iStep)
        {
          this->updateTime(iStep, _dt);

          // Advance BDF1
          feInfoCond(FE_VERBOSE > 0, "\t\t\tAdvancing BDF1 solution to t = %1.4e", _currentTime);
          _scBDF1->rotate(_dt);
          _scBDF1->computeBDFCoefficients(1, _deltaT);
          _currentSolution->setC0(_scBDF1->getC0());
          _currentSolution->setSolFromContainer(_scBDF1);
          _currentSolution->initializeEssentialBC(_mesh, _scBDF1);
          feCheckReturn(solveNewtonRaphson(_linearSystem, _currentSolution, _scBDF1, _NLoptions));

          // Advance DC2/BDF1
          feInfoCond(FE_VERBOSE > 0, "\t\t\tAdvancing DC2/BDF1 solution to t = %1.4e", _currentTime);
          DCBDFContainer *solDC = _scDC[1];
          solDC->rotate(_dt);
          solDC->computeBDFCoefficients(1, _deltaT);
          solDC->computeDCCoefficients(_deltaT);
          _currentSolution->setSolFromContainer(solDC);
          _currentSolution->initializeEssentialBC(_mesh, solDC);
          solDC->computeDCCorrection(*_scBDF1, _deltaT);
          feCheckReturn(solveNewtonRaphson(_linearSystem, _currentSolution, solDC, _NLoptions));
        }

        // Compute first step for DC3/BDF1
        this->resetTime();
        this->updateTime(0, _dt);
        DCBDFContainer *solDC = _scDC[0];
        solDC->rotate(_dt);
        feInfoCond(FE_VERBOSE > 0, "\t\t\tAdvancing DC3/BDF1 solution to t = %1.4e", _currentTime);
        solDC->computeFirstDCCoefficients(_deltaT);
        solDC->computeDCCorrection(*_scDC[1], _deltaT);
        solDC->computeBDFCoefficients(1, _deltaT);
        _currentSolution->setSolFromContainer(solDC);
        _currentSolution->initializeEssentialBC(_mesh, solDC);
        feCheckReturn(solveNewtonRaphson(_linearSystem, _currentSolution, solDC, _NLoptions));

        _scBDF1->invRotate(_dt);
        _scDC[1]->invRotate(_dt);
      }

      feCheckReturn(postProcessingAndVisualization(1, 1, _t[0], _exportData));
      _currentStep++;
    }

    feInfoCond(FE_VERBOSE > 0, "\t\tEnd of DC%d/BDF1 initialization procedure", _order);
    feInfoCond(FE_VERBOSE > 0, "");
  }

  //
  // Regular time steps
  //
  for(int iStep = start; iStep < nSteps; ++iStep, ++_currentStep)
  {
    feInfoCond(FE_VERBOSE > 0, "\t\tDC%d/BDF1 - Overall step = %d - Step %d/%d from t = %1.4e to t = %1.4e", 
      _order, _currentStep+1, iStep+1, nSteps, _currentTime, _currentTime + _dt);

    this->updateTime(iStep, _dt);

    //
    // Advance the BDF1 solution
    //
    feInfoCond(FE_VERBOSE > 0, "\t\t\tAdvancing BDF1 solution to t = %1.4e", _currentTime);
    _scBDF1->rotate(_dt);
    _scBDF1->computeBDFCoefficients(1, _deltaT);
    _currentSolution->setC0(_scBDF1->getC0());
    // feInfo("c0 = %f", _scBDF1->getC0());
    // feInfo("c0 = %f", _scBDF1->getBDFCoefficients()[0]);
    // feInfo("c1 = %f", _scBDF1->getBDFCoefficients()[1]);
    _currentSolution->setSolFromContainer(_scBDF1);
    _currentSolution->initializeEssentialBC(_mesh, _scBDF1);

    // Set nan in udot[0], which will be updated

    // feInfo("Entering BDF1 solve with containers:");
    // for(size_t i = 0; i < _scBDF1->_sol[0].size(); ++i) {
    //   // _scBDF1->_solDot[0][i] = nan("");
    //   feInfo("BDF1: s[0] = %+-1.2e s[1] = %+-1.2e [[d[0] = %+-1.2e]] d[1] = %+-1.2e  "
    //            "DC2: s[0] = %+-1.2e s[1] = %+-1.2e d[0] = %+-1.2e d[1] = %+-1.2e c = %+-1.3e ",
    //       _scBDF1->_sol[0][i], _scBDF1->_sol[1][i], _scBDF1->_solDot[0][i], _scBDF1->_solDot[1][i],
    //       _scDC[0]->_sol[0][i], _scDC[0]->_sol[1][i], _scDC[0]->_solDot[0][i], _scDC[0]->_solDot[1][i], _scDC[0]->_d[i]);
    //   if((int) i+1 == _linearSystem->getSystemSize())
    //     feInfo("=====================================================================================");
    // }

    feCheckReturn(solveNewtonRaphson(_linearSystem, _currentSolution, _scBDF1, _NLoptions));

    // feInfo("Leaving BDF1 solve with containers:");
    // for(size_t i = 0; i < _scBDF1->_sol[0].size(); ++i) {
    //   feInfo("BDF1: s[0] = %+-1.2e s[1] = %+-1.2e [[d[0] = %+-1.2e]] d[1] = %+-1.2e  "
    //            "DC2: s[0] = %+-1.2e s[1] = %+-1.2e d[0] = %+-1.2e d[1] = %+-1.2e c = %+-1.3e ",
    //       _scBDF1->_sol[0][i], _scBDF1->_sol[1][i], _scBDF1->_solDot[0][i], _scBDF1->_solDot[1][i],
    //       _scDC[0]->_sol[0][i], _scDC[0]->_sol[1][i], _scDC[0]->_solDot[0][i], _scDC[0]->_solDot[1][i], _scDC[0]->_d[i]);
    //   if((int) i+1 == _linearSystem->getSystemSize())
    //     feInfo("=====================================================================================");
    // }

    //
    // Advance the DC solutions
    // At each step, use the last solution to improve accuracy by one order
    // DC solutions are stored in reversed order (last is the least accurate)
    BDFContainer *scForCorrection = _scBDF1;
    for(size_t iDC = 0; iDC < _scDC.size(); ++iDC)
    {
      // feInfo("Solving for DC #%d", _scDC.size() - 1 - iDC);
      feInfoCond(FE_VERBOSE > 0, "\t\t\tAdvancing DC%d/BDF1 solution to t = %1.4e", _order - iDC, _currentTime);
      DCBDFContainer *solDC = _scDC[_scDC.size() - 1 - iDC];

      solDC->rotate(_dt);
      solDC->computeBDFCoefficients(1, _deltaT);
      solDC->computeDCCoefficients(_deltaT);
      _currentSolution->setSolFromContainer(solDC);
      _currentSolution->initializeEssentialBC(_mesh, solDC);
      solDC->computeDCCorrection(*scForCorrection, _deltaT);

      // feInfo("Entering DC2 solve with containers:");
      // for(size_t i = 0; i < _scBDF1->_sol[0].size(); ++i) {
      //   // _scBDF1->_solDot[0][i] = nan("");
      //   feInfo("BDF1: s[0] = %+-1.3e s[1] = %+-1.2e [[d[0] = %+-1.2e]] d[1] = %+-1.2e  "
      //          "DC2: s[0] = %+-1.3e s[1] = %+-1.2e d[0] = %+-1.2e d[1] = %+-1.2e c = %+-1.3e ",
      //     _scBDF1->_sol[0][i], _scBDF1->_sol[1][i], _scBDF1->_solDot[0][i], _scBDF1->_solDot[1][i],
      //     _scDC[0]->_sol[0][i], _scDC[0]->_sol[1][i], _scDC[0]->_solDot[0][i], _scDC[0]->_solDot[1][i], _scDC[0]->_d[i]);
      //   if((int) i+1 == _linearSystem->getSystemSize())
      //   feInfo("=====================================================================================");
      // }

      feCheckReturn(solveNewtonRaphson(_linearSystem, _currentSolution, solDC, _NLoptions));

      // feInfo("Leaving DC2 solve with containers:");
      // for(size_t i = 0; i < _scBDF1->_sol[0].size(); ++i) {
      //   // _scBDF1->_solDot[0][i] = nan("");
      //   feInfo("BDF1: s[0] = %+-1.3e s[1] = %+-1.2e [[d[0] = %+-1.2e]] d[1] = %+-1.2e  "
      //          "DC2: s[0] = %+-1.3e s[1] = %+-1.2e d[0] = %+-1.2e d[1] = %+-1.2e c = %+-1.3e ",
      //     _scBDF1->_sol[0][i], _scBDF1->_sol[1][i], _scBDF1->_solDot[0][i], _scBDF1->_solDot[1][i],
      //     _scDC[0]->_sol[0][i], _scDC[0]->_sol[1][i], _scDC[0]->_solDot[0][i], _scDC[0]->_solDot[1][i], _scDC[0]->_d[i]);
      //   if((int) i+1 == _linearSystem->getSystemSize())
      //   feInfo("=====================================================================================");
      // }

      scForCorrection = solDC;
    }

    feCheckReturn(postProcessingAndVisualization(iStep+1, _currentStep+1, _t[0], _exportData));
  }

  return FE_STATUS_OK;
}

//
// Deferred correction DC3/BDF2 time integrator
//
DC3BDF2Integrator::DC3BDF2Integrator(feLinearSystem *linearSystem,
                                     feSolution *sol,
                                     feMesh *mesh,
                                     feNLSolverOptions &NLoptions,
                                     std::vector<feNorm*> &postProcessing,
                                     feExportData &exportData,
                                     double t0, double tEnd, int nTimeSteps)
  : TimeIntegrator(linearSystem, sol, mesh, NLoptions,
    postProcessing, exportData, 3, t0, tEnd, nTimeSteps)
{
  _integratorName = "DC3BDF2";
  _nSol = 3;
  _t.resize(_nSol, 0.);
  _deltaT.resize(_nSol, 0.);

  _sCBDF2 = new BDFContainer(_nSol, _currentTime, _currentSolution->getNumDOFs());
  _sC     = new DCBDFContainer(_order, 2, _nSol, _currentTime, _currentSolution->getNumDOFs());

  this->initialize();

  // Copy the initialized solution into the containers
  _sCBDF2->setCurrentSolution(_currentSolution);
  _sCBDF2->setCurrentSolutionDot(_currentSolution);

  _sC->setCurrentSolution(_currentSolution);
  _sC->setCurrentSolutionDot(_currentSolution);
}

feStatus DC3BDF2Integrator::restartFromContainer(const feSolutionContainer &container)
{
  UNUSED(container);
  return feErrorMsg(FE_STATUS_ERROR, "Reimplement restart for DC3BDF2Integrator");
  // if(_nSol != container.getNbSol()) {
  //   return feErrorMsg(FE_STATUS_ERROR, "Could not restart %s time integrator from solution container: "
  //     "integrator needs %d solution arrays but container stores %d",
  //     _integratorName.data(),
  //     _nSol,
  //     container.getNbSol());
  // }

  // int nDOF = _currentSolution->getNumDOFs();
  // if(nDOF != container.getNbDOFs()) {
  //   return feErrorMsg(FE_STATUS_ERROR, "Could not restart %s time integrator from solution container: "
  //     "current solution has %d DOF but container stores %d DOF per solution",
  //     _integratorName.data(),
  //     nDOF,
  //     container.getNbDOFs());
  // }

  // // Copy data from container :
  // const std::vector<double> &tC = container.getTime();
  // // const std::vector<double> &deltaTC = container.getTimeDifferences();

  // // Call custom copy constructor for the internal solution container
  // *_sC = container;

  // for(int i = 0; i < _nSol; ++i) {
  //   _t[i] = tC[i];
  // }
  // for(int i = 0; i < _nSol-1; ++i) {
  //   _deltaT[i] = _t[i] - _t[i+1];
  // }

  // return FE_STATUS_OK;
}

feStatus DC3BDF2Integrator::makeStep(double /* dt */)
{
  return feErrorMsg(FE_STATUS_ERROR, "Reimplement single makeStep");
}

feStatus DC3BDF2Integrator::makeSteps(int nSteps)
{
  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "DC3/BDF2 SOLVER:");
  feInfoCond(FE_VERBOSE > 0, "\t\tAdvancing from t = %1.4e to t = %1.4e in %d step(s)", 
    _currentTime, _currentTime + nSteps * _dt, nSteps);

  int start = 0;
  // int bdf2step = -1;
  // double bdf2time = -1.;

  //
  // Starting the DC3/BDF2 requires the BDF2 solutions up to t = t_2
  //
  if(_currentStep == 0)
  {
    start = 1;

    computePostProcessing(0, _t[0]);
    // Export initial solution
    feCheckReturn(this->writeVisualizationFile(0, _exportData));

    ///////////////////////////////////////////////////////////
    // // Compute BDF2 solution at time t = t_1 using BDF1 to start
    // double timeStepRatio = 0.1;
    // feCheckReturn(_scBDF2->startWithBDF1(timeStepRatio));

    BDF2Integrator solver(_linearSystem, _currentSolution, 
    _mesh, _NLoptions, _postProcessing, _exportData, _t0, _tEnd, _nTimeSteps);
    feCheckReturn(solver.makeSteps(2));
    *_sCBDF2 = solver.getSolutionContainer();
    // bdf2step = solver.getCurrentStep();
    // bdf2time = solver.getCurrentTime();

    std::vector<double> fakeDeltaT(3, _dt);
    _sC->computeDCCoefficients(fakeDeltaT);
    _sC->computeBDFCoefficients(2, fakeDeltaT);

    this->updateTime(0, _dt);
    _currentSolution->setC0(_sC->getC0());
    _currentSolution->setSolFromContainer(_sC);
    _currentSolution->initializeEssentialBC(_mesh, _sC);
    _sC->computeDCCorrection(*_sCBDF2, fakeDeltaT);
    feCheckReturn(solveNewtonRaphson(_linearSystem, _currentSolution, _sC, _NLoptions));

    // ///////////////////////////////////////////////////////////
    // // Start with exact solution
    // this->updateTime(0, _dt);

    // _sCBDF2->rotate(_dt);
    // _sC->rotate(_dt);

    // _currentSolution->initialize(_mesh);
    // _sCBDF2->setCurrentSolution(_currentSolution);
    // _sC->setCurrentSolution(_currentSolution);

    // // _sCBDF2->computeBDFCoefficients(2, _deltaT);
    // // _sC->computeBDFCoefficients(2, _deltaT);
    // // _sC->computeDCCoefficients(_deltaT);

    // // // Compute BDF2 time derivative
    // // _sCBDF2->computeSolTimeDerivative(_currentSolution, _linearSystem);
    // // _sC->computeSolTimeDerivative(_currentSolution, _linearSystem);

    // // Set the exact derivative at time t_1
    // for(size_t ii = 0; ii < _sC->_solDot[0].size(); ++ii) {
    //   _sCBDF2->_solDot[0][ii] = -sin(_currentTime); // <=== This is critical for order 3
    //   // _sCBDF2->_solDot[0][ii] = nan(""); // <=== This is critical for order 3
    // }
    // /////////////////////////////////////////////////////////////

    computePostProcessing(1, _t[0]);
    // Export solution for visualization
    _currentStep++;
    feCheckReturn(writeVisualizationFile(_currentStep, _exportData));
  }

  // At this point:
  // - the scBDF2 container contains the BDF2 solution at times 0,1,2
  // - the sc     container contains the DC3  solution at times 0,1

  for(int iStep = start; iStep < nSteps; ++iStep, ++_currentStep)
  {
    feInfoCond(FE_VERBOSE > 0, "\t\tDC3/BDF2 - Overall step = %d - Step %d/%d from t = %1.4e to t = %1.4e", 
      _currentStep+1, iStep+1, nSteps, _currentTime, _currentTime + _dt);

    this->updateTime(iStep, _dt);

    if(_currentStep >= 2)
      _sCBDF2->rotate(_dt);
    _sC->rotate(_dt);

    _sCBDF2->computeBDFCoefficients(2, _deltaT);
    _sC->computeBDFCoefficients(2, _deltaT);
    _sC->computeDCCoefficients(_deltaT);

    if(_currentStep >= 2) {
      //
      // Solve nonlinear problem for the BDF2 solution
      //
      _currentSolution->setC0(_sCBDF2->getC0());
      _currentSolution->setSolFromContainer(_sCBDF2);
      _currentSolution->initializeEssentialBC(_mesh, _sCBDF2);
      feCheckReturn(solveNewtonRaphson(_linearSystem, _currentSolution, _sCBDF2, _NLoptions));
    }

    //
    // Solve nonlinear problem for the DC3 solution
    //
    _currentSolution->setC0(_sC->getC0());
    // feInfo("c0 = %f", _sC->getBDFCoefficients()[0]);
    // feInfo("c1 = %f", _sC->getBDFCoefficients()[1]);
    // feInfo("c2 = %f", _sC->getBDFCoefficients()[2]);
    _currentSolution->setSolFromContainer(_sC);
    _currentSolution->initializeEssentialBC(_mesh, _sC);
    _sC->computeDCCorrection(*_sCBDF2, _deltaT);
    feCheckReturn(solveNewtonRaphson(_linearSystem, _currentSolution, _sC, _NLoptions));

    // _currentSolution->setSolFromContainer(_sCBDF2);

    // Compute post-processing indicators (L2 norms, etc.)
    computePostProcessing(iStep+1, _t[0]);

    // Export solution for visualization
    feCheckReturn(this->writeVisualizationFile(_currentStep+1, _exportData));
  }

  return FE_STATUS_OK;
}