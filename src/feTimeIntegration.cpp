
#include "feTimeIntegration.h"
#include "feMessage.h"

extern int FE_VERBOSE;

feStatus createTimeIntegrator2(TimeIntegrator2 *&solver, timeIntegratorScheme2 scheme,
                              feNLSolverOptions NLoptions, feLinearSystem *system,
                              feSolution *solution, feMesh *mesh,
                              std::vector<feNorm*> &norms, feExportData exportData, double tBegin,
                              double tEnd, int nTimeSteps)
{
  UNUSED(mesh, tBegin, tEnd, nTimeSteps);
  switch(scheme) {
    case timeIntegratorScheme2::STATIONARY:
      solver = new StationaryIntegrator(system, solution, mesh, NLoptions, norms, exportData);
      break;
    case timeIntegratorScheme2::BDF1:
      solver = new BDF1Integrator(system, solution, mesh, NLoptions, norms, exportData, tBegin, tEnd, nTimeSteps);
      break;
    case timeIntegratorScheme2::BDF2:
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

TimeIntegrator2::TimeIntegrator2(feLinearSystem *linearSystem,
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
};

//
// Initialize the solution array
//
void TimeIntegrator2::initialize()
{
  _t[0] = _currentTime;
  _currentSolution->setCurrentTime(_currentTime);
  _currentSolution->initializeUnknowns(_mesh);
  _currentSolution->initializeEssentialBC(_mesh);
}

feStatus TimeIntegrator2::restartFromContainer(const feSolutionContainer &container)
{
  if(_nSol != container.getNbSol()) {
    return feErrorMsg(FE_STATUS_ERROR, "Could not restart %s time integrator from solution container: "
      "integrator needs %d solution arrays but container stores %d",
      _integratorName.data(),
      _nSol,
      container.getNbSol());
  }

  const std::vector<double> &tC = container.getTime();
  const std::vector<double> &deltaTC = container.getTimeDifferences();

  for(int i = 0; i < _nSol; ++i) {
    _t[i] = tC[i];
    _deltaT[i] = deltaTC[i];
  }

  return FE_STATUS_OK;
}

void TimeIntegrator2::updateTime(const int iStep, const double dt)
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

  feInfo("t = [%f,%f,%f], deltaT = [%f,%f,%f], ct = %f",
    _t[0], _t[1], _t[2],
    _deltaT[0], _deltaT[1], _deltaT[2],
    _currentTime);
}

feStatus TimeIntegrator2::writeVisualizationFile(const int iStep, const feExportData &data)
{
  if(data.exporter == nullptr)
    return FE_STATUS_OK;

  if((iStep % data.exportEveryNSteps) == 0) {
    std::string fileName = data.fileNameRoot + "_" + _integratorName + "_" + std::to_string(iStep) + ".vtk";
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
  : TimeIntegrator2(linearSystem, sol, mesh, NLoptions, postProcessing, exportData, 0)
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
  : TimeIntegrator2(linearSystem, sol, mesh, NLoptions,
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

feStatus BDF1Integrator::makeStep(double /* dt */)
{
  return feErrorMsg(FE_STATUS_ERROR, "Reimplement single makeStep");
}

feStatus BDF1Integrator::makeSteps(int nSteps)
{
  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "BDF1 SOLVER:");
  feInfoCond(FE_VERBOSE > 0, "\t\tAdvancing %d steps from t = %f to t = %f", 
    nSteps, _currentTime, _currentTime + nSteps * _dt);

  // Export initial solution
  feStatus s = this->writeVisualizationFile(0, _exportData);
  if(s != FE_STATUS_OK) return s;

  for(int iStep = 0; iStep < nSteps; ++iStep)
  {
    feInfoCond(FE_VERBOSE > 0, "\t\tBDF1 - Step %d/%d from t = %1.4e to t = %1.4e", 
      iStep+1, nSteps, _currentTime, _currentTime + _dt);

    this->updateTime(iStep, _dt);

    _sC->rotateWithoutTime();
    _sC->computeBDFCoefficients(_order, _deltaT);
    _currentSolution->setC0(_sC->getC0());
    _currentSolution->initializeEssentialBC(_mesh, _sC);

    // Solve nonlinear problem
    s = solveNewtonRaphson(_linearSystem, _currentSolution, _sC, _NLoptions);
    if(s == FE_STATUS_ERROR) { return s; }

    // Compute post-processing indicators (L2 norms, etc.)
    for(size_t i = 0; i < _postProcessing.size(); ++i) {
      _postProcessingData[i][0] = _postProcessing[i]->compute();
    }

    // Export solution for visualization
    s = this->writeVisualizationFile(iStep + 1, _exportData);
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
  : TimeIntegrator2(linearSystem, sol, mesh, NLoptions,
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

feStatus BDF2Integrator::makeStep(double /* dt */)
{
  return feErrorMsg(FE_STATUS_ERROR, "Reimplement single makeStep");
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
  feInfoCond(FE_VERBOSE > 0, "\t\tAdvancing %d steps from t = %f to t = %f", 
    nSteps, _currentTime, _currentTime + nSteps * _dt);

  if(_currentStep == 0)
  {
    // Export initial solution
    feStatus s = this->writeVisualizationFile(0, _exportData);
    if(s != FE_STATUS_OK) return s;

    // Start the BDF2 by computing initial conditions at time n-1 and n-2
    this->updateTime(0, _dt);
    _sC->rotate(_dt);
    _currentSolution->initializeUnknowns(_mesh);
    _currentSolution->initializeEssentialBC(_mesh);
    _sC->setCurrentSolution(_currentSolution);
    _currentStep++;
    nSteps--;
  }

  for(int iStep = 0; iStep < nSteps; ++iStep, ++_currentStep)
  {
    feInfoCond(FE_VERBOSE > 0, "\t\tBDF2 - Overall step = %d - Step %d/%d from t = %1.4e to t = %1.4e", 
      _currentStep, iStep+1, nSteps, _currentTime, _currentTime + _dt);

    this->updateTime(iStep, _dt);
    _sC->rotate(_dt);
    _sC->computeBDFCoefficients(_order, _deltaT);
    _currentSolution->setC0(_sC->getC0());
    _currentSolution->setCurrentTime(_currentTime);
    _currentSolution->initializeEssentialBC(_mesh, _sC);

    // Solve nonlinear problem
    feStatus s = solveNewtonRaphson(_linearSystem, _currentSolution, _sC, _NLoptions);
    if(s == FE_STATUS_ERROR) { return s; }

    // Compute post-processing indicators (L2 norms, etc.)
    for(size_t i = 0; i < _postProcessing.size(); ++i) {
      _postProcessingData[i][0] = _postProcessing[i]->compute();
    }

    // Export solution for visualization
    if(_exportData.exporter != nullptr && (_currentStep % _exportData.exportEveryNSteps) == 0) {
      std::string fileName = _exportData.fileNameRoot + "_BDF2_" + std::to_string(_currentStep) + ".vtk";
      s = _exportData.exporter->writeStep(fileName);
      if(s != FE_STATUS_OK) return s;
    }
  }

  return FE_STATUS_OK;
}