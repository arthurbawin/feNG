#ifndef _FETIMEINTEGRATION_
#define _FETIMEINTEGRATION_

#include "feNonLinearSolver.h"
#include "feLinearSystem.h"
#include "feSolution.h"
#include "feNorm.h"
#include "feExporter.h"
#include "feSolutionContainer.h"

//
// Supported time integrators schemes
//
enum class timeIntegratorScheme {
  STATIONARY,
  BDF1,
  BDF2,
  DC2F,
  DC3,
  DC3F,
  EIGENSOLVER
};

enum class BDF2Starter {
  BDF1,
  InitialCondition
};

class TimeIntegrator;

feStatus createTimeIntegrator(TimeIntegrator *&solver, timeIntegratorScheme scheme,
                              feNLSolverOptions tolerances, feLinearSystem *system,
                              feSolution *solution, feMesh *mesh,
                              std::vector<feNorm *> &norms, feExportData exportData = {},
                              double tBegin = 0., double tEnd = 0., int nTimeSteps = 1);

//
// Abstract class for a time integration scheme
//
class TimeIntegrator
{
protected:
  feLinearSystem       *_linearSystem;
  feSolution           *_currentSolution;
  feMesh               *_mesh;
  feNLSolverOptions     _NLoptions;
  std::vector<feNorm*> &_postProcessing;
  feExportData          _exportData;

  // Number of stored previous solutions to approximate time derivative
  // (= number of steps of the method)
  int _nSol;
  int _order;
  double _t0;
  double _tEnd;
  double _dt;
  int _nTimeSteps;
  double _currentTime;
  int _currentStep;

  // Time of the N last steps (for multistep methods)
  std::vector<double> _t;
  std::vector<double> _deltaT;
  std::vector<double> _timeHistory;

  // Name of the time integration scheme to write visualization files
  std::string _integratorName;

  bool _quietStart = false;
  bool _exportInitialCondition = true;

  // Scalar-valued post-processing operators computed at each step
  // First vector is the time
  std::vector<std::vector<double>> _postProcessingData;

protected:
  TimeIntegrator(feLinearSystem *linearSystem,
                 feSolution *sol,
                 feMesh *mesh,
                 feNLSolverOptions &NLoptions,
                 std::vector<feNorm*> &postProcessing,
                 feExportData &exportData,
                 int order,
                 double t0 = 0.,
                 double tEnd = 0.,
                 int nTimeSteps = 1);
public:
  virtual ~TimeIntegrator(){};
  virtual feStatus makeStep(double dt) = 0;
  virtual feStatus makeSteps(int nSteps) = 0;
  int getCurrentStep() { return _currentStep; };
  void setCurrentStep(int step) { _currentStep = step; };

  const std::vector<double> &getTime() const { return _t; }
  double getCurrentTime() const { return _t[0]; }

  // Set the previous times and time differences from existing container,
  // e.g. when continuing integration from a previous integrator
  virtual feStatus restartFromContainer(const feSolutionContainer &container) = 0;

  std::vector<std::vector<double>> getPostProcessingData(){ return _postProcessingData; };
  virtual feSolutionContainer getSolutionContainer() = 0;

  void setQuietStart(bool flag) { _quietStart = flag; }
  void setExportInitialCondition(bool flag) { _exportInitialCondition = flag; }

protected:
  void initialize();
  void updateTime(const int iStep, const double dt);
  void computePostProcessing(const int iStep, const double t);
  feStatus writeVisualizationFile(const int iStep, const feExportData &data);
};

class StationaryIntegrator : public TimeIntegrator
{
protected:
  // A specialized solutionContainer
  BDFContainer *_sC;

public:
  StationaryIntegrator(feLinearSystem *linearSystem,
                       feSolution *sol,
                       feMesh *mesh,
                       feNLSolverOptions &NLoptions,
                       std::vector<feNorm*> &postProcessing,
                       feExportData &exportData);
  ~StationaryIntegrator()
  {
    delete _sC;
  };
  feStatus makeStep(double dt = 0.);
  feStatus makeSteps(int nSteps);
  feSolutionContainer getSolutionContainer() { return *_sC; };
  feStatus restartFromContainer(const feSolutionContainer &container);
};

class BDF1Integrator : public TimeIntegrator
{
protected:
  BDFContainer *_sC;

public:
  BDF1Integrator(feLinearSystem *linearSystem,
                 feSolution *sol,
                 feMesh *mesh,
                 feNLSolverOptions &NLoptions,
                 std::vector<feNorm*> &postProcessing,
                 feExportData &exportData,
                 double t0, double tEnd, int nTimeSteps);
  ~BDF1Integrator()
  {
    delete _sC;
  };
  feStatus makeStep(double dt);
  feStatus makeSteps(int nSteps);
  feSolutionContainer getSolutionContainer() { return *_sC; };
  feStatus restartFromContainer(const feSolutionContainer &container);
};

class BDF2Integrator : public TimeIntegrator
{
protected:
  BDF2Starter _startMethod = BDF2Starter::BDF1;
  BDFContainer *_sC;

public:
  BDF2Integrator(feLinearSystem *linearSystem,
                 feSolution *sol,
                 feMesh *mesh,
                 feNLSolverOptions &NLoptions,
                 std::vector<feNorm*> &postProcessing,
                 feExportData &exportData,
                 double t0, double tEnd, int nTimeSteps);
  ~BDF2Integrator()
  {
    delete _sC;
  };
  feStatus makeStep(double dt);
  feStatus makeSteps(int nSteps);
  feSolutionContainer getSolutionContainer() { return *_sC; };
  feStatus restartFromContainer(const feSolutionContainer &container);
private:
  feStatus startWithBDF1(const double timeStepRatio);
};

#endif