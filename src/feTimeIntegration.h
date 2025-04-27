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
  DC2BDF1,
  DC3BDF1,
  DC3BDF2,
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
  virtual feSolutionContainer getSolutionContainer() const = 0;

  void setQuietStart(bool flag) { _quietStart = flag; }
  void setExportInitialCondition(bool flag) { _exportInitialCondition = flag; }

protected:
  void initialize();
  void updateTime(const int iStep, const double dt);
  void resetTime();
  void computePostProcessing(const int iStep, const double t);
  feStatus writeVisualizationFile(const int iStep, const feExportData &data);
  feStatus postProcessingAndVisualization(const int localStep, const int globalStep, const double t, const feExportData &data);
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
  feSolutionContainer getSolutionContainer() const { return *_sC; };
  feStatus restartFromContainer(const feSolutionContainer &container);
};

class BDF1Integrator : public TimeIntegrator
{
public:
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
  feSolutionContainer getSolutionContainer() const { return *_sC; };
  feStatus restartFromContainer(const feSolutionContainer &container);
};

class BDF2Integrator : public TimeIntegrator
{
public:
  // BDF2Starter _startMethod = BDF2Starter::BDF1;
  BDF2Starter _startMethod = BDF2Starter::InitialCondition;
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
  feSolutionContainer getSolutionContainer() const { return *_sC; };
  feStatus restartFromContainer(const feSolutionContainer &container);
private:
  feStatus startWithBDF1(const double timeStepRatio);
};

//
// Deferred correction implicit methods based on BDF1
//
class DCBDF1Integrator : public TimeIntegrator
{
public:
  // A container for the BDF1 solution
  BDFContainer *_scBDF1;
  // Containers for the DC solutions of increasing order.
  // A DCn/BDF1 method needs the DC solutions of order 2 to n. 
  std::vector<DCBDFContainer*> _scDC;

public:
  DCBDF1Integrator(const int order,
                   feLinearSystem *linearSystem,
                   feSolution *sol,
                   feMesh *mesh,
                   feNLSolverOptions &NLoptions,
                   std::vector<feNorm*> &postProcessing,
                   feExportData &exportData,
                   double t0, double tEnd, int nTimeSteps);
  ~DCBDF1Integrator()
  {
    delete _scBDF1;
    for(auto *ptr : _scDC) delete ptr;
  };
  feStatus makeStep(double dt);
  feStatus makeSteps(int nSteps);
  feSolutionContainer getSolutionContainer() const { return *_scDC[0]; };
  feStatus restartFromContainer(const feSolutionContainer &container);
};

class DC3BDF2Integrator : public TimeIntegrator
{
public:
  BDF2Starter _startMethod = BDF2Starter::BDF1;
  BDFContainer *_sCBDF2;
  DCBDFContainer *_sC;

public:
  DC3BDF2Integrator(feLinearSystem *linearSystem,
                    feSolution *sol,
                    feMesh *mesh,
                    feNLSolverOptions &NLoptions,
                    std::vector<feNorm*> &postProcessing,
                    feExportData &exportData,
                    double t0, double tEnd, int nTimeSteps);
  ~DC3BDF2Integrator()
  {
    delete _sCBDF2;
    delete _sC;
  };
  feStatus makeStep(double dt);
  feStatus makeSteps(int nSteps);
  feSolutionContainer getSolutionContainer() const { return *_sC; };
  feStatus restartFromContainer(const feSolutionContainer &container);
private:
  feStatus startWithBDF1(const double timeStepRatio);
};

#endif