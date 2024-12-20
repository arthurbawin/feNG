#ifndef _FESOLVER_
#define _FESOLVER_

#include "feMessage.h"
#include "feMesh.h"
#include "feLinearSystem.h"
#include "feBilinearForm.h"
#include "feNumber.h"
#include "feNorm.h"
#include "feSpace.h"
#include "feExporter.h"

typedef struct feTol {
  double tolDx;
  double tolResidual;
  double maxIter;
  // int maxIterWithoutRecompute = 5;
  // double ratioUntilRecomputing= 1e-2;
} feTolerances;

/* Supported time integrators schemes */
typedef enum { EIGENSOLVER, STATIONARY, BDF1, BDF2, DC2F, DC3, DC3F } timeIntegratorScheme;

class TimeIntegrator
{
protected:
  feTolerances _tol;
  feMetaNumber *_metaNumber;
  feLinearSystem *_linearSystem;
  feSolution *_sol;
  std::vector<feNorm *> &_norms;
  feMesh *_mesh;

  feExportData _exportData;

  feSolutionContainer *_solutionContainer;

  double _t0;
  double _tEnd;
  double _dt;
  int _nTimeSteps;

  double _tCurrent;
  int _currentStep;

  bool _shouldInitialize = true;

  std::vector<std::vector<double> > _normL2;
  std::vector<std::vector<double> > _IntBoundary;

  std::string _CodeIni; // Code pour connaitre le mode d'initialisation    |"SolExacte"->avec la
                        // solution exacte|  |"BDF1/DCF"->en partant d'un BDF1 puis DCF|

public:
  TimeIntegrator(feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem,
                 feSolution *sol, std::vector<feNorm *> &norms, feMesh *mesh,
                 feExportData exportData, double t0 = 0., double tEnd = 0., int nTimeSteps = 1,
                 std::string CodeIni = "")
    : _tol(tol), _metaNumber(metaNumber), _linearSystem(linearSystem), _sol(sol), _norms(norms),
      _mesh(mesh), _exportData(exportData), _t0(t0), _tEnd(tEnd), _nTimeSteps(nTimeSteps),
      _tCurrent(t0), _currentStep(0), _CodeIni(CodeIni)
  {
    _dt = (tEnd - t0) / (double)nTimeSteps;
    _sol->initializeTemporalSolution(_t0, _tEnd, _nTimeSteps);
  };
  virtual ~TimeIntegrator()
  {
    delete _solutionContainer;
  };

  // Bad: the solution container may be deleted
  feSolutionContainer *getSolutionContainer() { return _solutionContainer; }
  // Better?
  void copySolutionContainerInto(feSolutionContainer &targetContainer)
  {
    targetContainer.copy(*(this->_solutionContainer));
  }

  void setMetaNumber(feMetaNumber *metaNumber) { _metaNumber = metaNumber; }
  void setLinearSystem(feLinearSystem *linearSystem) { _linearSystem = linearSystem; }
  void setSolution(feSolution *sol) { _sol = sol; }
  void setMesh(feMesh *mesh) { _mesh = mesh; }
  std::vector<double> &getNorm(int iNorm) { return _normL2[iNorm]; };

  int getCurrentStep() { return _currentStep; };
  void setCurrentStep(int step) { _currentStep = step; };

  virtual feStatus makeStep() = 0;
  virtual feStatus makeSteps(int nSteps) = 0;
};

class StationarySolver : public TimeIntegrator
{
protected:
public:
  StationarySolver(feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem,
                   feSolution *sol, std::vector<feNorm *> &norms, feMesh *mesh,
                   feExportData exportData);
  ~StationarySolver() {}

  double getNorm(int iNorm) { return _normL2[iNorm][0]; };

  feStatus makeStep();
  feStatus makeSteps(int nSteps);
};

class BDF2Solver : public TimeIntegrator
{
protected:
public:
  BDF2Solver(feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem,
             feSolution *sol, std::vector<feNorm *> &norms, feMesh *mesh, feExportData exportData,
             double t0, double tEnd, int nTimeSteps, std::string CodeIni);
  ~BDF2Solver() {}

  std::vector<double> &getNorm(int iNorm) { return _normL2[iNorm]; };

  feStatus makeStep();
  feStatus makeSteps(int nSteps);
};

class BDF1Solver : public TimeIntegrator
{
protected:
public:
  BDF1Solver(feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem,
             feSolution *sol, std::vector<feNorm *> &norms, feMesh *mesh, feExportData exportData,
             double t0, double tEnd, int nTimeSteps);
  ~BDF1Solver() {}

  std::vector<double> &getNorm(int iNorm) { return _normL2[iNorm]; };

  feStatus makeStep();
  feStatus makeSteps(int nSteps);
};

class DC2FSolver : public TimeIntegrator
{
protected:
  feSolutionContainer *_solutionContainerBDF1;
  std::vector<std::vector<double> > _normL2BDF1;

public:
  DC2FSolver(feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem,
             feSolution *sol, std::vector<feNorm *> &norms, feMesh *mesh, feExportData exportData,
             double t0, double tEnd, int nTimeSteps);
  ~DC2FSolver() { delete _solutionContainerBDF1; }

  std::vector<double> &getNorm(int iNorm) { return _normL2[iNorm]; };

  feStatus makeStep() { return FE_STATUS_OK; };
  feStatus makeSteps(int nSteps);
};

class DC3FSolver : public TimeIntegrator
{
protected:
  feSolutionContainer *_solutionContainerBDF1;
  feSolutionContainer *_solutionContainerDC2F;

public:
  DC3FSolver(feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem,
             feSolution *sol, std::vector<feNorm *> &norms, feMesh *mesh, feExportData exportData,
             double t0, double tEnd, int nTimeSteps);
  ~DC3FSolver()
  {
    delete _solutionContainerBDF1;
    delete _solutionContainerDC2F;
  }

  std::vector<double> &getNorm(int iNorm) { return _normL2[iNorm]; };

  feStatus makeStep() { return FE_STATUS_OK; };
  feStatus makeSteps(int nSteps);
};

class DC3FSolver_centered : public TimeIntegrator
{
protected:
  feSolutionContainer *_solutionContainerBDF1;
  feSolutionContainer *_solutionContainerDC2F;

public:
  DC3FSolver_centered(feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem,
                      feSolution *sol, std::vector<feNorm *> &norms, feMesh *mesh,
                      feExportData exportData, double t0, double tEnd, int nTimeSteps);
  ~DC3FSolver_centered()
  {
    delete _solutionContainerBDF1;
    delete _solutionContainerDC2F;
  }

  std::vector<double> &getNorm(int iNorm) { return _normL2[iNorm]; };

  feStatus makeStep() { return FE_STATUS_OK; };
  feStatus makeSteps(int nSteps);
};

class DC3Solver : public TimeIntegrator
{
protected:
  feSolutionContainer *_solutionContainerBDF2;
  feSolutionContainer *_solutionContainerDC3F;

public:
  DC3Solver(feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem,
            feSolution *sol, std::vector<feNorm *> &norms, feMesh *mesh, feExportData exportData,
            double t0, double tEnd, int nTimeSteps, std::string CodeIni);
  ~DC3Solver()
  {
    delete _solutionContainerBDF2;
    delete _solutionContainerDC3F;
  }

  std::vector<double> &getNorm(int iNorm) { return _normL2[iNorm]; };
  std::vector<double> &getInt(int iNorm) { return _normL2[iNorm]; };

  feStatus makeStep() { return FE_STATUS_OK; };
  feStatus makeSteps(int nSteps);
};

feStatus createTimeIntegrator(TimeIntegrator *&solver, timeIntegratorScheme scheme,
                              feTolerances tolerances, feLinearSystem *system,
                              feMetaNumber *metaNumber, feSolution *solution, feMesh *mesh,
                              std::vector<feNorm *> &norms, feExportData exportData = {},
                              double tBegin = 0., double tEnd = 0., int nTimeSteps = 1,
                              std::string initializationCode = "");

#endif