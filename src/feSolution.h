#ifndef _FESOLUTION_
#define _FESOLUTION_

#include "feMesh.h"
#include "feSpace.h"
#include "feNumber.h"
#include "feSolutionContainer.h"

class feSolutionContainer;

class feSolution {
protected:
  int _dim;

  std::vector<double> _sol;
  std::vector<double> _dsoldt;

  const std::vector<feSpace *> &_space;
  const std::vector<feSpace *> &_essBC;

  double _c0;
  double _tn; // The current time

  double _t0;
  double _t1;
  int _nTimeSteps;
  double _dt;

public:
  feSolution(feMesh *mesh, const std::vector<feSpace *> &space, const std::vector<feSpace *> &essBC, feMetaNumber *metaNumber);
  feSolution(std::string solutionFile);
  ~feSolution() {}

  void initializeTemporalSolution(double t0, double t1, int nTimeSteps);
  void initializeUnknowns(feMesh *mesh, feMetaNumber *metaNumber);
  void initializeEssentialBC(feMesh *mesh, feMetaNumber *metaNumber,
                             feSolutionContainer *solContainer = nullptr);

  double getC0() { return _c0; }
  void setC0(double c0) { _c0 = c0; }
  double getCurrentTime() { return _tn; }
  void setCurrentTime(double t) { _tn = t; }
  double getTimeStep() { return _dt; }
  double getNbTimeSteps() { return _nTimeSteps; }
  double getInitialTime() { return _t0; }
  double getFinalTime() { return _t1; }

  std::vector<double> getSolutionCopy() { return _sol; }
  std::vector<double> &getSolutionReference() { return _sol; }
  double getSolAtDOF(int iDOF) { return _sol[iDOF]; }
  void setSolAtDOF(int iDOF, double val) { _sol[iDOF] = val; }
  void incrementSolAtDOF(int iDOF, double val) { _sol[iDOF] += val; }
  void setSolFromContainer(feSolutionContainer *solContainer, int iSol = 0);
  double getSolDotAtDOF(int iDOF) { return _dsoldt[iDOF]; }
  void setSolDotAtDOF(int iDOF, double val) { _dsoldt[iDOF] = val; }
  void setSolDotToZero();

  void printSol(std::string file = "");
};

#endif