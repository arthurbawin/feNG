#ifndef _FESOLUTION_
#define _FESOLUTION_

#include "feMesh.h"
#include "feSpace.h"
#include "feNumber.h"
#include "feSolutionContainer.h"
#include "feSolutionContainerV2.h"
#include "omp.h"

class feSolutionContainer;
class feSolutionContainerV2;

class feSolution
{
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
  feSolution(feMesh *mesh, const std::vector<feSpace *> &space, const std::vector<feSpace *> &essBC,
             feMetaNumber *metaNumber);
  feSolution(std::string solutionFile);
  ~feSolution() {}

  void initializeTemporalSolution(double t0, double t1, int nTimeSteps);
  void initializeUnknowns(feMesh *mesh, feMetaNumber *metaNumber);
  void initializeEssentialBC(feMesh *mesh, feMetaNumber *metaNumber, feSolutionContainer *solContainer = nullptr);
  void initializeContainerEssentialBC(feMesh *mesh, feMetaNumber *metaNumber, feSolutionContainerV2 *solContainer = nullptr);
  
  double getC0() { return _c0; }
  void setC0(double c0) { _c0 = c0; }
  double getCurrentTime() { return _tn; }
  void setCurrentTime(double t) { _tn = t; }
  double getTimeStep() { return _dt; }
  double getNbTimeSteps() { return _nTimeSteps; }
  double getInitialTime() { return _t0; }
  double getFinalTime() { return _t1; }

  std::vector<double> getSolutionCopy() { return _sol; }
  std::vector<double> getSolutionDotCopy() { return _dsoldt; }
  std::vector<double> &getSolutionReference() { return _sol; }
  std::vector<double> &getSolutionReferenceDot() { return _dsoldt; }
  double getSolAtDOF(int iDOF) { return _sol[iDOF]; }
  void setSolAtDOF(int iDOF, double val) { _sol[iDOF] = val; }
  void incrementSolAtDOF(int iDOF, double val) { _sol[iDOF] += val; }
  void incrementSolDotAtDOF(int iDOF, double val) { _dsoldt[iDOF] += val; }
  void setSolFromContainer(feSolutionContainer *solContainer, int iSol = 0);
  void setSolFromContainer(feSolutionContainerV2 *solContainer, int iSol = 0);
  void setSolDotFromContainer(feSolutionContainerV2 *solContainer, int iSolDot = 0);
  double getSolDotAtDOF(int iDOF) { return _dsoldt[iDOF]; }
  void setSolDotAtDOF(int iDOF, double val) { _dsoldt[iDOF] = val; }
  void setSolDotToZero();

  void copySpace(feMesh *mesh, feMetaNumber *metaNumber, feSpace *s1, feSpace *s2);

  void printSol(std::string file = "");
  void displaySolDot();

};

void myFun();

#endif