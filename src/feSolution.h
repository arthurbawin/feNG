#ifndef _FESOLUTION_
#define _FESOLUTION_

#include "feMesh.h"
#include "feSpace.h"
#include "feNumber.h"
#include "feSolutionContainer.h"
#include "omp.h"

class feSolutionContainer;

//
// A wrapper containing the FE solution (DOF values) (and 
// time integration info that should be moved to the solver).
//
class feSolution
{
protected:
  
  // The number of DOFs
  int _nDOF;

  // The FE solution and its time derivative in DOF number order.
  std::vector<double> _sol;
  std::vector<double> _dsoldt;

  // Convenience references to the FE spaces
  const std::vector<feSpace *> &_spaces;
  const std::vector<feSpace *> &_essentialSpaces;

  // Time integration info will be moved to the solver later
  // First coefficient of the BDF approximation, multiplying the
  // newest solution.
  double _c0;
  // Current time
  double _tn;
  // Initial and final time of integration, (number of) time steps
  double _t0;
  double _t1;
  int _nTimeSteps;
  double _dt;

public:
  // Create a solution wrapper
  feSolution(int numDOF, const std::vector<feSpace *> &spaces, const std::vector<feSpace *> &essentialSpaces);
  // Create 
  feSolution(std::string solutionFile);
  ~feSolution() {}

  double getC0() { return _c0; }
  void setC0(double c0) { _c0 = c0; }
  double getCurrentTime() const { return _tn; }
  void setCurrentTime(double t) { _tn = t; }
  double getTimeStep() { return _dt; }
  double getNbTimeSteps() { return _nTimeSteps; }
  double getInitialTime() { return _t0; }
  double getFinalTime() { return _t1; }

  // Return a copy of the solution vector
  std::vector<double> getSolutionCopy() { return _sol; }
  // Return a reference to the solution vector
  std::vector<double> &getSolutionReference() { return _sol; }
  // Return a reference to the time derivative of the solution
  std::vector<double> &getSolutionReferenceDot() { return _dsoldt; }

  double getSolAtDOF(int iDOF) { return _sol[iDOF]; }
  void getSolAtDOF(const std::vector<feInt> &addressing, std::vector<double> &sol) const;
  void setSolAtDOF(int iDOF, double val) { _sol[iDOF] = val; }
  void incrementSolAtDOF(int iDOF, double val) { _sol[iDOF] += val; }
  void setSolFromContainer(feSolutionContainer *solContainer, int iSol = 0);
  double getSolDotAtDOF(int iDOF) { return _dsoldt[iDOF]; }
  void setSolDotAtDOF(int iDOF, double val) { _dsoldt[iDOF] = val; }

  // Reset the time derivative of the solution
  void setSolDotToZero();

  // Initialize the unknown and essential DOFs
  void initializeUnknowns(feMesh *mesh);
  void initializeEssentialBC(feMesh *mesh, feSolutionContainer *solContainer = nullptr);

  // Should be moved elsewhere
  void initializeTemporalSolution(double t0, double t1, int nTimeSteps);

  void copySpace(feMesh *mesh, feSpace *s1, feSpace *s2);

  void printSol(std::string file = "");
};

#endif