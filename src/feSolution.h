#ifndef _FESOLUTION_
#define _FESOLUTION_

#include "feMesh.h"
#include "feSpace.h"
#include "feNumber.h"
#include "feSolutionContainer.h"
#if defined(HAVE_OMP)
#include "omp.h"
#endif

class feSolutionContainer;

//
// A wrapper containing the FE solution (DOF values) (and
// time integration info that should be moved to the solver).
//
class feSolution
{
protected:
  // The number of DOFs (size of the arrays)
  int _nDOF;

  // The FE solution and its time derivative in DOF number order.
  std::vector<double> _sol;
  std::vector<double> _dsoldt;

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
  // Copy of the FE spaces ptrs
  std::vector<feSpace *> _spaces;
  std::vector<feSpace *> _essentialSpaces;

public:
  // Create a solution wrapper holding numDOF degrees of freedom and their time derivative.
  feSolution(int numDOF, const std::vector<feSpace *> &spaces,
             const std::vector<feSpace *> &essentialSpaces);
  feSolution(const feSolutionContainer &container,
             const int solutionIndex,
             const std::vector<feSpace*> &spaces,
             const std::vector<feSpace*> &essentialSpaces);
  // Create a solution wrapper from a solution file. See printSol function for parsing.
  feSolution(){};
  feSolution(std::string solutionFile);
  ~feSolution() {}

  // Trivial getters and setters
  int getNumDOFs() { return _nDOF; }
  double getC0() const { return _c0; }
  void setC0(double c0) { _c0 = c0; }
  double getCurrentTime() const { return _tn; }
  void setCurrentTime(double t) { _tn = t; }
  double getTimeStep() const { return _dt; }
  double getNbTimeSteps() const { return _nTimeSteps; }
  double getInitialTime() const { return _t0; }
  double getFinalTime() const { return _t1; }

  const std::vector<double> &getSolution() const { return _sol; }
  const std::vector<double> &getSolutionDot() const { return _dsoldt; }
  std::vector<double> getSolutionCopy() const { return _sol; }
  std::vector<double> getSolutionDotCopy() const { return _dsoldt; }
  std::vector<double> &getSolutionReference() { return _sol; }
  std::vector<double> &getSolutionReferenceDot() { return _dsoldt; }

  // Get the solution (or time derivative) at dof with number 'iDOF'
  double getSolAtDOF(int iDOF) const { return _sol[iDOF]; }
  double getSolDotAtDOF(int iDOF) const { return _dsoldt[iDOF]; }
  // Fill 'sol' with the solution at dofs with numbers stored in 'addressing'.
  // Sizes of the vector arguments must match.
  void getSolAtDOF(const std::vector<feInt> &addressing, std::vector<double> &sol) const;
  void setSolAtDOF(int iDOF, double val) { _sol[iDOF] = val; }
  void setSolDotAtDOF(int iDOF, double val) { _dsoldt[iDOF] = val; }

  void incrementSolAtDOF(int iDOF, double val) { _sol[iDOF] += val; }

  // Set the solution array to the iSol-th solution stored in solContainer
  void setSolFromContainer(feSolutionContainer *solContainer, int iSol = 0);

  // Reset the time derivative of the solution
  void setSolDotToZero();

  // Initialize the unknown and essential DOFs using each FE space's feFunction
  // and dof initialization mode.
  void initialize(feMesh *mesh);
  void initializeUnknowns(feMesh *mesh);
  void initializeEssentialBC(feMesh *mesh, feSolutionContainer *solContainer = nullptr);

  // Should be moved elsewhere
  void initializeTemporalSolution(double t0, double t1, int nTimeSteps);

  void copySpace(feMesh *mesh, feSpace *s1, feSpace *s2);
  feStatus addConstantTimesSpace(feMesh *mesh, 
  const double coeff, feSpace *sourceSpace, feSpace *targetSpace);
  feStatus addSquaredNormOfVectorSpace(feMesh *mesh, 
    const double coeff, feSpace *sourceSpace, feSpace *targetSpace);

  void printSol(std::string file = "");
};

#endif