#ifndef _FELINEARSYSTEM_
#define _FELINEARSYSTEM_

#include "feMessage.h"
#include "feNumber.h"
#include "feMesh.h"
#include "feBilinearForm.h"
#include "feSolution.h"
#include "feSolutionContainer.h"

/* Supported linear solvers */
typedef enum { MKLPARDISO, PETSC } linearSolverType;

// To call at the very beginning and the very end of the program
void petscInitialize(int argc, char **argv);
void petscFinalize();

// Create a linear system and perform safety checks.
// This is the recommended way of creating a linear system.
// Call within feCheck( ... ) to safely exit if a problem arises.
//
//        system: pointer to the linear system, initially undefined and assigned
//                during the call.
// bilinearForms: vector of (bi-)linear forms to assemble
//   numUnknowns: total number of unknowns (dimension of the system)

// argc and argv are provided for PETSc command line options
feStatus createLinearSystem(feLinearSystem *&system,
  linearSolverType type,
  std::vector<feBilinearForm *> bilinearForms,
  int numUnknowns,
  int argc = 0, char **argv = nullptr);

//
// Abstract class handling the linear system to be solved at each step
// of the nonlinear solver.
//
class feLinearSystem
{
protected:
  // The number of (bi-)linear forms that require to 
  // assemble a residual or a residual and a matrix
  int _numMatrixForms;
  int _numResidualForms;

  // The linear forms with a residual and a matrix
  std::vector<feBilinearForm *> _formMatrices;
  // The linear forms with only a residual
  std::vector<feBilinearForm *> _formResiduals;

  // Recompute the jacobian matrix?
  bool recomputeMatrix;

public:
  // Create an abstract linear system. Do not call directly, 
  // call the derived constructors instead.
  feLinearSystem(std::vector<feBilinearForm *> bilinearForms);
  virtual ~feLinearSystem() {}

  // Return the size m of the linear system (dimension of the square matrix m x m)
  virtual feInt getSystemSize() = 0;

  bool getRecomputeStatus() { return recomputeMatrix; }
  void setRecomputeStatus(bool status) { recomputeMatrix = status; }

  // Reset the matrix and/or the right-hand side
  virtual void setToZero() = 0;
  virtual void setMatrixToZero() = 0;
  virtual void setResidualToZero() = 0;

  // Assemble the matrix and/or the right-hand side
  virtual void assemble(feSolution *sol) = 0;
  virtual void assembleMatrices(feSolution *sol) = 0;
  virtual void assembleResiduals(feSolution *sol) = 0;

  // Solve the linear system Ax-b
  //
  //       normDx: norm of the solution vector x (correction in the Newton-Rapshon iteration)
  // normResidual: norm of the RHS b (residual in the Newton-Raphson iteration)
  //      normAxb: norm of the residual Ax-b
  //        nIter: number of iteration used to solve (0 for direct solver)
  virtual void solve(double *normDx, double *normResidual, double *normAxb, int *nIter) = 0;

  // Apply the Newton-Raphson correction u_new = u_old + du
  virtual void correctSolution(feSolution *sol) = 0;

  // Assign the RHS to the residual vector of the solution container.
  // Used in DC (deferred correction) time integration schemes to
  // improve the time derivative.
  virtual void assignResidualToDCResidual(feSolutionContainer *solContainer) = 0;

  // Perform RHS += coeff * d
  // Used in DC schemes.
  virtual void applyCorrectionToResidual(double coeff, std::vector<double> &d) = 0;

  // Print the matrix to the console or PETSc viewer
  virtual void viewMatrix() = 0;
  // Print the RHS to the console
  virtual void printRHS() = 0;
};

#endif