#ifndef _FELINEARSYSTEM_
#define _FELINEARSYSTEM_

#include "feMessage.h"
#include "feNumber.h"
#include "feMesh.h"
#include "feBilinearForm.h"
#include "feSolution.h"
#include "feSolutionContainer.h"
#include "feCompressedRowStorage.h"

#if defined(HAVE_PETSC)
#include "petscksp.h"
#endif

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
feStatus createLinearSystem(feLinearSystem *&system, linearSolverType type,
                            std::vector<feBilinearForm *> bilinearForms, int numUnknowns,
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

  // Recompute the jacobian matrix at current Newton-Raphson step?
  bool recomputeMatrix;

  // Options for iterative solver:
  double _rel_tol = 1e-8;
  double _abs_tol = 1e-14;
  double _div_tol = 1e6;
  int _max_iter = 5e3;

  // Display and export options:
  bool _displayMatrixInConsole = false;
  bool _displayMatrixInWindow = false;
  bool _displayRHSInConsole = false;
  bool _exportMatrixMatlab = false;
  bool _exportRHSMatlab = false;

public:
  // Create an abstract linear system. Do not call directly,
  // call the derived constructors instead.
  feLinearSystem(std::vector<feBilinearForm *> bilinearForms);
  virtual ~feLinearSystem() {}

  // Return the size m of the linear system (dimension of the square matrix m x m)
  virtual feInt getSystemSize() = 0;

  void setAbsoluteTol(double absTol) { _abs_tol = absTol; };
  void setRelativeTol(double relTol) { _rel_tol = relTol; };
  void setDivergenceTol(double divTol) { _div_tol = divTol; };
  void setMaxIter(int maxIter) { _max_iter = maxIter; };

  void setDisplayMatrixInConsole(bool flag) { _displayMatrixInConsole = flag; };
  void setDisplayMatrixInWindow(bool flag) { _displayMatrixInWindow = flag; };
  void setDisplayRHSInConsole(bool flag) { _displayRHSInConsole = flag; };
  void setExportMatrixMatlab(bool flag) { _exportMatrixMatlab = flag; };
  void setExportRHSMatlab(bool flag) { _exportRHSMatlab = flag; };

  bool getRecomputeStatus() { return recomputeMatrix; }
  void setRecomputeStatus(bool status) { recomputeMatrix = status; }

  virtual void getResidualMaxNorm(double *norm) = 0;

  // Reset the matrix and/or the right-hand side
  virtual void setToZero() = 0;
  virtual void setMatrixToZero() = 0;
  virtual void setResidualToZero() = 0;

  // Assemble the matrix and/or the right-hand side
  virtual void assemble(feSolution *sol) = 0;
  virtual void assembleMatrices(feSolution *sol) = 0;
  virtual void assembleResiduals(feSolution *sol) = 0;

  virtual void constraintEssentialComponents(feSolution *sol) = 0;

  // Solve the linear system Ax = b, which is J*du = -NL(u) (or -R(u)) in a Newton step
  //
  //       normDx: norm of the solution vector x (= du, correction in the Newton-Rapshon iteration)
  // normResidual: norm of the RHS b (= -NL(u)//R(u), residual in the Newton-Raphson iteration)
  //      normAxb: norm of the residual Ax-b
  //        nIter: number of iteration used to solve (0 for direct solver)
  virtual bool solve(double *normDx, double *normResidual, double *normAxb, int *nIter) = 0;

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
  virtual void viewRHS() = 0;
};

//
// PETSc iterative solvers
//
class feLinearSystemPETSc : public feLinearSystem
{
protected:
  int _argc;
  char **_argv;

#if defined(HAVE_PETSC)
  // Number of unknown DOFs (size of the system)
  PetscInt _nInc;
  // Global finite element matrix
  Mat _A;
  // Reordered global matrix and row/col indices
  Mat _Ap;
  IS _rowMap;
  IS _colMap;

  // Global finite element RHS
  Vec _rhs;
  // Solution vector du
  Vec _du;
  // Residual of the linear system : A*du - RHS
  Vec _linSysRes;
  // Values of constrained DOF (usually 0)
  Vec _constrainedDOFValue;

  // Krylov solver (default is GMRES)
  KSP ksp;
  // Preconditioner (default is ILU(0))
  PC preconditioner;
#endif

  feEZCompressedRowStorage *_EZCRS;

public:
  feLinearSystemPETSc(int argc, char **argv, std::vector<feBilinearForm *> bilinearForms,
                      int numUnknowns);
  ~feLinearSystemPETSc();

  feInt getSystemSize() { return (feInt)_nInc; };

  void getResidualMaxNorm(double *norm);

  // See doc in feLinearSystem.h
  void setToZero();
  void setMatrixToZero();
  void setResidualToZero();
  void assemble(feSolution *sol);
  void assembleMatrices(feSolution *sol);
  void assembleResiduals(feSolution *sol);
  void constraintEssentialComponents(feSolution *sol);
  bool solve(double *normSolution, double *normRHS, double *normResidualAxMinusb, int *nIter);
  void correctSolution(feSolution *sol);
  void assignResidualToDCResidual(feSolutionContainer *solContainer);
  void applyCorrectionToResidual(double coeff, std::vector<double> &d);
  void viewMatrix();
  void viewRHS();

private:
  void initialize();
  void finalize();
};

//
// MKL Pardiso direct solver
//
#if defined(HAVE_MKL)
class feLinearSystemMklPardiso : public feLinearSystem
{
protected:
  feCompressedRowStorageMklPardiso *crsMklPardiso;
  feMKLPardisoInt matrixOrder;
  double *du;
  double *residu;
  feInt iparm12 = 1; // Modification de IPARM[12]
  bool symbolicFactorization = true;

  void *PT[64];
  feMKLPardisoInt MYTPE;
  feMKLPardisoInt IPARM[64];
  double DPARM[64];

  feMKLPardisoInt MAXFCT;
  feMKLPardisoInt MNUM;
  feMKLPardisoInt MTYPE;
  feMKLPardisoInt SOLVER;
  feMKLPardisoInt PHASE;
  feMKLPardisoInt N;
  feMKLPardisoInt NRHS;
  feMKLPardisoInt ERROR;
  feMKLPardisoInt MSGLVL;
  feMKLPardisoInt IDUM;
  double DDUM;
  feMKLPardisoInt IPIVOT;

  feMKLPardisoInt nz;
  feMKLPardisoInt *Ap; // dimension ordre+1, IA
  feMKLPardisoInt *Aj; // dimension nz, JA
  double *Ax; // dimension nz

public:
  feLinearSystemMklPardiso(std::vector<feBilinearForm *> bilinearForms, int numUnknowns);
  ~feLinearSystemMklPardiso();

  // ====================================================================
  // Pour modifier les paramètres de Pardiso
  // pivot : valeur de l'exposant du nombre réel 10^(-pivot)
  //         valeur par défaut pivot := 13
  // ====================================================================
  void setPivot(int pivot);
  void setPardisoMsglvlHigh() { MSGLVL = 1; };
  void setPardisoMsglvlLow() { MSGLVL = 0; };
  bool getSymbolicFactorizationStatus() { return symbolicFactorization; }
  void setSymbolicFactorizationStatus(bool status)
  {
    symbolicFactorization = status;
    if(symbolicFactorization) recomputeMatrix = true;
  }
  void toggleSymbolicFactorizationStatus()
  {
    symbolicFactorization = !symbolicFactorization;
    if(symbolicFactorization) recomputeMatrix = true;
  }

  void setToZero();
  void setMatrixToZero();
  void setResidualToZero();
  void assemble(feSolution *sol);
  void assembleMatrices(feSolution *sol);
  void assembleResiduals(feSolution *sol);
  void constraintEssentialComponents(feSolution *sol){};
  bool solve(double *normDx, double *normResidual, double *normAxb, int *nIter);
  void assignResidualToDCResidual(feSolutionContainer *solContainer);
  void applyCorrectionToResidual(double coeff, std::vector<double> &d);
  void correctSolution(feSolution *sol);
  void correctSolution(double *sol);
  void viewMatrix(){};
  void viewRHS(){};

private:
  void mklSymbolicFactorization(void);
  void mklFactorization(void);
  void mklSolve(void);
};
#endif // HAVE_MKL

#endif