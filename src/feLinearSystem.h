#ifndef _FELINEARSYSTEM_
#define _FELINEARSYSTEM_

#include "feBilinearForm.h"
#include "feCompressedRowStorage.h"
#include "feMesh.h"
#include "feMessage.h"
#include "feNumber.h"
#include "feSolution.h"
#include "feSolutionContainer.h"

#if defined(HAVE_PETSC)
#  include "petscksp.h"
#endif

/* Supported linear solvers */
typedef enum
{
  MKLPARDISO,
  PETSC,
  PETSC_MUMPS
} linearSolverType;

// Create a linear system and perform safety checks.
// This is the recommended way of creating a linear system.
// Call within feCheck( ... ) to safely exit if a problem arises.
//
//        system: pointer to the linear system, initially undefined and assigned
//                during the call.
// bilinearForms: vector of (bi-)linear forms to assemble
//   numUnknowns: total number of unknowns (dimension of the system)

feStatus createLinearSystem(feLinearSystem                    *&system,
                            const linearSolverType              type,
                            const std::vector<feBilinearForm *> bilinearForms,
                            const feMetaNumber                 *numbering,
                            const int ownershipSplit = -1);

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

  // A copy of the DOF numbering
  const feMetaNumber *_numbering;

  // Recompute the jacobian matrix at current Newton-Raphson step?
  bool _recomputeMatrix = true;

  // Reorder the matrix?
  bool _permute = false;

  // Options for iterative solver:
  double _rel_tol  = 1e-8;
  double _abs_tol  = 1e-14;
  double _div_tol  = 1e6;
  int    _max_iter = 1e4;

  // Display and export options:
  bool _displayMatrixInConsole = false;
  bool _displayMatrixInWindow  = false;
  bool _displayRHSInConsole    = false;
  bool _exportMatrixMatlab     = false;
  bool _exportRHSMatlab        = false;

public:
  // Create an abstract linear system. Do not call directly,
  // call the derived constructors instead.
  feLinearSystem(const std::vector<feBilinearForm *> bilinearForms,
                 const feMetaNumber                 *numbering);
  virtual ~feLinearSystem();

  // Return the size m of the linear system (dimension of the square matrix m x
  // m)
  virtual feInt getSystemSize() const = 0;

  void setAbsoluteTol(double absTol) { _abs_tol = absTol; };
  void setRelativeTol(double relTol) { _rel_tol = relTol; };
  void setDivergenceTol(double divTol) { _div_tol = divTol; };
  void setMaxIter(int maxIter) { _max_iter = maxIter; };

  void setDisplayMatrixInConsole(bool flag) { _displayMatrixInConsole = flag; };
  void setDisplayMatrixInWindow(bool flag) { _displayMatrixInWindow = flag; };
  void setDisplayRHSInConsole(bool flag) { _displayRHSInConsole = flag; };
  void setExportMatrixMatlab(bool flag) { _exportMatrixMatlab = flag; };
  void setExportRHSMatlab(bool flag) { _exportRHSMatlab = flag; };

  bool getRecomputeStatus() const { return _recomputeMatrix; }
  void setRecomputeStatus(bool status) { _recomputeMatrix = status; }

  bool getReorderingStatus() const { return _permute; }
  void setReorderingStatus(bool status) { _permute = status; }

  virtual void getRHSMaxNorm(double *norm) const      = 0;
  virtual void getResidualMaxNorm(double *norm) const = 0;

  // Reset the matrix and/or the right-hand side
  virtual void setToZero()         = 0;
  virtual void setMatrixToZero()   = 0;
  virtual void setResidualToZero() = 0;

  // Assemble the matrix and/or the right-hand side
  virtual void assemble(const feSolution *sol,
                        const bool assembleOnlyTransientMatrices = false) = 0;
  virtual void
               assembleMatrices(const feSolution *sol,
                                const bool        assembleOnlyTransientMatrices = false) = 0;
  virtual void assembleResiduals(const feSolution *sol) = 0;

  virtual void constrainEssentialComponents(const feSolution *sol) = 0;

  virtual void applyPeriodicity() = 0;

  virtual void permute() = 0;

  // Solve the linear system Ax = b, which is J*du = -NL(u) (or -R(u)) in a
  // Newton step
  //
  //       normDx: norm of the solution vector x (= du, correction in the
  //       Newton-Rapshon iteration)
  // normResidual: norm of the RHS b (= -NL(u)//R(u), residual in the
  // Newton-Raphson iteration)
  //      normAxb: norm of the residual Ax-b
  //        nIter: number of iteration used to solve (0 for direct solver)
  virtual bool
  solve(double *normDx, double *normResidual, double *normAxb, int *nIter) = 0;

  // Apply the Newton-Raphson correction u_new = u_old + du
  virtual void correctSolution(feSolution *sol,
                               const bool  correctSolutionDot = false) = 0;

  // Assign the RHS to the residual vector of the solution container.
  // Used in DC (deferred correction) time integration schemes to
  // improve the time derivative.
  virtual void
  assignResidualToDCResidual(feSolutionContainer *solContainer) = 0;

  // Perform RHS += coeff * d
  // Used in DC schemes.
  virtual void applyCorrectionToResidual(double               coeff,
                                         std::vector<double> &d) = 0;

  // Print the matrix to the console or PETSc viewer
  virtual void viewMatrix() const = 0;
  // Print the RHS to the console
  virtual void viewRHS() const = 0;
  // Print the current solution vector (Newton residual)
  virtual void viewResidual() const = 0;
  // Export raw matrix to file
  virtual void writeMatrix(const std::string fileName,
                           const double      time = 0.) = 0;
  // Export raw RHS to file
  virtual void writeRHS(const std::string fileName, const double time = 0.) = 0;
  // Export raw residual du to file
  virtual void writeResidual(const std::string fileName,
                             const double      time = 0.) = 0;
};

//
// PETSc solvers
//
class feLinearSystemPETSc : public feLinearSystem
{
protected:
#if defined(HAVE_PETSC)
  // Number of unknown DOFs (size of the system)
  PetscInt _nInc;
  // Global finite element matrix
  Mat _A;
  // Reordered global matrix and row/col indices
  Mat _Ap;
  IS  _rowMap;
  IS  _colMap;

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

  // Row ownership arrays (bounds and number of rows owned by each MPI process)
  int *_ownedLowerBounds, *_ownedUpperBounds, *_numOwnedRows;
  // Owned solution buffer
  // Currently each process knows the whole solution vector,
  // but only updates its owned portion, then the vector is Allgatherv'ed
  double *_ownedSolution;

  // To solve with MUMPS without runtime options:
  bool _solveWithMUMPS = false;
  Mat  _factoredMatrix;
  // Selected MUMPS options, see MUMPS user guide
  int _mumps_icntl7 =
    2; // Sequential ordering (if icntl28 = 1 or only 1 MPI proc)
#  if defined(PETSC_HAVE_PTSCOTCH) || defined(PETSC_HAVE_PARMETIS)
  int _mumps_icntl28 = 2; // Sequential (1) or parallel (2) ordering
#  else
  int _mumps_icntl28 = 1;
#  endif
  int _mumps_icntl29 = 1; // Parallel ordering (if icntl28 = 2)
#endif

public:
  feLinearSystemPETSc(const std::vector<feBilinearForm *> bilinearForms,
                      const feMetaNumber                 *numbering,
                      const linearSolverType              type);
  ~feLinearSystemPETSc();

  feInt getSystemSize() const
  {
#if defined(HAVE_PETSC)
    return (feInt)_nInc;
#else
    return 0;
#endif
  };

  void getRHSMaxNorm(double *norm) const;
  void getResidualMaxNorm(double *norm) const;

  // See doc in feLinearSystem.h
  void setToZero();
  void setMatrixToZero();
  void setResidualToZero();
  void assemble(const feSolution *sol,
                const bool        assembleOnlyTransientMatrices);
  void assembleMatrices(const feSolution *sol,
                        const bool        assembleOnlyTransientMatrices);
  void assembleResiduals(const feSolution *sol);
  void constrainEssentialComponents(const feSolution *sol);
  // void applyPeriodicity();
  void permute();
  bool solve(double *normSolution,
             double *normRHS,
             double *normResidualAxMinusb,
             int    *nIter);
  void correctSolution(feSolution *sol, const bool correctSolutionDot);
  void assignResidualToDCResidual(feSolutionContainer *solContainer);
  void applyCorrectionToResidual(double coeff, std::vector<double> &d);
  void viewMatrix() const;
  void viewRHS() const;
  void viewResidual() const;
  void writeMatrix(const std::string /*fileName*/, const double /* time */);
  void writeRHS(const std::string /*fileName*/, const double /* time */);
  void writeResidual(const std::string /*fileName*/, const double /* time */);

private:
  void initialize();
  void initializeSequential();
  int  initializeMPI();
};

//
// MKL Pardiso direct solver
//
#if defined(HAVE_MKL)
class feLinearSystemMklPardiso : public feLinearSystem
{
protected:
  // Number of unknown DOFs (size of the system)
  PardisoInt _nInc;

  // FE matrix stored in CRS (compressed row storage) format
  PardisoInt  _nnz;
  PardisoInt *_mat_ia;     // dimension _nInc+1, IA
  PardisoInt *_mat_ja;     // dimension _nnz, JA
  double     *_mat_values; // dimension _nnz

  // Global finite element RHS
  double *_rhs;
  // Solution vector du
  double *du;

  // Compute symbolic factorization (matrix reordering)
  bool symbolicFactorization = true;

  // Row ownership arrays (bounds and number of rows owned by each MPI process)
  // For single process run, these are the data for the whole matrix
  // stored on the process.
  // Ownership can be determined as follows:
  //  - -1 : Matrix is stored sequentially on process 0
  //  -  0 : Even split of rows between processes according to
  //  PetscSplitOwnership()'s formula:
  //         n = N/size + ((N % size) > rank)
  //         Independent of the number of nnz per row, thus
  //         does not guarantee a load-balancing distribution.
  //  -  1 : Even split of nonzero entries in the matrix (how is the RHS split
  //  then?)
  //         Causes overlapping rows in general between consecutive processes,
  //         as the split will probably happen in the middle of a row
  //  -  2 : Even split of nonzero entries without overlapping rows
  int         _ownershipSplit = -1;
  PardisoInt  _numOwnedRows, _numOwnedNNZ;
  PardisoInt *_ownedLowerBounds, *_ownedUpperBounds, *_numOwnedRowsOnAllProc;
  // Owned solution buffer
  double *_ownedSolution;

  // Data to speedup DOF constraining:
  // For each row, store the number and positions of its occurences in JA
  bool                         _initializeConstrainData = true;
  std::vector<feInt>           _rowsToConstrain;
  std::vector<size_t>          _numOccurencesInJa;
  std::vector<std::set<feInt>> _posOccurencesInJa;

  // Pardiso options
  void       *PT[64];
  PardisoInt  MYTPE;
  PardisoInt  IPARM[64];
  double      DPARM[64];
  PardisoInt  MAXFCT;
  PardisoInt  MNUM;
  PardisoInt  MTYPE;
  PardisoInt  SOLVER;
  PardisoInt  PHASE;
  PardisoInt  NRHS;
  PardisoInt  ERROR;
  PardisoInt  MSGLVL;
  PardisoInt *IDUM;
  double      DDUM;
  PardisoInt  IPIVOT;

public:
  feLinearSystemMklPardiso(const std::vector<feBilinearForm *> bilinearForms,
                           const feMetaNumber                 *numbering,
                           const int                           ownershipSplit);
  ~feLinearSystemMklPardiso();

  void setPivot(int pivot);
  void setPardisoMsglvlHigh() { MSGLVL = 1; };
  void setPardisoMsglvlLow() { MSGLVL = 0; };

  feInt getSystemSize() const { return _nInc; };

  void getRHSMaxNorm(double *norm) const;
  void getResidualMaxNorm(double *norm) const;

  void setToZero();
  void setMatrixToZero();
  void setResidualToZero();
  void assemble(const feSolution *sol,
                const bool        assembleOnlyTransientMatrices);
  void assembleMatrices(const feSolution *sol,
                        const bool        assembleOnlyTransientMatrices);
  void assembleResiduals(const feSolution *sol);
  void constrainEssentialComponents(const feSolution *sol);
  void applyPeriodicity();
  void permute();
  bool solve(double *normDx, double *normResidual, double *normAxb, int *nIter);
  void assignResidualToDCResidual(feSolutionContainer *solContainer);
  void applyCorrectionToResidual(double coeff, std::vector<double> &d);
  void correctSolution(feSolution *sol, const bool correctSolutionDot = false);
  void viewMatrix() const;
  void viewRHS() const;
  void viewResidual() const;
  void writeMatrix(const std::string fileName, const double time = 0.);
  void writeRHS(const std::string fileName, const double time = 0.);
  void writeResidual(const std::string fileName, const double time = 0.);

private:
  void initialize(void);
  void mklSymbolicFactorization(void);
  void mklFactorization(void);
  void mklSolve(void);
  // void mklSolveWithPhase(const int phase);
};
#endif // HAVE_MKL

#endif