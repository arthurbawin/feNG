#ifndef _FELINEARSYSTEMPETSC_
#define _FELINEARSYSTEMPETSC_

#include "feLinearSystem.h"
#include "feCompressedRowStorage.h"

#if defined(HAVE_PETSC)
#include "petscksp.h"
#endif

//
// PETSc iterative solver
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

  // Krylov solver (default is GMRES)
  KSP ksp;
  // Preconditioner (default is ILU(0))
  PC preconditioner;
#endif

  feEZCompressedRowStorage *_EZCRS;

public:
  feLinearSystemPETSc(int argc, char **argv, std::vector<feBilinearForm *> bilinearForms, int numUnknowns);
  ~feLinearSystemPETSc();

  feInt getSystemSize() {return (feInt) _nInc; };

  // See doc in feLinearSystem.h
  void setToZero();
  void setMatrixToZero();
  void setResidualToZero();
  void assemble(feSolution *sol);
  void assembleMatrices(feSolution *sol);
  void assembleResiduals(feSolution *sol);
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

#endif