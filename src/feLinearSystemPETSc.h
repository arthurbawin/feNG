#ifndef _FELINEARSYSTEMPETSC_
#define _FELINEARSYSTEMPETSC_

#include "feLinearSystem.h"

#if defined(HAVE_PETSC)
#include "petscksp.h"
#endif

class feLinearSystemPETSc : public feLinearSystem
{
protected:
  int _argc;
  char **_argv;
#if defined(HAVE_PETSC)
  PetscInt _nInc;
  Mat _A;
  Vec _res; // True residual
  Vec _dx;
  Vec _linSysRes; // Residual of the linear system : A*dx - res

  KSP ksp;
  PC preconditioner;
#endif

public:
  feLinearSystemPETSc(int argc, char **argv, std::vector<feBilinearForm *> bilinearForms, int numUnknowns);
  ~feLinearSystemPETSc();

  feInt getSystemSize() {return (feInt) _nInc; };

  void setToZero();
  void setMatrixToZero();
  void setResidualToZero();
  void assemble(feSolution *sol);
  void assembleMatrices(feSolution *sol);
  void assembleResiduals(feSolution *sol);
  void solve(double *normDx, double *normResidual, double *normAxb, int *nIter);
  void correctSolution(feSolution *sol);
  void assignResidualToDCResidual(feSolutionContainer *solContainer);
  void applyCorrectionToResidual(double coeff, std::vector<double> &d);
  void viewMatrix();
  void printRHS();

private:
  void initialize();
  void finalize();
  
};

#endif