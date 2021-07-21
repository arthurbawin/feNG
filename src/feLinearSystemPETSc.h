#ifndef _FELINEARSYSTEMPETSC_
#define _FELINEARSYSTEMPETSC_

#include "feLinearSystem.h"

#include "petscksp.h"

class feLinearSystemPETSc : public feLinearSystem {

protected:
  int _argc;
  char** _argv;
  PetscInt _nInc;
  PetscInt _nDofs;
  Mat _A;
  Vec _res; // True residual
  Vec _dx;
  Vec _linSysRes; // Residual of the linear system : A*dx - res

  KSP ksp;
  PC preconditioner;
public:
  feLinearSystemPETSc(int argc, char** argv, std::vector<feBilinearForm*> &formMatrices, 
    std::vector<feBilinearForm*> &formResiduals, feMetaNumber *metaNumber, feMesh *mesh);
	virtual ~feLinearSystemPETSc();

  virtual void initialize();
  virtual void finalize();
  virtual void setToZero();
  virtual void setMatrixToZero();
  virtual void setResidualToZero();
  virtual void assembleMatrices(feSolution *sol);
  virtual void assembleResiduals(feSolution *sol);
  virtual void assemble(feSolution *sol);
  virtual void solve(double *normDx, double *normResidual, double *normAxb, int *nIter);
  virtual void correctSolution(feSolution *sol);
  virtual void assignResidualToDCResidual(feSolutionContainer *solContainer);
  virtual void applyCorrectionToResidual(double coeff, std::vector<double> &d);
  virtual void viewMatrix();
  virtual void printResidual();
};

// To call at the very beginning and the very end of the program
void petscInitialize(int argc, char** argv);
void petscFinalize();

#endif