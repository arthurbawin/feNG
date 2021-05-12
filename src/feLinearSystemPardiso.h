#ifndef _FELINEARSYSTEMPARDISO_
#define _FELINEARSYSTEMPARDISO_

#include "feLinearSystem.h"

class feLinearSystemPardiso : public feLinearSystem {

protected:

public:
  feLinearSystemPardiso(std::vector<feBilinearForm*> &formMatrices, 
    std::vector<feBilinearForm*> &formResiduals, feMetaNumber *metaNumber, feMesh *mesh)
  : feLinearSystem(formMatrices, formResiduals, metaNumber, mesh)
  {
  }
	virtual ~feLinearSystemPardiso() {
  }

  virtual void initialize();
  virtual void finalize();
  virtual void setToZero();
  virtual void setMatrixToZero();
  virtual void setResidualToZero();
  virtual void assembleMatrices(feSolution *sol);
  virtual void assembleResiduals(feSolution *sol);
  virtual void assemble(feSolution *sol);
  virtual void solve(double *normDx, double *normResidual);
  virtual void correctSolution(feSolution *sol);
  virtual void assignResidualToDCResidual(feSolutionContainer *solContainer);
  virtual void applyCorrectionToResidual(feSolutionContainer *solContainer);
  virtual void viewMatrix();
  virtual void printResidual();
};

#endif