#ifndef _FELINEARSYSTEM_
#define _FELINEARSYSTEM_

#include "feNumber.h"
#include "feMesh.h"
#include "feBilinearForm.h"
#include "feSolution.h"
#include "feSolutionContainer.h"

class feLinearSystem {
protected:
  std::vector<feBilinearForm *> &_formMatrices;
  std::vector<feBilinearForm *> &_formResiduals;
  feMetaNumber *_metaNumber;
  feMesh *_mesh;

  bool recomputeMatrix;

public:
  feLinearSystem(std::vector<feBilinearForm *> &formMatrices,
                 std::vector<feBilinearForm *> &formResiduals, feMetaNumber *metaNumber,
                 feMesh *mesh)
    : _formMatrices(formMatrices), _formResiduals(formResiduals), _metaNumber(metaNumber),
      _mesh(mesh), recomputeMatrix(false){};
  virtual ~feLinearSystem() {}

  bool getRecomputeStatus() { return recomputeMatrix; }
  void setRecomputeStatus(bool status) { recomputeMatrix = status; }

  virtual void initialize(){};
  // virtual void finalize();
  virtual void setToZero() = 0;
  virtual void setMatrixToZero(){};
  virtual void setResidualToZero(){};
  virtual void assembleMatrices(feSolution *sol){};
  virtual void assembleResiduals(feSolution *sol){};
  virtual void assemble(feSolution *sol){};
  virtual void solve(double *normDx, double *normResidual, double *normAxb, int *nIter){};
  virtual void correctSolution(feSolution *sol){};
  virtual void assignResidualToDCResidual(feSolutionContainer *solContainer){};
  virtual void applyCorrectionToResidual(double coeff, std::vector<double> &d){};
  virtual void viewMatrix(){};
  virtual void printResidual(){};
};

#endif