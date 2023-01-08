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

class feLinearSystem
{
protected:
  int _numMatrixForms;
  int _numResidualForms;
  std::vector<feBilinearForm *> _formMatrices;
  std::vector<feBilinearForm *> _formResiduals;
  feMetaNumber *_metaNumber;
  feMesh *_mesh;

  bool recomputeMatrix;

public:
  feLinearSystem(std::vector<feBilinearForm *> bilinearForms, feMetaNumber *metaNumber,
                 feMesh *mesh)
    : _metaNumber(metaNumber), _mesh(mesh), recomputeMatrix(false)
  {
    _numMatrixForms = 0;
    _numResidualForms = bilinearForms.size();

#if defined(HAVE_OMP)
    int nThreads = omp_get_max_threads();
#else
    int nThreads = 1;
#endif

    for(int i = 0; i < nThreads; ++i) {
      for(feBilinearForm *f : bilinearForms) {
        #if defined(HAVE_OMP)
            feBilinearForm *fCpy = new feBilinearForm(*f);
            _formResiduals.push_back(fCpy);
            if(f->hasMatrix()) _formMatrices.push_back(fCpy);
        #else
            _formResiduals.push_back(f);
            if(f->hasMatrix()) _formMatrices.push_back(f);
        #endif
        if(f->hasMatrix() && i == 0) _numMatrixForms++;
      }
    }
  };

  virtual ~feLinearSystem() {}

  bool getRecomputeStatus() { return recomputeMatrix; }
  void setRecomputeStatus(bool status) { recomputeMatrix = status; }

  virtual void initialize(){};
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
  virtual std::vector<feBilinearForm *> getformMatrices() { return _formMatrices; };
  virtual std::vector<feBilinearForm *> getformResiduals() { return _formResiduals; };
};

feStatus createLinearSystem(feLinearSystem *&system, linearSolverType type,
                            std::vector<feSpace *> allFESpaces,
                            std::vector<feBilinearForm *> bilinearForms, feMetaNumber *metaNumber,
                            feMesh *mesh, int argc = 0, char **argv = nullptr);

#endif