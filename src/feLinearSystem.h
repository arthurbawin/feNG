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
  std::vector<feBilinearForm *> _formMatrices;
  std::vector<feBilinearForm *> _formResiduals;
  feMetaNumber *_metaNumber;
  feMesh *_mesh;

  int _nbThreads;
  std::vector<feBilinearForm *> _formMatricesOMP;
  std::vector<feBilinearForm *> _formResidualsOMP;

  bool recomputeMatrix;

public:
  feLinearSystem(std::vector<feBilinearForm *> bilinearForms, feMetaNumber *metaNumber,
                 feMesh *mesh)
    : _metaNumber(metaNumber), _mesh(mesh), recomputeMatrix(false)
  {
// <<<<<<< HEAD
//     #if defined(HAVE_OMP)
//    _nbThreads = omp_get_max_threads();
//    #else
//    _nbThreads = 1;
//    #endif
//     _formResiduals.resize(_nbThreads);
//     _formMatrices.resize(_nbThreads);
//     for(int i=0; i<_nbThreads;++i){
// =======
    for(feBilinearForm *f : bilinearForms) {
      feBilinearForm *fCpy = new feBilinearForm(*f);
      _formResiduals.push_back(fCpy);
      if(f->hasMatrix()) _formMatrices.push_back(fCpy);
    }

    _nbThreads = 1;
#if defined(HAVE_OMP)
    _nbThreads = omp_get_max_threads();
    feInfo("Nombre Max Threads : %d", _nbThreads);

    for(int i = 0; i < _nbThreads; ++i) {
// >>>>>>> origin/baptiste
      for(feBilinearForm *f : bilinearForms) {
        feBilinearForm *fCpy = new feBilinearForm(*f);
        _formResidualsOMP.push_back(fCpy);
        if(f->hasMatrix()) _formMatricesOMP.push_back(fCpy);
      }
    }
#endif

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
  virtual int getNbThreads() { return _nbThreads; };
  virtual std::vector<feBilinearForm *> getformMatrices() { return _formMatrices; };
  virtual std::vector<feBilinearForm *> getformResiduals() { return _formResiduals; };
};

feStatus createLinearSystem(feLinearSystem *&system, linearSolverType type,
                            std::vector<feBilinearForm *> bilinearForms, feMetaNumber *metaNumber,
                            feMesh *mesh, int argc = 0, char **argv = nullptr);

#endif