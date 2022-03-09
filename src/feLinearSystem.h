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
  std::vector<std::vector <feBilinearForm *>> _formMatrices;
  std::vector<std::vector <feBilinearForm *>> _formResiduals;
  feMetaNumber *_metaNumber;
  feMesh *_mesh;

  int _nbThreads;

  bool recomputeMatrix;

public:
  feLinearSystem(std::vector<feBilinearForm *> bilinearForms, feMetaNumber *metaNumber,
                 feMesh *mesh)
    : _metaNumber(metaNumber), _mesh(mesh), recomputeMatrix(false)
  {
    
    _nbThreads=1;
    #if defined(HAVE_OMP)
     _nbThreads=omp_get_max_threads();
    #endif

    _formResiduals.resize(_nbThreads);
    _formMatrices.resize(_nbThreads);
    for (int i=0; i<_nbThreads;++i){
      for(feBilinearForm *f : bilinearForms) {
        feBilinearForm *fCpy = new feBilinearForm(*f);
        feInfo("Copie #%d est a l'adresse", i);
        std::cout<<fCpy<<std::endl;
        _formResiduals[i].push_back(fCpy);
        if(f->hasMatrix()) _formMatrices[i].push_back(fCpy);
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
  virtual int getNbThreads(){return _nbThreads;};
  virtual std::vector<std::vector <feBilinearForm *>> getformMatrices(){return _formMatrices;};
  virtual std::vector<std::vector <feBilinearForm *>> getformResiduals(){return _formResiduals;};
};

feStatus createLinearSystem(feLinearSystem *&system, linearSolverType type,
                            std::vector<feBilinearForm *> bilinearForms, feMetaNumber *metaNumber,
                            feMesh *mesh, int argc = 0, char **argv = nullptr);

#endif