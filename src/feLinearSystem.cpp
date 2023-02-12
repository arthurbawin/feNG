#include "feLinearSystem.h"
#include "feLinearSystemMklPardiso.h"
#include "feLinearSystemPETSc.h"

extern int FE_VERBOSE;

extern bool petscWasInitialized;

feStatus createLinearSystem(feLinearSystem *&system,
  linearSolverType type, 
  std::vector<feBilinearForm *> bilinearForms,
  int numUnknowns,
  int argc, char **argv)
{
  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "LINEAR SYSTEM:");

  for(auto *form : bilinearForms){
    if(form == nullptr){
      return feErrorMsg(FE_STATUS_ERROR, "Null pointer in vector of bilinear forms, maybe you forgot to initialize it.");
    }
  }

  if(numUnknowns == 0)
    return feErrorMsg(FE_STATUS_ERROR, "0 unknowns : attempting to create a linear system of size 0.");

  switch(type) {
    case MKLPARDISO:
#if defined(HAVE_MKL)
      system = new feLinearSystemMklPardiso(bilinearForms, numUnknowns);
      break;
#else
      return feErrorMsg(FE_STATUS_ERROR,
                        "feNG must be compiled with Intel MKL to solve with MKL Pardiso.");
#endif
    case PETSC:
#if defined(HAVE_PETSC)
      if(argc == 0 || argv == nullptr) {
        return feErrorMsg(FE_STATUS_ERROR,
                          "Please provide argc and argv to create a PETSc linear system.");
      }
      if(!petscWasInitialized) {
        return feErrorMsg(FE_STATUS_ERROR,
                          "PETSc was not initialized : please initialize PETSc first by calling "
                          "petscInitialize(argc, argv) as the very first line of your program.");
      }
      system = new feLinearSystemPETSc(argc, argv, bilinearForms, numUnknowns);
      break;
#else
      return feErrorMsg(FE_STATUS_ERROR,
                        "feNG must be compiled with MPI and PETSc to solve with PETSc.");
#endif
    default:
      return feErrorMsg(FE_STATUS_ERROR, "Unsupported linear solver type.");
  }

  return FE_STATUS_OK;
}

feLinearSystem::feLinearSystem(std::vector<feBilinearForm *> bilinearForms)
  : recomputeMatrix(false)
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