#include "feLinearSystem.h"

extern int FE_VERBOSE;

extern bool wasInitialized;

feStatus createLinearSystem(feLinearSystem *&system,
                            const linearSolverType type,
                            const std::vector<feBilinearForm*> bilinearForms,
                            const int numUnknowns,
                            const int ownershipSplit)
{
#if !defined(HAVE_MKL)
  UNUSED(ownershipSplit);
#endif

#if !defined(HAVE_PETSC) && !defined(HAVE_MKL)
  UNUSED(system);
#endif

  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "LINEAR SYSTEM:");

  for(auto *form : bilinearForms) {
    if(form == nullptr) {
      return feErrorMsg(
        FE_STATUS_ERROR,
        "Null pointer in vector of bilinear forms, maybe you forgot to initialize it.");
    }
  }

  if(numUnknowns == 0)
    return feErrorMsg(FE_STATUS_ERROR,
                      "0 unknowns : attempting to create a linear system of size 0.");

  if(type == MKLPARDISO) {
#if defined(HAVE_MKL)
      system = new feLinearSystemMklPardiso(bilinearForms, numUnknowns, ownershipSplit);
#else
      return feErrorMsg(FE_STATUS_ERROR,
                        "feNG must be compiled with Intel MKL to solve with MKL Pardiso.");
#endif
  } else if(type == PETSC || type == PETSC_MUMPS) {
#if defined(HAVE_PETSC)
    if(!wasInitialized) {
      return feErrorMsg(FE_STATUS_ERROR,
                        "PETSc was not initialized : please initialize PETSc first by calling "
                        "initialize(argc, argv) as the very first line of your program.");
    }
    if(type == PETSC_MUMPS) {
      #if !defined(PETSC_HAVE_MUMPS)
        return feErrorMsg(FE_STATUS_ERROR, "Cannot solve with MUMPS through PETSc"
          " because PETSc was not configured with MUMPS.\n"
          " Add the option --with-mumps=1 or --download-mumps=1 when configuring PETSc.");
      #endif
    }
    system = new feLinearSystemPETSc(bilinearForms, numUnknowns, type);
#else
    return feErrorMsg(FE_STATUS_ERROR,
                        "feNG must be compiled with MPI and PETSc to solve with PETSc.");
#endif
  } else {
    return feErrorMsg(FE_STATUS_ERROR, "Unsupported choice of linear solver.");
  }
  
  return FE_STATUS_OK;
}

feLinearSystem::feLinearSystem(const std::vector<feBilinearForm*> bilinearForms)
{
  _recomputeMatrix = false;
  _numMatrixForms = 0;
  _numResidualForms = bilinearForms.size();

#if defined(HAVE_OMP)
  int nThreads = omp_get_max_threads();
#else
  int nThreads = 1;
#endif

  for(int i = 0; i < nThreads; ++i) {
    for(const feBilinearForm *f : bilinearForms) {
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
}

feLinearSystem::~feLinearSystem()
{
  // Delete the copies of the (bi)linear forms
#if defined(HAVE_OMP)
  for(auto *fCpy : _formResiduals)
    delete fCpy;
#endif
}