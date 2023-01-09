#include "feLinearSystem.h"
#include "feLinearSystemMklPardiso.h"
#include "feLinearSystemPETSc.h"

extern int FE_VERBOSE;

extern bool petscWasInitialized;

feStatus createLinearSystem(feLinearSystem *&system, linearSolverType type,
                            std::vector<feSpace *> allFESpaces,
                            std::vector<feBilinearForm *> bilinearForms, feMetaNumber *metaNumber,
                            feMesh *mesh, int argc, char **argv)
{
  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "LINEAR SYSTEM:");

  if(metaNumber->getNbUnknowns() == 0)
    return feErrorMsg(FE_STATUS_ERROR, "0 unknowns : attempting to create a linear system of size 0.");

  // Check that all (bi-)linear forms are defined on existing connectivities.
  // The finite element spaces of a single form must be defined on the same geometric connectivity.
  for(feBilinearForm *form : bilinearForms) {
    bool OK = false;
    for(feSpace *space : allFESpaces) {
      if(space->getCncGeoTag() == form->getCncGeoTag()) OK = true;
    }
    if(!OK) {
      return feErrorMsg(
        FE_STATUS_ERROR,
        "(Bi-)linear form %s is not defined on any known geometric connectivity (named subdomain).",
        form->getIDName().c_str());
    }
  }

  switch(type) {
    case MKLPARDISO:
#if defined(HAVE_MKL)
      system = new feLinearSystemMklPardiso(bilinearForms, metaNumber, mesh);
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
      system = new feLinearSystemPETSc(argc, argv, bilinearForms, metaNumber, mesh);
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