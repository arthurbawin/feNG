#include "feLinearSystem.h"
#include "feLinearSystemMklPardiso.h"
#include "feLinearSystemPETSc.h"

extern bool petscWasInitialized;

feStatus createLinearSystem(feLinearSystem *&system,
                            linearSolverType type,
                            std::vector<feBilinearForm*> bilinearForms,
                            feMetaNumber *metaNumber,
                            feMesh *mesh,
                            int argc,
                            char **argv)
{
	if(metaNumber->getNbUnknowns() == 0)
		return feErrorMsg(FE_STATUS_ERROR, "0 unknowns : attempting to create a 0 x 0 linear system.");
	if(type == MKLPARDISO){
#if defined(HAVE_MKL)
		system = new feLinearSystemMklPardiso(bilinearForms, metaNumber, mesh);
#else
		return feErrorMsg(FE_STATUS_ERROR, "feNG must be compiled with Intel MKL to solve with Pardiso.");
#endif
	} else if(type == PETSC){
#if defined(HAVE_PETSC)
		if(argc == 0 || argv == nullptr)
			return feErrorMsg(FE_STATUS_ERROR, "Please provide argc and argv to create a PETSc linear system.");
		if(!petscWasInitialized)
			return feErrorMsg(FE_STATUS_ERROR, "PETSc was not initialized : please initialize PETSc first by calling petscInitialize(argc, argv) as the very first line of your program.");
		system = new feLinearSystemPETSc(argc, argv, bilinearForms, metaNumber, mesh);
#else
		return feErrorMsg(FE_STATUS_ERROR, "feNG must be compiled with MPI and PETSc to solve with PETSc.");
#endif
	} else{
		return feErrorMsg(FE_STATUS_ERROR, "Unsupported linear solver type.");
	}
	return FE_STATUS_OK;
}