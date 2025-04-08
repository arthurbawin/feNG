#include "feNG.h"
#include "feMessage.h"
#if defined(HAVE_PETSC)
  #include "petscsys.h"
#endif
#if defined(HAVE_MPI)
  #include "mpi.h"
#endif
#if defined(HAVE_GMSH)
 #include "gmsh.h"
#endif

using namespace std::chrono_literals;

double tic(int mode)
{
  static std::chrono::steady_clock::time_point begin;

  if(mode == 0){
    begin = std::chrono::steady_clock::now();
    return 0.;
  }
  else {
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::milli> d = end - begin;
    double timeInSeconds  = d / 1.0s;
    return timeInSeconds;
  }
}

double toc() { return tic(1); }

bool wasInitialized = false;
#if defined(HAVE_GMSH)
bool gmshWasInitialized = false;
extern int FE_VERBOSE;
#endif

void initialize(int argc, char **argv)
{
#if defined(HAVE_PETSC)
  if(!wasInitialized) {
    wasInitialized = true;
    PetscErrorCode ierr = PetscInitialize(&argc, &argv, (char *)0, nullptr);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  } else {
    feErrorMsg(FE_STATUS_ERROR, "PETSc was already initialized\n");
    return;
  }
#elif defined(HAVE_MPI)
  if(!wasInitialized) {
    wasInitialized = true;
    MPI_Init(&argc, &argv);
  } else {
    feErrorMsg(FE_STATUS_ERROR, "MPI was already initialized\n");
    return;
  }
#else
  // Do nothing
  UNUSED(argc, argv);
#endif
}

void finalize()
{
#if defined(HAVE_PETSC)
  PetscFinalize();
#elif defined(HAVE_MPI)
  MPI_Finalize();
#else
  // Do nothing
#endif
}

#if defined(HAVE_GMSH)
void initializeGmsh()
{
  if(!gmshWasInitialized) {

    // Get max threads *before* initializing gmsh
    #if defined(HAVE_OMP)
    int maxNumThreads = omp_get_max_threads();
    #endif

    gmsh::initialize();

    #if defined(HAVE_OMP)
    gmsh::option::setNumber("General.NumThreads", maxNumThreads);
    #endif

    gmshWasInitialized = true;
  }

  if(FE_VERBOSE == VERBOSE_NONE) gmsh::option::setNumber("General.Verbosity", 2);
}

void finalizeGmsh()
{
  if(gmshWasInitialized) {
    gmsh::clear();
    gmsh::finalize();
    gmshWasInitialized = false;
  }
}
#endif