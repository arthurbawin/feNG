#ifndef _FENG_
#define _FENG_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <wchar.h>
#include <float.h>
#include <time.h>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <complex>
#include <assert.h>
#include <sstream>
#include <map>
#include <unordered_map>
#include <set>
#include <iomanip>
#include <limits>
#include <tuple>
#include <utility>
#include <chrono>
#include <thread>

#if defined(HAVE_OMP)
#include "omp.h"
#endif

#if defined(HAVE_MPI)
  // This is advised here: https://github.com/open-mpi/ompi/issues/5157
  // to silence cast-function-type warnings when compiling with -Wextra
  #define OMPI_SKIP_MPICXX 1
#endif

#if defined(HAVE_PETSC)
#include "petscsystypes.h"
#endif

// Macro to specify unused parameters
// and avoid compiler warnings
// #define UNUSED(...) (void)sizeof(__VA_ARGS__)
template<class... T> void UNUSED(T&&...){}

#if defined(HAVE_PETSC)
typedef PetscInt feInt;
#elif defined(HAVE_MKL)
typedef long int feInt;
#else
typedef long int feInt;
#endif

#if defined(HAVE_MKL)
#include "mkl_types.h"

// typedef long long int PardisoInt;
// typedef int PardisoInt;
typedef MKL_INT PardisoInt;

// Define format specifier for printf
#ifdef MKL_ILP64
/* oneMKL ILP64 integer types */
#define MKL_INT_FORMAT "%lld"
#else
/* oneMKL LP64 integer types */
#define MKL_INT_FORMAT "%d"
#endif
#endif

double tic(int mode = 0);
double toc();

//
// Initialize/finalize routines, to call at the beginning/end of a main function.
// If compiled with PETSc, they call PetscInitialize() and PetscFinalize().
// If not, but if compiled with MPI, they call MPI_Init() and MPI_Finalize() instead.
// If compiled with neither PETSc nor MPI, these are dummy functions.
//
void initialize(int argc, char **argv);
void finalize();

#if defined(HAVE_GMSH)
void initializeGmsh();
void openGmshModel(std::string &filename);
void finalizeGmsh();
#endif

#endif
