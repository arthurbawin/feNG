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
#include <iomanip>
#include <limits>
#include <tuple>
#include <utility>
#include <chrono>
#include <thread>

#include "omp.h"

#if defined(HAVE_PETSC)
#include "petscsystypes.h"
#endif

#if defined(HAVE_PETSC)
typedef PetscInt feInt;
#elif defined(HAVE_MKL)
typedef long int feInt;
#else
typedef long int feInt;
#endif

#if defined(HAVE_MKL)
typedef int fePardisoInt;
typedef long long int feMKLPardisoInt;
#endif

void tic(int mode = 0);
void toc();

#endif
