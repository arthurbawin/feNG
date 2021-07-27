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

#if defined(HAVE_PETSC)
#include "petscsystypes.h"
#endif

#if defined(HAVE_PETSC)
typedef PetscInt64 feInt;
#elif defined(HAVE_MKL)
typedef long int feInt;
typedef int fePardisoInt;
typedef long long int feMKLPardisoInt;
#else
typedef long int feInt;
#endif

// void log(){}

// template<typename First, typename ...Rest>
// void log(First && first, Rest && ...rest)
// {
//     std::cout << std::forward<First>(first);
//     log(std::forward<Rest>(rest)...);
// }

#endif
