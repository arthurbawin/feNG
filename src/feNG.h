#ifndef _EF5_
#define _EF5_

#define _TRACE_

using namespace std;

// LES .h DU SYSTEME
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
// #include <mpi.h>
//#include "mkl.h"

const int maxchar=1024;

// LES MACROS LES PLUS UTILES
#ifndef signe
#define signe(a)   ( ((a)>=0.0) ? (1.0) : (-1.0) )
#endif

// NAMING CONVENTION POUR INTERFACE AVEC FORTRAN SOUS UNIX
#define ef5_abs  labs

// LES VARIABLES LOGIQUES
#define EFint           long int 
#define EFPARDISOint    int
#define EFMKLPARDISOint long long int
// LES VARIABLES LOGIQUES
#define VRAI       (1)
#define FAUX       (0)

#endif

// LES ENUMEROTATIONS

#ifdef _OPENMP_
extern "C" {
#include "omp.h"
}
#endif
