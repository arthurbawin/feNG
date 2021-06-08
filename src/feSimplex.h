#ifndef _FESIMPLEX_
#define _FESIMPLEX_

#include <list>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <cmath>
#include <assert.h>

#define MAX_N 10001
#define MAX_M 10001

typedef long long lld;
typedef unsigned long long llu;


void pivot(int x, int y);
int iterate_simplex();
int initialise_simplex();
std::pair<std::vector<double>, double> simplex(int nInput, int mInput, std::vector<double> AInput, std::vector<double> bInput, std::vector<double> cInput, double vInput);


#endif