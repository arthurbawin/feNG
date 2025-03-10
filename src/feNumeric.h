
#ifndef _FENUMERIC_
#define _FENUMERIC_

//
// Adapted from Gmsh's Numeric.h
// Toolbox of various operations on e.g. vectors and matrices
//

inline double dotProdN(int n, double *a, double *b)
{
	double res = 0.;
  for(int i = 0; i < n; ++i) {
    res += a[i]*b[i];
  }
  return res;
}
double det2x2(double mat[2][2]);
double det3x3(double mat[3][3]);
double inv3x3(double mat[3][3], double inv[3][3]);
double inv2x2(double mat[2][2], double inv[2][2]);

//
// Compute the convergence rate for the given error values
// assuming isotropic meshes, i.e., involving log((nElm_i/nElm_i-1)^(1/dim))
//
void computeAndPrintConvergence(const int dim,
                                const int nMesh,
                                const std::vector<double> &error,
                                const std::vector<int> &nElm,
                                const int nSignificantDigits = 8,
                                std::ostream &os = std::cout);

#endif
