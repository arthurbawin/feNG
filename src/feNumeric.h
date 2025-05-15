
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

template <typename T>
bool solveNewtonRaphson1D(const T &f,
                          const T &dfdx,
                          const int maxiter,
                          const double tol,
                          double &x)
{
  bool success = false, stop = false;
  int iter = 0;
  double fx = f(x);

  while(!stop)
  {
    x -= fx / dfdx(x);

    fx = f(x);
    iter++;
    feInfo("NR - Iter %d: x = %+-1.6e - fx = %+-1.6e", iter, x, fx);

    if(fabs(fx) <= tol || iter > maxiter) {
      stop = true;
      if(fabs(fx) <= tol) success = true;
    }
  }

  return success;
}

// Compute matrix exponential of 2x2 symmetric matrix [a b; b c]
template <class MetricType>
void exponentialMatrix(const double a, const double b, const double c, MetricType &expm)
{
  // Check if a = c
  if(fabs(a - c) < 1e-14) {
    double expamb = exp(a - b), exp2b = exp(2. * b);
    expm(0, 0) = (expamb * (exp2b + 1.0)) / 2.0;
    expm(0, 1) = (expamb * (exp2b - 1.0)) / 2.0;
    expm(1, 0) = expm(0, 1);
    expm(1, 1) = (expamb * (exp2b + 1.0)) / 2.0;
  } else {
    double S = sqrt((a - c) * (a - c) + 4. * b * b);
    double oneOverS = 1. / S;
    double expapcs2 = exp(0.5 * (a + c));
    double sinhss2 = sinh(0.5 * S);
    double coshss2 = cosh(0.5 * S);
    expm(0, 0) =
      coshss2 * expapcs2 + a * sinhss2 * expapcs2 * oneOverS - c * sinhss2 * expapcs2 * oneOverS;
    expm(0, 1) = b * sinhss2 * expapcs2 * oneOverS * 2.0;
    expm(1, 0) = expm(0, 1);
    expm(1, 1) =
      coshss2 * expapcs2 - a * sinhss2 * expapcs2 * oneOverS + c * sinhss2 * expapcs2 * oneOverS;
  }
}

template <class MetricType>
void classicalP1Interpolation(const double *xsi, const MetricType &M0, const MetricType &M1,
                              const MetricType &M2, MetricType &result)
{
  // Hardcoded P1 Lagrange basis functions for simplicity
  double phi[3] = {1. - xsi[0] - xsi[1], xsi[0], xsi[1]};
  result = M0 * phi[0] + M1 * phi[1] + M2 * phi[2];
}

template <class MetricType>
void logEuclidianP1Interpolation(const double *xsi, const MetricType &logM0,
                                 const MetricType &logM1, const MetricType &logM2,
                                 MetricType &result)
{
  // Hardcoded P1 Lagrange basis functions for simplicity
  double phi[3] = {1. - xsi[0] - xsi[1], xsi[0], xsi[1]};
  // Interpolate log(M) then take exponential
  result = logM0 * phi[0] + logM1 * phi[1] + logM2 * phi[2];
  result = result.exp();
}

// Hardcoded symbolic expression of the matrix exponential
template <class MetricType>
void logEuclidianP1InterpolationExplicit(const double *xsi, const MetricType &logM0,
                                         const MetricType &logM1, const MetricType &logM2,
                                         MetricType &result)
{
  double phi[3] = {1. - xsi[0] - xsi[1], xsi[0], xsi[1]};
  double L11[3] = {logM0(0, 0), logM1(0, 0), logM2(0, 0)};
  double L12[3] = {logM0(0, 1), logM1(0, 1), logM2(0, 1)};
  double L22[3] = {logM0(1, 1), logM1(1, 1), logM2(1, 1)};

  double a = 0., b = 0., c = 0.;
  for(size_t i = 0; i < 3; ++i) {
    a += phi[i] * L11[i];
    b += phi[i] * L12[i];
    c += phi[i] * L22[i];
  }

  // Hardcoded exponential matrix
  exponentialMatrix(a, b, c, result);

  if(isnan(result(0, 0)) || isnan(result(0, 1)) || isnan(result(1, 1))) {
    feInfo("metric is nan :%+-1.10e - %+-1.10e - %+-1.10e", result(0, 0), result(0, 1),
           result(1, 1));
    exit(-1);
  }
}

template <class MetricType>
void logEuclidianP2Interpolation(const double *xsi, const MetricType &logM0,
                                 const MetricType &logM1, const MetricType &logM2,
                                 const MetricType &logM3, const MetricType &logM4,
                                 const MetricType &logM5, MetricType &result)
{
  // Hardcoded P2 Lagrange basis functions for simplicity
  double phi[6] = {(1. - xsi[0] - xsi[1]) * (1. - 2. * xsi[0] - 2. * xsi[1]),
                   xsi[0] * (2. * xsi[0] - 1.),
                   xsi[1] * (2. * xsi[1] - 1.),
                   4. * xsi[0] * (1. - xsi[0] - xsi[1]),
                   4. * xsi[0] * xsi[1],
                   4. * xsi[1] * (1. - xsi[0] - xsi[1])};
  // Interpolate log(M) then take exponential
  result = (logM0 * phi[0] + logM1 * phi[1] + logM2 * phi[2] + logM3 * phi[3] + logM4 * phi[4] +
            logM5 * phi[5])
             .exp();
}

// Hardcoded symbolic expression of the matrix exponential
template <class MetricType>
void logEuclidianP2InterpolationExplicit(const double *xsi, const MetricType &logM0,
                                         const MetricType &logM1, const MetricType &logM2,
                                         const MetricType &logM3, const MetricType &logM4,
                                         const MetricType &logM5, MetricType &result)
{
  double phi[6] = {(1. - xsi[0] - xsi[1]) * (1. - 2. * xsi[0] - 2. * xsi[1]),
                   xsi[0] * (2. * xsi[0] - 1.),
                   xsi[1] * (2. * xsi[1] - 1.),
                   4. * xsi[0] * (1. - xsi[0] - xsi[1]),
                   4. * xsi[0] * xsi[1],
                   4. * xsi[1] * (1. - xsi[0] - xsi[1])};
  double L11[6] = {logM0(0, 0), logM1(0, 0), logM2(0, 0), logM3(0, 0), logM4(0, 0), logM5(0, 0)};
  double L12[6] = {logM0(0, 1), logM1(0, 1), logM2(0, 1), logM3(0, 1), logM4(0, 1), logM5(0, 1)};
  double L22[6] = {logM0(1, 1), logM1(1, 1), logM2(1, 1), logM3(1, 1), logM4(1, 1), logM5(1, 1)};

  double a = 0., b = 0., c = 0.;
  for(size_t i = 0; i < 6; ++i) {
    a += phi[i] * L11[i];
    b += phi[i] * L12[i];
    c += phi[i] * L22[i];
  }

  // Hardcoded exponential matrix
  exponentialMatrix(a, b, c, result);

  if(isnan(result(0, 0)) || isnan(result(0, 1)) || isnan(result(1, 1))) {
    feInfo("metric is nan :%+-1.10e - %+-1.10e - %+-1.10e", result(0, 0), result(0, 1),
           result(1, 1));
    exit(-1);
  }
}

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
