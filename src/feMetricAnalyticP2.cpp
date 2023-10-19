#include <complex>
#include "feMetric.h"
#include "feMatrixInterface.h"
#include "../contrib/root_finder/root_finder.hpp"

extern int FE_VERBOSE;

#define TOL_DISC 1e-14
#define TOL_HESS 1e-14
#define TOL_LEADING_COEFF 1e-14
#define TOL_ZERO 1e-14
#define TOL_CHECK 1e-6

double LARGE_VALUE = 1e16;

#define PRINT(name) printScalar(#name, (name))
#define PRINTVEC(name) printVec(#name, (name))
#define PRINTMAT(name) printMat(#name, (name))
void printScalar(const char *name, double val) {
  feInfo("%s = %+-1.16e", name, val);
} 
void printVec(const char *name, double vec[2]) {
  feInfo("%s = %+-1.16e - %+-1.16e", name, vec[0], vec[1]);
} 
void printMat(const char *name, const SquareMatrix &m) {
  feInfo("%s = %+-1.16e - %+-1.16e - %+-1.16e - %+-1.16e", name, m(0,0), m(0,1), m(1,0), m(1,1));
}

// To quickly sort 3 double
static void sort3(double *d){
#define SWAP(x, y)                                                             \
  if(d[y] < d[x]) {                                                            \
    double tmp = d[x];                                                         \
    d[x] = d[y];                                                               \
    d[y] = tmp;                                                                \
  }
  SWAP(0, 1);
  SWAP(1, 2);
  SWAP(0, 1);
#undef SWAP
}

double getFirstDegreeRoot(const double a, const double b)
{
  if(fabs(a) < TOL_LEADING_COEFF) return LARGE_VALUE;
  return -b/a;
}

int getSecondDegreeRoots(const double a, const double b, const double c, double roots[2])
{
  if(fabs(a) < TOL_LEADING_COEFF) {
    // Polynomial is actually linear, but we still keep a large root
    roots[0] = LARGE_VALUE;
    roots[1] = getFirstDegreeRoot(b,c);
    return 1;
  }
  double rho = b*b - 4.*a*c;
  if(rho < 0) {
    // Return real and imaginary parts of the two complex conjugate roots
    roots[0] = -b/(2*a);
    roots[1] = sqrt(-rho)/(2*a);
    return -1;
  } else {
    // Return the two real roots
    roots[0] = (-b - sqrt(rho))/(2*a);
    roots[1] = (-b + sqrt(rho))/(2*a);
    return 0;
  }
}

inline double sign(double x){return x > 0 ? 1. : (x < 0 ? -1. : 0.); };
inline double sqrt3(double x){return x >= 0 ? pow(x,1./3.) : -pow(-x, 1./3.); };
inline double maxAbs(double x, double y) { return fabs(x) > fabs(y) ? x : y; };

// Solve ax^3 + bx^2 + cx + d = 0.
// Adapted from J-M Mirebeau's routines.
int getCubicPolynomialRoots(const double a, const double b, const double c, const double d, double roots[3])
{
  feInfo("Coeff = %+-1.10e - %+-1.10e - %+-1.10e - %+-1.10e", a,b,c,d);

  if(fabs(a) < TOL_LEADING_COEFF) {
    // Cubic is actually a second degree polynomial, but we still need a large first root
    roots[0] = LARGE_VALUE;
    return getSecondDegreeRoots(b,c,d,roots+1);
  }

  double disc = b*b*c*c - 4.*a*c*c*c - 4.*b*b*b*d + 18.*a*b*c*d - 27.*a*a*d*d;
  double d0 = b*b - 3.*a*c;

  if(fabs(disc) < TOL_DISC) {
    // Multiple real roots
    if(fabs(d0) < TOL_ZERO) {
      // One triple real root
      roots[0] = roots[1] = roots[2] = -b/(3*a);
      return 3;
    } else {
      // One simple real root and one double real root
      roots[0] = (4.*a*b*c - 9.*a*a*d - b*b*b)/(a*d0);
      roots[1] = roots[2] = (9.*a*d-b*c)/(2.*d0);
      return 2;
    }
  } else {

    const double d1 = 2.*b*b*b - 9.*a*b*c + 27.*a*a*d;
    const double delta = d1*d1 - 4.*d0*d0*d0;

    PRINT(d0);
    PRINT(d1);
    PRINT(delta);
    PRINT(disc);

    if(disc > 0) {
      // Three distinct real roots (but complex parameters in the process)
      if(delta > 0) {
        feInfo("Discriminant and delta should have different sign! 1");
        feInfo("disc = %+-1.10e == delta = %+-1.10e", disc, delta);
        exit(-1);
      }

      std::complex<double> xi1(-0.5,  sqrt(3.)/2.);
      std::complex<double> xi2(-0.5, -sqrt(3.)/2.);
      std::complex<double> C = pow( std::complex<double>(d1/2., sqrt(-delta)/2.) , 1./3.);

      roots[0] = std::real(-(b + C + d0/C)/(3.*a));
      roots[1] = std::real(-(b + xi1*C + xi2*d0/C)/(3.*a));
      roots[2] = std::real(-(b + xi2*C + xi1*d0/C)/(3.*a));
      return 0;
    } else {
      // One real root and two complex conjugate roots
      if(delta < 0) {
        if(delta / disc > 1e-10) {
          feInfo("Discriminant and delta should have different sign! 2");
          feInfo("disc = %+-1.10e == delta = %+-1.10e", disc, delta);
          exit(-1);
        }
      }
      
      double C = sqrt3(maxAbs(-d1 + sqrt(fabs(delta)), -d1 - sqrt(fabs(delta)))/2.);
      // double C = sqrt3(-d1 - sqrt(fabs(delta))/2.);
      roots[0] = (-b + C + d0/C)/(3.*a);
      roots[1] = (-b - 0.5*d0/C - 0.5*C)/(3*a); // Real part of the two complex roots
      roots[2] = fabs((C - d0/C)/(2.*sqrt(3.)*a)); // Positive imaginary part
      // PRINT(roots[0]);
      // PRINT(roots[1]);
      // PRINT(roots[2]);
      // if(delta < 0) exit(-1);
      return -1;
    }
  }
  return -2;
}

// Get the type of factorization for cubic homogeneous polynomial
// based on hessian (H) and discriminant (D).
// Returns:
//  - 1 if H = 0: factors into a triple root
//  - 2 if H != 0 and D = 0: factors into a double root and a simple root
//  - 3 if H != 0 and D < 0: factors into a real root and two complex conjugate roots
//  - 4 if H != 0 and D > 0: factors into three real roots
int getFactorizationCase(const double a, const double b, const double c, const double d, double &disc)
{
  double hess1 = b*b - 3.*a*c;
  double hess2 = b*c - 9.*a*d;
  double hess3 = c*c - 3.*b*d;
  disc = b*b*c*c - 4.*a*c*c*c - 4.*b*b*b*d + 18.*a*b*c*d - 27.*a*a*d*d;
  if(fabs(hess1) < TOL_HESS && fabs(hess2) <= TOL_HESS && fabs(hess3) <= TOL_HESS)
    return 1;
  if(fabs(disc) < TOL_DISC)
    return 2;
  if(disc < 0)
    return 3;
  return 4;
}

// Computes the transformation phi and its inverse
void getPhiAndInverse(const double a, const double b, const double c, const double d, const int factorizationCase,
  const double disc, SquareMatrix &phi, SquareMatrix &invPhi)
{
  switch(factorizationCase) {
    case 1 :
      // Phi is unused for pi o phi = x^3
      return;
    case 2 :
      {
        // disc = 0: pi o phi = x^2y
        double d0 = -(2.*sqrt(pow(b*b-3.*a*c,3)) - 9.*a*b*c + 2.*b*b*b)/(27.*a*a);
        double r1, r3;
        if(fabs(d-d0) < TOL_ZERO) {
          r1 = -(      sqrt(pow(b*b-3.*a*c, 3)) - 3.*a*b*c + b*b*b)/(3.*a*b*b - 9.*a*a*c);
          r3 = -(-2. * sqrt(pow(b*b-3.*a*c, 3)) - 3.*a*b*c + b*b*b)/(3.*a*b*b - 9.*a*a*c);
        } else {
          r1 = -(-     sqrt(pow(b*b-3.*a*c, 3)) - 3.*a*b*c + b*b*b)/(3.*a*b*b - 9.*a*a*c);
          r3 = -( 2. * sqrt(pow(b*b-3.*a*c, 3)) - 3.*a*b*c + b*b*b)/(3.*a*b*b - 9.*a*a*c);
        }

        // PRINT(r1);
        // PRINT(r3);

        // double roots[3];
        // int res = getCubicPolynomialRoots(a,b,c,d,roots);
        // if(res != 2) {
        //   feInfo("There should be a double root");
        //   exit(-1);
        // }
        // PRINT(roots[0]);
        // PRINT(roots[1]);
        // PRINT(roots[2]);


        invPhi(0,0) = 1.;
        invPhi(0,1) = -r1;
        invPhi(1,0) = a;
        invPhi(1,1) = -a*r3;
        invPhi.inverse(phi);
        return;
      }
    case 3 :
      {
        // disc < 0: pi o phi = x(x^2 + 3y^2)
        // double cbrt2 = pow(2., 1./3.);
        // double v1 = pow(sqrt( pow(-27. *a*a *d + 9. *a *b *c - 2 *b*b*b, 2) + 4. * pow(3. *a *c - b*b, 3)) - 27. *a*a *d + 9. *a *b *c - 2. *b*b*b, 1./3.);
        // double v2 = cbrt2 * (3. * a *c - b*b);
        // double v3 = 3. * a * cbrt2;

        // double r1 = v1/v3 - v2/(3. *a *v1) - b/(3. *a);
        // double v5 = -2. *a *b *r1 - 3. *a*a *r1*r1 - 4. *a *c + b*b;
        // double real = -(b + a*r1)/(2.*a);
        // double imag = sqrt(fabs(v5)) / (2.*a);

        double roots[3];
        int res = getCubicPolynomialRoots(a,b,c,d,roots);
        if(res != -1) {
          feErrorMsg(FE_STATUS_ERROR, "Expected three real roots.");
          exit(-1);
        }

        double r1 = roots[0];
        double real = roots[1];
        double imag = roots[2];

        double factor = a / std::cbrt(2.*disc);
        phi(0,0) = factor * (r1 * 2. * real - 2.* (real*real + imag*imag));
        phi(0,1) = factor * (-2.*imag*r1*sqrt(3.));
        phi(1,0) = factor * (2.*r1 - 2*real);
        phi(1,1) = factor * (-2.*imag*sqrt(3.));
        phi.inverse(invPhi);
        return;
      }
    case 4 :
      {
        // disc > 0: pi o phi = x(x^2 - 3y^2)
        double roots[3];
        int res = getCubicPolynomialRoots(a,b,c,d,roots);
        if(res != 0) {
          feErrorMsg(FE_STATUS_ERROR, "Expected three real roots.");
          exit(-1);
        }

        sort3(roots);

        double r1 = roots[0];
        double r2 = roots[1];
        double r3 = roots[2];

        double factor = a / sqrt3(2.*disc);
        phi(0,0) = factor * (r1 * (r2+r3) - 2.* r2*r3);
        phi(0,1) = factor * ((r2-r3)*r1*sqrt(3.));
        phi(1,0) = factor * (2.*r1 - (r2+r3));
        phi(1,1) = factor * ((r2-r3)*sqrt(3.));
        phi.inverse(invPhi);
        return;
      }
      return;
    default:
      feErrorMsg(FE_STATUS_ERROR, "Unexpected factorizationCase. Aborting.");
      exit(-1);
  }
}

// Computes the matrix hPi associated to the largest ellipse included in the error
// curve 1 of pi(x,y) = ax^3 + bx^2y + cxy^2 + dy^3 without constraint on the diameter.
// Only valid if the discriminant of pi is nonzero since there is no unconstrained ellipse if disc = 0.
// Also assigns betaPi, related to the diameter of the ellipse.
void getUnconstrainedOptimalMetric(const double a, const double b, const double c, const double d, const int factorizationCase, const double disc, 
  const SquareMatrix &phi, const SquareMatrix &invPhi, MetricTensor &hPi, double &betaPi)
{
  if(factorizationCase == 1 || factorizationCase == 2)
    return;

  double cbrt2 = pow(2., 1./3.);
  if(factorizationCase == 3) {
    // hPi = 2^(1/3) * invPhi^T * invPhi
    hPi(0,0) = cbrt2 * (invPhi(0,0)*invPhi(0,0) + invPhi(1,0)*invPhi(1,0));
    hPi(0,1) = cbrt2 * (invPhi(0,0)*invPhi(0,1) + invPhi(1,0)*invPhi(1,1));
    hPi(1,0) = hPi(0,1);
    hPi(1,1) = cbrt2 * (invPhi(0,1)*invPhi(0,1) + invPhi(1,1)*invPhi(1,1));
  } else if(factorizationCase == 4) {
    double B = b/3.;
    double C = c/3.;
    double factor = 3. / cbrt2 * pow(disc, -1./3.);
    hPi(0,0) = factor * (2. * (B*B- a *C));
    hPi(0,1) = factor * (B*C - a*d);
    hPi(1,0) = factor * (B*C - a*d);
    hPi(1,1) = factor * (2. * (C*C - B*d));
  }
  betaPi = hPi.getMinEigenvalue();
}

// Computes res = m*x
void multiplyMat2(const SquareMatrix &m, const double x[2], double res[2]) {
  res[0] = m(0,0) * x[0] + m(0,1) * x[1];
  res[1] = m(1,0) * x[0] + m(1,1) * x[1];
}

thread_local std::set<double> ROOTS;
thread_local Eigen::VectorXd COEFF_DEG3 = Eigen::VectorXd::Zero(4);

// Computes:
// - mu, the diameter of the largest disc centered at 0 and inscribed in the level curve 1 of pi
// - alpha, the diameter of the largest ellipse inscribed in the level curve and containing the disc of radius mu
// - U, the rotation matrix associated to zPi the contact point of the largest circle inscribed in the level curve
void computeDiameters(const double a, const double b, const double c, const double d, const double disc,
  const int factorizationCase, const SquareMatrix &phi, const SquareMatrix &invPhi, SquareMatrix &U, double &mu, double &alpha)
{
  double roots[3];

  //////////////////////////
  // COEFF_DEG3(0) = -b;
  // COEFF_DEG3(1) =  3.*a-2.*c;
  // COEFF_DEG3(2) =  2.*b-3.*d;
  // COEFF_DEG3(3) =  c;
  // ROOTS = RootFinder::solvePolynomial(COEFF_DEG3, 0., INFINITY, 1e-14);

  // // if(ROOTS.size() == 0) {
  // //   feErrorMsg(FE_STATUS_ERROR, "No root found.");
  // //   exit(-1);
  // // }

  // for(auto it = ROOTS.begin(); it != ROOTS.end(); it++) {
  //   std::cout << *it << std::endl;
  // }
  //////////////////////////

  int res = getCubicPolynomialRoots(-b,3.*a-2.*c,2.*b-3.*d,c,roots);
  int nRealRoots = (res == -1) ? 1 : 3;
  double zPi[2], max_pz = 0.;
  for(int i = 0; i < nRealRoots; ++i) {
    double z = roots[i];
    PRINT(z);
    double pz = fabs( (a*z*z*z + b*z*z + c*z + d)/pow(sqrt(z*z+1.),3) );
    if(pz > max_pz) {
      zPi[0] = z;
      max_pz = pz;
    }
  }

  // if(zPi[0] >= LARGE_VALUE) {
  //   feInfo("Ici");
  //   zPi[0] = 1.;
  //   zPi[1] = 0.;
  // } else {
    double normzPi = sqrt(zPi[0]*zPi[0] + 1.);
    zPi[0] /= normzPi;
    zPi[1] = 1./normzPi;
  // }

  mu = pow(sqrt3(max_pz), 2);

  // zPi orthogonal
  double zPiO[2];
  zPiO[0] = -zPi[1];
  zPiO[1] =  zPi[0];
  U(0,0) = zPi[0];  U(0,1) = zPi[1];
  U(1,0) = zPiO[0]; U(1,1) = zPiO[1];

  if(factorizationCase == 1) {
    // Triple root, there exists phi s.t. pi o phi = x^3.
    alpha = 0.;
    return;
  }

  // invPhi*zPi, phi^T*zPi, then same for zPiO 
  double zI[2], zT[2], zIO[2], zTO[2];
  multiplyMat2(invPhi, zPi, zI);
  multiplyMat2(invPhi, zPiO, zIO);
  zT[0] = phi(0,0) * zPi[0] + phi(1,0) * zPi[1];
  zT[1] = phi(0,1) * zPi[0] + phi(1,1) * zPi[1];
  zTO[0] = phi(0,0) * zPiO[0] + phi(1,0) * zPiO[1];
  zTO[1] = phi(0,1) * zPiO[0] + phi(1,1) * zPiO[1];
  PRINTVEC(zPi);
  PRINTVEC(zPiO);
  PRINTVEC(zI);
  PRINTVEC(zIO);
  PRINTMAT(phi);
  PRINTMAT(invPhi);

  switch(factorizationCase) {
    case 2 :
      {
        // disc = 0: pi o phi = x^2y
        if(fabs(zI[0]) < TOL_ZERO || fabs(zTO[0]) < TOL_ZERO) {
          feErrorMsg(FE_STATUS_ERROR, "FIXME: Unexpected zero component after rotation. Case 2.");
          exit(-1);
        }
        double lambdaC = mu * zT[0]/zI[0];
        alpha = lambdaC * zIO[0]/zTO[0];
        double check1 = 4./(27.*lambdaC*lambdaC) * zI[1] - mu*zT[1];
        double check2 = 4./(27.*lambdaC*lambdaC) * zIO[1] - alpha*zTO[1];
        if(fabs(check1) > TOL_CHECK || fabs(check2) > TOL_CHECK) {
          feErrorMsg(FE_STATUS_ERROR, "Check failed for case 2: check1 = %+-1.4e - check2 = %+-1.4e",
            check1, check2);
          exit(-1);
        }
        return;
      }
    case 3 :
      {
        // disc < 0: pi o phi = x(x^2 + 3y^2)
        double lambdaC, upperBound = 2.;
        if(fabs(zI[1]) < TOL_ZERO) {
          lambdaC = upperBound;
        } else {
          lambdaC = mu * zT[1]/zI[1];
        }
        if(fabs(zTO[1]) < TOL_ZERO) {
          // Should not happen even if invPhi = phi^T (rotation) because zPi and zPiO are orthogonal
          feErrorMsg(FE_STATUS_ERROR, "Simultaneous zero component in rotated eigenvectors. Should not happen?");
          exit(-1);
        }
        alpha = lambdaC * zIO[1]/zTO[1];
        double check1 = (4. - sign(disc) * lambdaC*lambdaC*lambdaC)/(3*lambdaC*lambdaC) * zI[0] - mu*zT[0];
        double check2 = (4. - sign(disc) * lambdaC*lambdaC*lambdaC)/(3*lambdaC*lambdaC) * zIO[0] - alpha*zTO[0];
        if(fabs(check1) > TOL_CHECK || fabs(check2) > TOL_CHECK) {
          PRINT(zI[1]);
          PRINT(zT[1]);
          feInfo("%+-1.8e", zT[1]/zI[1]);
          PRINT(a);
          PRINT(b);
          PRINT(c);
          PRINT(d);
          PRINT(mu);
          PRINT(lambdaC);
          PRINT(alpha);
          feErrorMsg(FE_STATUS_ERROR, "Check failed for case 3: check1 = %+-1.4e - check2 = %+-1.4e",
            check1, check2);
          exit(-1);
        }
        return;
      }
    case 4 :
      // disc > 0: pi o phi = x(x^2 - 3y^2)
      {
        double V1[2] = {-8.660254037844386e-01, 0.5}; // Second row of rotation of pi/3
        double V2[2] = {-8.660254037844386e-01,-0.5}; // Second row of rotation of 2pi/3

        double zItmp_y1 = V1[0] * zI[0] + V1[1] * zI[1];
        double zTtmp_y1 = V1[0] * zT[0] + V1[1] * zT[1];
        double zItmp_y2 = V2[0] * zI[0] + V2[1] * zI[1];
        double zTtmp_y2 = V2[0] * zT[0] + V2[1] * zT[1];

        if(fabs(zI[1]) < TOL_ZERO || fabs(zItmp_y1) < TOL_ZERO || fabs(zItmp_y2) < TOL_ZERO) {
          feErrorMsg(FE_STATUS_ERROR, "FIXME: Unexpected zero component in rotated eigenvectors. Case 4.1.");
          exit(-1);
        }

        double lambdaC;
        double lambdaC0 = mu * zT[1]/zI[1];
        double lambdaC1 = mu * zTtmp_y1/zItmp_y1;
        double lambdaC2 = mu * zTtmp_y2/zItmp_y2;

        bool OK = false;
        if(lambdaC0 > 0. && lambdaC0 < 1.) { lambdaC = lambdaC0; OK = true; }
        if(lambdaC1 > 0. && lambdaC1 < 1.) {
          lambdaC = lambdaC1;
          OK = true;
          zI[0] = V1[1] * zI[0] - V1[0] * zI[1];
          zT[0] = V1[1] * zT[0] - V1[0] * zT[1];
          double zIOtmp = zIO[0], zTOtmp = zTO[0];
          zIO[0] = V1[1] * zIO[0] - V1[0] * zIO[1];
          zIO[1] = V1[0] * zIOtmp + V1[1] * zIO[1];
          zTO[0] = V1[1] * zTO[0] - V1[0] * zTO[1];
          zTO[1] = V1[0] * zTOtmp + V1[1] * zTO[1];
        }
        if(lambdaC2 > 0. && lambdaC2 < 1.) {
          lambdaC = lambdaC2;
          OK = true;
          zI[0] = V2[1] * zI[0] - V2[0] * zI[1];
          zT[0] = V2[1] * zT[0] - V2[0] * zT[1];
          double zIOtmp = zIO[0], zTOtmp = zTO[0];
          zIO[0] = V2[1] * zIO[0] - V2[0] * zIO[1];
          zIO[1] = V2[0] * zIOtmp + V2[1] * zIO[1];
          zTO[0] = V2[1] * zTO[0] - V2[0] * zTO[1];
          zTO[1] = V2[0] * zTOtmp + V2[1] * zTO[1];
        }

        if(!OK) {
          PRINT(lambdaC0);
          PRINT(lambdaC1);
          PRINT(lambdaC2);
          feErrorMsg(FE_STATUS_ERROR, "No rotation provide an acceptable lambdaC.");
          exit(-1);
        }

        if(fabs(zTO[1]) < TOL_ZERO) {
          feErrorMsg(FE_STATUS_ERROR, "FIXME: Unexpected zero component in rotated eigenvectors. Case 4.2.");
          exit(-1);
        }

        alpha = lambdaC * zIO[1]/zTO[1];
        double check1 = (4. - sign(disc) * lambdaC*lambdaC*lambdaC)/(3*lambdaC*lambdaC) * zI[0] - mu*zT[0];
        double check2 = (4. - sign(disc) * lambdaC*lambdaC*lambdaC)/(3*lambdaC*lambdaC) * zIO[0] - alpha*zTO[0];
        if(fabs(check1) > TOL_CHECK || fabs(check2) > TOL_CHECK) {
          feErrorMsg(FE_STATUS_ERROR, "Check failed for case 4: check1 = %+-1.4e - check2 = %+-1.4e",
            check1, check2);
          exit(-1);
        }

        return;
      }
    default:
      feErrorMsg(FE_STATUS_ERROR, "Unexpected factorizationCase. Aborting.");
      exit(-1);
  }
}

// Computes res = A^T * M * A with M = diag(v1,v2)
void multiplyATMAdiag(const SquareMatrix &A, const double v1, const double v2, MetricTensor &res) {
  res(0,0) = A(0,0)*A(0,0)*v1 + A(1,0)*A(1,0)*v2;
  res(0,1) = A(0,1)*A(0,0)*v1 + A(1,1)*A(1,0)*v2;
  res(1,0) = A(0,0)*A(0,1)*v1 + A(1,0)*A(1,1)*v2;
  res(1,1) = A(0,1)*A(0,1)*v1 + A(1,1)*A(1,1)*v2;
}

// Computes res = A^T * M * A
void multiplyATMA(const SquareMatrix &A, const SquareMatrix &M, MetricTensor &res) {
  res(0,0) = A(0,0)*(A(0,0)*M(0,0) + A(1,0)*M(1,0)) + A(1,0)*(A(0,0)*M(0,1) + A(1,0)*M(1,1));
  res(0,1) = A(0,1)*(A(0,0)*M(0,0) + A(1,0)*M(1,0)) + A(1,1)*(A(0,0)*M(0,1) + A(1,0)*M(1,1));
  res(1,0) = A(0,0)*(A(0,1)*M(0,0) + A(1,1)*M(1,0)) + A(1,0)*(A(0,1)*M(0,1) + A(1,1)*M(1,1));
  res(1,1) = A(0,1)*(A(0,1)*M(0,0) + A(1,1)*M(1,0)) + A(1,1)*(A(0,1)*M(0,1) + A(1,1)*M(1,1));
}

// Roots of the polynomial of degree 3 or 4 in lambdaC in the transition metric
// Obtained by imposing det(h(lambdaC) - alpha*I) = 0.
thread_local std::set<double> LAMBDA_ROOTS;
thread_local Eigen::VectorXd LAMBDA_COEFF_DEG3 = Eigen::VectorXd::Zero(4);
thread_local Eigen::VectorXd LAMBDA_COEFF_DEG4 = Eigen::VectorXd::Zero(5);

// Compute the transition metric between the largest ellipse containing Dpi (disc of largest radius included
// in the level curve 1) and the unconstrained optimal metric.
void getTransitionMetric(const int factorizationCase, const double disc, const double alpha, 
  const SquareMatrix &phi, SquareMatrix &invPhi, MetricTensor &Q)
{
  switch(factorizationCase) {
    case 1 :
      // Does not apply
      return;
    case 2 :
      {
        // Compute lambda (L) such that alpha is an eigenvalue of h = invPhi^T * H * invPhi,
        // with H = diag(L, f(L)). Solve det(h(L) - alpha*I) = 0 by applying the identity
        // det(A + B) = det(A) + det(B) + det(A) * tr(A^-1 * B) with A = -alpha*I and B = h.
        // Thus we solve alpha^2 + det(invPhi)^2 * L * f(L) - alpha * tr(h).
        // In case 2, disc = 0 thus f(L) = 4/(27*L^2).
        double C1 = invPhi(0,0)*invPhi(0,0) + invPhi(0,1)*invPhi(0,1);
        double C2 = invPhi(1,0)*invPhi(1,0) + invPhi(1,1)*invPhi(1,1);
        double det = invPhi.determinant();
        double D = det*det;
        LAMBDA_COEFF_DEG3(0) = -27.*alpha*C1;
        LAMBDA_COEFF_DEG3(1) =  27.*alpha*alpha;
        LAMBDA_COEFF_DEG3(2) =  4.*D;
        LAMBDA_COEFF_DEG3(3) = -4.*alpha*C2;
        LAMBDA_ROOTS = RootFinder::solvePolynomial(LAMBDA_COEFF_DEG3, 0., INFINITY, 1e-8);

        if(LAMBDA_ROOTS.size() == 0) {
          feErrorMsg(FE_STATUS_ERROR, "No root found for transition metric in case 2.");
          exit(-1);
        }

        // Keep the root that minimizes the determinant.
        // The transformed metric is invPhi^T * H * invPhi.
        // We have det(H) = L*f(L) and det(h) = det(invPhi)^2 * det(H),
        // so its sufficient to check det(H) without assembling the metric.
        double L, minDet = 1e22;
        for(auto it = LAMBDA_ROOTS.begin(); it != LAMBDA_ROOTS.end(); it++) {
          double l = *it;
          double detH = l * 4./(27.*l*l);
          if(detH < minDet) {
            minDet = detH;
            L = l;
          }
        }

        multiplyATMAdiag(invPhi, L, 4./(27.*L*L), Q);
        return;
      }
    case 3 :
      {
        // Here H = diag(f(L), L) and f(L) = (4 + L^3)/(3 * L^2).
        double C1 = invPhi(0,0)*invPhi(0,0) + invPhi(0,1)*invPhi(0,1);
        double C2 = invPhi(1,0)*invPhi(1,0) + invPhi(1,1)*invPhi(1,1);
        double det = invPhi.determinant();
        double D = det*det;
        LAMBDA_COEFF_DEG4(0) = D;
        LAMBDA_COEFF_DEG4(1) = -alpha*(C1 + 3.*C2);
        LAMBDA_COEFF_DEG4(2) = 3.*alpha*alpha;
        LAMBDA_COEFF_DEG4(3) = 4.*D;
        LAMBDA_COEFF_DEG4(4) = -4.*alpha*C1;
        LAMBDA_ROOTS = RootFinder::solvePolynomial(LAMBDA_COEFF_DEG4, 0., 2., 1e-8);

        if(LAMBDA_ROOTS.size() == 0) {
          feErrorMsg(FE_STATUS_ERROR, "No root found for transition metric in case 3.");
          exit(-1);
        }

        double L, minDet = 1e22;
        for(auto it = LAMBDA_ROOTS.begin(); it != LAMBDA_ROOTS.end(); it++) {
          double l = *it;
          double detH = l * (4. + l*l*l)/(3.*l*l);
          if(detH < minDet) {
            minDet = detH;
            L = l;
          }
        }

        multiplyATMAdiag(invPhi, (4. + L*L*L)/(3.*L*L), L, Q);
        return;
      }
    case 4 :
      {
        // Same as case 3, but f(L) = (4 - L^3)/(3 * L^2)
        // and we need to check the three rotations 0, pi/3, 2pi/3.
        double L, minDet = 1e22;
        size_t nAcceptableRoots = 0;
        int whichTheta;
        double allV[6] = {1., 0., 0.5, 8.660254037844386e-01, -0.5, 8.660254037844386e-01};
        for(int i = 0; i < 3; ++i) {
          double V[2] = {allV[2*i], allV[2*i+1]};
          // V*invPhi
          double m1 = invPhi(0,0)*V[0] + invPhi(1,0)*V[1];
          double m2 = invPhi(0,1)*V[0] + invPhi(1,1)*V[1];
          double m3 = invPhi(1,0)*V[0] - invPhi(0,0)*V[1];
          double m4 = invPhi(1,1)*V[0] - invPhi(0,1)*V[1];
          double C1 = m1*m1 + m2*m2;
          double C2 = m3*m3 + m4*m4;
          double det = m1*m4 - m2*m3;
          double D = det*det; // D is invariant by rotation though (-:
          LAMBDA_COEFF_DEG4(0) = -D;
          LAMBDA_COEFF_DEG4(1) = alpha*(C1 - 3.*C2);
          LAMBDA_COEFF_DEG4(2) = 3.*alpha*alpha;
          LAMBDA_COEFF_DEG4(3) = 4.*D;
          LAMBDA_COEFF_DEG4(4) = -4.*alpha*C1;
          LAMBDA_ROOTS = RootFinder::solvePolynomial(LAMBDA_COEFF_DEG4, 0., 1., 1e-8);
          nAcceptableRoots += LAMBDA_ROOTS.size();
          
          for(auto it = LAMBDA_ROOTS.begin(); it != LAMBDA_ROOTS.end(); it++) {
            double l = *it;
            double detH = l * (4. - l*l*l)/(3.*l*l);
            if(detH < minDet) {
              minDet = detH;
              L = l;
              whichTheta = i;
            }
          }
        }

        if(nAcceptableRoots == 0) {
          feErrorMsg(FE_STATUS_ERROR, "No root found for transition metric in case 4.");
          exit(-1);
        }

        if(whichTheta != 0) {
          // Apply rotation to invPhi
          double V[2] = {allV[2*whichTheta], allV[2*whichTheta+1]};
          double m1 = invPhi(0,0)*V[0] + invPhi(1,0)*V[1];
          double m2 = invPhi(0,1)*V[0] + invPhi(1,1)*V[1];
          double m3 = invPhi(1,0)*V[0] - invPhi(0,0)*V[1];
          double m4 = invPhi(1,1)*V[0] - invPhi(0,1)*V[1];
          invPhi(0,0) = m1; invPhi(0,1) = m2;
          invPhi(1,0) = m3; invPhi(1,1) = m4;
        }

        multiplyATMAdiag(invPhi, (4. - L*L*L)/(3.*L*L), L, Q);
        return;
      }
  }

}
                             
void feMetric::computeAnalyticMetricP2(std::vector<double> &errorCoeff, MetricTensor &Q, const double maxDiameter)
{
  double a = errorCoeff[0];
  double b = errorCoeff[1];
  double c = errorCoeff[2];
  double d = errorCoeff[3];

  // a = 1.;
  // b = 3.;
  // c = -1.;
  // d = -(2.*sqrt(pow(b*b-3.*a*c,3)) - 9.*a*b*c + 2.*b*b*b)/(27.*a*a);

  // a = -10.;
  // b = 5.;
  // c = 2.;
  // d = -1.;

  // a = 6.;
  // b = d = 0.;
  // c = 2.;

  // a = 1.;
  // b = 3.;
  // c = -1.;
  // d = -(2.*sqrt(pow(b*b-3.*a*c,3)) - 9.*a*b*c + 2.*b*b*b)/(27.*a*a);

  // a = 1.;
  // b = 1.;
  // c = 0.;
  // d = 0.;

  // feInfo("a = %f", a);
  // feInfo("b = %f", b);
  // feInfo("c = %f", c);
  // feInfo("d = %f", d);
  PRINT(a);
  PRINT(b);
  PRINT(c);
  PRINT(d);

  double disc;
  int cas = getFactorizationCase(a,b,c,d,disc);

  feInfo("cas = %d", cas);
  feInfo("disc = %+-1.16e", disc);

  SquareMatrix phi(2), invPhi(2), U(2), diag(2);
  feInfo("Computing phi");
  getPhiAndInverse(a,b,c,d, cas, disc, phi, invPhi);
  feInfo("Done phi");

  double muPi, alphaPi, betaPi = 0.;
  feInfo("Computing diameters");
  computeDiameters(a, b, c, d, disc, cas, phi, invPhi, U, muPi, alphaPi);
  feInfo("Done diameters");

  MetricTensor hPi;
  getUnconstrainedOptimalMetric(a,b,c,d,cas,disc,phi,invPhi,hPi,betaPi);

  double alpha = 4. / (maxDiameter*maxDiameter);

  // alpha = (muPi+alphaPi)/2.;
  // alpha = (alphaPi+betaPi)/2.;

  if(alpha >= muPi)
  {
    // Q = Dpi, the disc of radius alpha^(-1/2)
    Q(0,0) = alpha;
    Q(0,1) = 0.;
    Q(1,0) = 0.;
    Q(1,1) = alpha;
  } else if(alpha >= alphaPi)
  {
    // Ellipse of diameter 2*alpha^(-1/2) inscribed in the level curve and containing Dpi
    // Ellipse is tangent at 2 points (4 points if alpha = alphaPi)
    multiplyATMAdiag(U, muPi, alpha, Q);

  } else if(alpha >= betaPi)
  {
    // Ellipse is tangent at 4 points and has diameter 2*alpha^(-1/2).
    // Transition between the previous case and the optimal unconstrained case
    // Cannot happen for case 1 since the level curve is the one of x^3 after some transformation
    getTransitionMetric(cas, disc, alpha, phi, invPhi, Q);

  } else {
    // Optimal unconstrained case
    if(cas == 1 || cas == 2) {
      feErrorMsg(FE_STATUS_ERROR, "Unexpected case %d: bounds are %+-1.4e - %+-1.4e - %+-1.4e", cas, betaPi, alphaPi, muPi);
        exit(-1);
    }
    Q(0,0) = hPi(0,0);
    Q(0,1) = hPi(0,1);
    Q(1,0) = hPi(1,0);
    Q(1,1) = hPi(1,1);
  }

  PRINT(betaPi);
  PRINT(alphaPi);
  PRINT(muPi);

  if(Q.determinant() <= 0) {
    Q.print();
    feInfo("det = %+-1.16e", Q.determinant());
    exit(-1);
  }
}