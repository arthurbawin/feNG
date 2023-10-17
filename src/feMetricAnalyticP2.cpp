#include "feMetric.h"
#include "feMatrixInterface.h"
#include <complex>

extern int FE_VERBOSE;

#define TOL_DISC 1e-14
#define TOL_HESS 1e-14
#define TOL_LEADING_COEFF 1e-14
#define TOL_ZERO 1e-14

double LARGE_VALUE = 1e10;

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
  // feInfo("Coeff = %f - %f - %f - %f", a,b,c,d);

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
      return 0;
    } else {
      // One simple real root and one double real root
      roots[0] = (4.*a*b*c - 9.*a*a*d - b*b*b)/(a*d0);
      roots[1] = roots[2] = (9.*a*d-b*c)/(2.*d0);
      return 0;
    }
  } else {

    const double d1 = 2.*b*b*b - 9.*a*b*c + 27.*a*a*d;
    const double delta = d1*d1 - 4.*d0*d0*d0;

    if(disc > 0) {
      // Three distinct real roots (but complex parameters in the process)
      if(delta > 0) {
        feInfo("Discriminant and delta should have different sign! 1");
        feInfo("disc = %+-f == delta = %+-f", disc, delta);
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
        feInfo("Discriminant and delta should have different sign! 2");
        feInfo("disc = %+-f == delta = %+-f", disc, delta);
        exit(-1);
      }
      
      double C = sqrt3(maxAbs(-d1 + sqrt(delta), -d1 - sqrt(delta))/2.);
      feInfo("C = %f",C);
      roots[0] = (-b + C + d0/C)/(3.*a);
      roots[1] = (-b -0.5*d0/C - 0.5*C)/(3*a); // Real part of the two complex roots
      roots[2] = (C - d0/C)/(2.*sqrt(3.)*a); // Imaginary part
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
  if(hess1 < TOL_HESS && hess2 <= TOL_HESS && hess3 <= TOL_HESS)
    return 1;
  disc = b*b*c*c - 4.*a*c*c*c - 4.*b*b*b*d + 18.*a*b*c*d - 27.*a*a*d*d;
  if(disc < TOL_DISC)
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
          feInfo("VAL = %f",pow(b*b-3.*a*c, 3));
          r1 = -(      sqrt(pow(b*b-3.*a*c, 3)) - 3.*a*b*c + b*b*b)/(3.*a*b*b - 9.*a*a*c);
          r3 = -(-2. * sqrt(pow(b*b-3.*a*c, 3)) - 3.*a*b*c + b*b*b)/(3.*a*b*b - 9.*a*a*c);
        } else {
          r1 = -(-     sqrt(pow(b*b-3.*a*c, 3)) - 3.*a*b*c + b*b*b)/(3.*a*b*b - 9.*a*a*c);
          r3 = -( 2. * sqrt(pow(b*b-3.*a*c, 3)) - 3.*a*b*c + b*b*b)/(3.*a*b*b - 9.*a*a*c);
        }
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
        double cbrt2 = pow(2., 1./3.);
        double v1 = pow(sqrt( pow(-27. *a*a *d + 9. *a *b *c - 2 *b*b*b, 2) + 4. * pow(3. *a *c - b*b, 3)) - 27. *a*a *d + 9. *a *b *c - 2. *b*b*b, 1./3.);
        double v2 = cbrt2 * (3. * a *c - b*b);
        double v3 = 3. * a * cbrt2;

        double r1 = v1/v3 - v2/(3. *a *v1) - b/(3. *a);
        double v5 = -2. *a *b *r1 - 3. *a*a *r1*r1 - 4. *a *c + b*b;
        double real = -(b + a*r1)/(2.*a);
        double imag = sqrt(fabs(v5)) / (2.*a);

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
        double cbrt2 = pow(2., 1./3.);
        double v1 = pow(sqrt( pow(-27. *a*a *d + 9. *a *b *c - 2 *b*b*b, 2) + 4. * pow(3. *a *c - b*b, 3)) - 27. *a*a *d + 9. *a *b *c - 2. *b*b*b, 1./3.);
        double v2 = cbrt2 * (3. * a *c - b*b);
        double v3 = 3. * a * cbrt2;

        double r1 = v1/v3 - v2/(3. *a *v1) - b/(3. *a);
        double v4 = -2. *a *b *r1 - 3. *a*a *r1*r1 - 4. *a *c + b*b;
        if(v4 < 0) {
          feErrorMsg(FE_STATUS_ERROR,"Unexpected negative value: %+-1.6e", v4);
          exit(-1);
        }
        double r2 = -(b + a*r1 + sqrt(v4))/(2.*a);
        double r3 = -(b + a*r1 - sqrt(v4))/(2.*a);

        double factor = a / std::cbrt(2.*disc);
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
// The value of disc is not checked.
void getUnconstrainedOptimalMetric(const double a, const double b, const double c, const double d, const int factorizationCase, const double disc, 
  const SquareMatrix &phi, const SquareMatrix &invPhi, MetricTensor &hPi)
{
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
}

#define PRINT(name) printScalar(#name, (name))
#define PRINTVEC(name) printVec(#name, (name))
#define PRINTMAT(name) printMat(#name, (name))
void printScalar(const char *name, double val) {
  feInfo("%s = %+-1.4f", name, val);
} 
void printVec(const char *name, double vec[2]) {
  feInfo("%s = %+-1.4f - %+-1.4f", name, vec[0], vec[1]);
} 
void printMat(const char *name, const SquareMatrix &m) {
  feInfo("%s = %+-1.4f - %+-1.4f - %+-1.4f - %+-1.4f", name, m(0,0), m(0,1), m(1,0), m(1,1));
}

// Computes res = m*x
void multiplyMat2(const SquareMatrix &m, const double x[2], double res[2]) {
  res[0] = m(0,0) * x[0] + m(0,1) * x[1];
  res[1] = m(1,0) * x[0] + m(1,1) * x[1];
}

// Computes:
// - mu, the diameter of the largest disc centered at 0 and inscribed in the level curve 1 of pi
// - alpha, the diameter of the largest ellipse inscribed in the level curve and containing the disc of radius mu
void computeDiameters(const double a, const double b, const double c, const double d, const double disc,
  const int factorizationCase, const SquareMatrix &phi, const SquareMatrix &invPhi, double &mu, double &alpha)
{
  double roots[3];
  int res = getCubicPolynomialRoots(-b,3.*a-2.*c,2.*b-3.*d,c,roots);
  int nRealRoots = (res == -1) ? 1 : 3;
  double zPi[2], max_pz = 0.;
  for(int i = 0; i < nRealRoots; ++i) {
    double z = roots[i];
    double pz = fabs( (a*z*z*z + b*z*z + c*z + d)/pow(sqrt(z*z+1.),3) );
    if(pz > max_pz) {
      zPi[0] = z;
      max_pz = pz;
    }
  }

  double normzPi = sqrt(zPi[0]*zPi[0] + 1.);
  zPi[0] /= normzPi;
  zPi[1] = 1./normzPi;
  mu = pow(sqrt3(max_pz), 2);

  // zPi orthogonal, invPhi*zPi, phi^T*zPi, then same for zPiO 
  double zPiO[2], zI[2], zT[2], zIO[2], zTO[2];
  if(factorizationCase != 1) {
    zPiO[0] = -zPi[1];
    zPiO[1] =  zPi[0];
    multiplyMat2(invPhi, zPi, zI);
    multiplyMat2(invPhi, zPiO, zIO);
    zT[0] = phi(0,0) * zPi[0] + phi(1,0) * zPi[1];
    zT[1] = phi(0,1) * zPi[0] + phi(1,1) * zPi[1];
    zTO[0] = phi(0,0) * zPiO[0] + phi(1,0) * zPiO[1];
    zTO[1] = phi(0,1) * zPiO[0] + phi(1,1) * zPiO[1];
    if(fabs(zI[0]) < TOL_ZERO || fabs(zTO[0]) < TOL_ZERO || fabs(zI[1]) < TOL_ZERO || fabs(zTO[1]) < TOL_ZERO) {
      feErrorMsg(FE_STATUS_ERROR, "FIXME: zero component in rotated eigenvectors.");
      exit(-1);
    }
  }

  switch(factorizationCase) {
    case 1 :
      // Triple root, there exists phi s.t. pi o phi = x^3.
      alpha = 0.;
      break;
    case 2 :
      {
        // disc = 0: pi o phi = x^2y
        double lambdaC = mu * zT[0]/zI[0];
        alpha = lambdaC * zIO[0]/zTO[0];
        double check1 = 4./(27.*lambdaC*lambdaC) * zI[1] - mu*zT[1];
        double check2 = 4./(27.*lambdaC*lambdaC) * zIO[1] - alpha*zTO[1];
        if(fabs(check1) > TOL_ZERO || fabs(check2) > TOL_ZERO) {
          feErrorMsg(FE_STATUS_ERROR, "Check failed for case 2: check1 = %+-1.4e - check2 = %+-1.4e",
            check1, check2);
          exit(-1);
        }
        return;
      }
    case 3 :
      {
        // disc < 0: pi o phi = x(x^2 + 3y^2)
        double upperBound = 2.;
        double lambdaC = mu * zT[1]/zI[1];
        alpha = lambdaC * zIO[1]/zTO[1];
        double check1 = (4. - sign(disc) * lambdaC*lambdaC*lambdaC)/(3*lambdaC*lambdaC) * zI[0] - mu*zT[0];
        double check2 = (4. - sign(disc) * lambdaC*lambdaC*lambdaC)/(3*lambdaC*lambdaC) * zIO[0] - alpha*zTO[0];
        if(fabs(check1) > TOL_ZERO || fabs(check2) > TOL_ZERO) {
          feErrorMsg(FE_STATUS_ERROR, "Check failed for case 3: check1 = %+-1.4e - check2 = %+-1.4e",
            check1, check2);
          exit(-1);
        }
        return;
      }
    case 4 :
      // disc > 0: pi o phi = x(x^2 - 3y^2)
      break;
    default:
      feErrorMsg(FE_STATUS_ERROR, "Unexpected factorizationCase. Aborting.");
      exit(-1);
  }
}
                             
void feMetric::computeAnalyticMetricP2(std::vector<double> &errorCoeff, MetricTensor &Q, const double maxDiameter)
{
  double a = errorCoeff[0];
  double b = errorCoeff[1];
  double c = errorCoeff[2];
  double d = errorCoeff[3];

  a = 1.;
  b = 3.;
  c = -1.;
  d = -(2.*sqrt(pow(b*b-3.*a*c,3)) - 9.*a*b*c + 2.*b*b*b)/(27.*a*a);
  feInfo("a = %f", a);
  feInfo("b = %f", b);
  feInfo("c = %f", c);
  feInfo("d = %f", d);

  double disc;
  int cas = getFactorizationCase(a,b,c,d,disc);

  feInfo("cas = %d", cas);
  feInfo("disc = %f", disc);

  SquareMatrix phi(2), invPhi(2);
  getPhiAndInverse(a,b,c,d, cas, disc, phi, invPhi);

  // MetricTensor hPi;
  // getUnconstrainedOptimalMetric(a,b,c,d,disc,hPi,phi,invPhi);

  double mu, alpha;
  computeDiameters(a, b, c, d, disc, cas, phi, invPhi, mu, alpha);

  PRINT(mu);
  PRINT(alpha);

  exit(-1);
}