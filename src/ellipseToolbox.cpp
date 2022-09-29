#include "ellipseToolbox.h"

#include <cstdio>
#include <iostream>
#include <cmath>
#include <initializer_list>

// Get discriminant of the conic
double discriminant(std::vector<double> &p) {
  double a = p[0];
  double b = p[1] / 2.0;
  double c = p[2];
  return a * c - b * b;
}

// Get radius for an implicit circle
double imCircle(std::vector<double> &p) {
  double a = p[0];
  double d = p[3];
  double e = p[4];
  double f = p[5];
  double val = (d * d + e * e) / (4.0 * a * a) - f / a;
  if(val < 0) { printf("In imCircle : Warning : circle radius is imaginary.\n"); }
  return sqrt(val);
}

// Get semi-axes for an implicit ellipse
void imEllipse(std::vector<double> &p, double *semiA, double *semiB) {
  double a = p[0];
  double b = p[1] / 2.0;
  double c = p[2];
  double d = p[3] / 2.0;
  double f = p[4] / 2.0;
  double g = p[5];
  double q =
    2.0 * (a * f * f + c * d * d + g * b * b - 2.0 * b * d * f - a * c * g) / (b * b - a * c);
  double r = (a - c) * (a - c) + 4.0 * b * b;
  if(r < 0) { printf("In imEllipse : Warning : r is negative.\n"); }

  r = sqrt(r);

  if(q / (r - (a + c)) < 0) { printf("In imEllipse : Warning : q/(r-(a+c)) is negative.\n"); }
  if(q / (-r - (a + c)) < 0) { printf("In imEllipse : Warning : q/(-r-(a+c)) is negative.\n"); }
  double semi_a = sqrt(q / (r - (a + c))); // major axis
  double semi_b = sqrt(q / (-r - (a + c))); // minor axis

  *semiA = semi_a;
  *semiB = semi_b;
}

// Get translation vector from the origin
void imConicTranslation(std::vector<double> &p, std::vector<double> &t) {
  double a = p[0];
  double b = p[1] / 2.0;
  double c = p[2];
  double d = p[3] / 2.0;
  double f = p[4] / 2.0;
  double D = discriminant(p);
  if(D > 1e-12) {
    t[0] = -(c * d - b * f) / D;
    t[1] = -(a * f - d * b) / D;
  }
}

// Translates a conic section given as an implicit equation
void imConicTranslate(std::vector<double> &p, std::vector<double> &t) {
  double a = p[0];
  double b = p[1];
  double c = p[2];
  double d = p[3];
  double f = p[4];
  double g = p[5];

  double x = t[0];
  double y = t[1];

  double dt = d - 2.0 * a * x - b * y;
  double ft = f - b * x - 2.0 * c * y;
  double gt = g + a * x * x + b * x * y + c * y * y - d * x - f * y;

  p[0] = a;
  p[1] = b;
  p[2] = c;
  p[3] = dt;
  p[4] = ft;
  p[5] = gt;
}

double acot(double x) {
  if(fabs(x) < 1e-14) { printf("In acot : Warning : x = %10.16e - 1/x = %10.16e\n", x, 1. / x); }
  return atan(1. / x);
}

// Get rotation matrix
void imConicRotation(std::vector<double> &p, std::vector<double> &R) {
  double a = p[0];
  double b = p[1];
  double c = p[2];

  double D = discriminant(p);
  double ra = a;
  double rc = c;

  double phi, sin_phi, cos_phi;

  if(fabs(b) < 1e-12) { // b == 0 : conic is aligned with axes
    if(ra > rc) {
      phi = M_PI / 2.0;
      sin_phi = 1.0;
      cos_phi = 0.0;
    } else {
      phi = 0.0;
      sin_phi = 0.0;
      cos_phi = 1.0;
    }
  } else { // conic is rotated
    if(fabs(a - c) > 1e-12) { // a != c
      phi = acot((a - c) / b) / 2.0;
      if(ra > rc) phi += M_PI / 2.0;
    } else {
      phi = -M_PI / 4.0 * sgn(b);
    }
    sin_phi = sin(phi);
    cos_phi = cos(phi);
  }

  // R is row major
  // R = [ cos_phi -sin_phi ; sin_phi cos_phi ]
  R[0] = cos_phi;
  R[1] = -sin_phi;
  R[2] = sin_phi;
  R[3] = cos_phi;
}

bool getExplicitEllipse(std::vector<double> &p, double *semiA, double *semiB,
                        std::vector<double> &R, std::vector<double> &t) {
  double a = p[0];
  double b = p[1] / 2.0;
  double c = p[2];
  double d = p[3] / 2.0;
  double f = p[4] / 2.0;
  double g = p[5];

  // std::cout << a << " - " << b << " - " << c << " - " << d << " - " << f << " - " << g <<
  // std::endl;

  if(fabs(d) > 1e-10 || fabs(f) > 1e-10)
    printf("In getExplicitEllipse : Warning : ignoring extra implicit coefficients beyond p[2] (d "
           "and/or f).\n");

  // Check that the coefficients are those of an ellipse
  double D = discriminant(p);
  double delta = -g * b * b + 2 * b * d * f - c * d * d - a * f * f + a * c * g;

  // std::cout << D << " - " << delta << std::endl;

  if(D > 1e-12) {
    if(fabs(delta) > 1e-12 && delta / (a + c) < 0) {
      if(fabs(a - c) < 1e-12 && fabs(b) < 1e-12) { // a == c and b == 0 : Circle
        *semiA = imCircle(p);
        *semiB = *semiA;
      } else { // Ellipse
        imEllipse(p, semiA, semiB);
      }
    } else {
      printf("In getExplicitEllipse : Warning : coefficients describe a degenerate ellipse.\n");
      return false;
    }
  } else {
    printf("In getExplicitEllipse : Warning : coefficients describe a hyperbola.\n");
    return false;
  }

  imConicTranslation(p, t);
  t[0] *= -1.0;
  t[1] *= -1.0;
  imConicTranslate(p, t);
  imConicRotation(p, R);
  return true;
}

/* Return the discretization (x,y) of the ellipse given by the implicit form ax² + b*x*y + cy² = 1
   The ellipse is centered in (xC,yC).
   This ellipse is the unit circle associated to the metric tensor :
    ( a b )
    ( b c )
*/

bool getEllipsePoints(double a, double b, double c, double xC, double yC, std::vector<double> &x,
                      std::vector<double> &y) {
  std::vector<double> p = {a, b, c, 0., 0., -1.};
  std::vector<double> t(2, 0.);
  std::vector<double> R(4, 0.);
  double semiA, semiB;

  bool res = getExplicitEllipse(p, &semiA, &semiB, R, t);

  if(res){
    size_t size = x.size();

    double xTmp, yTmp;

    for(size_t i = 0; i < size; ++i) {
      xTmp = semiA * cos(i * 2.0 * M_PI / size);
      yTmp = semiB * sin(i * 2.0 * M_PI / size);

      x[i] = xC + R[0] * xTmp + R[1] * yTmp;
      y[i] = yC + R[2] * xTmp + R[3] * yTmp;
    }
    return true;
  } else{
    return false;
  }
}
