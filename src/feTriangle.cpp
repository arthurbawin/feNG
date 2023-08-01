#include "feTriangle.h"
#include <iostream>

bool TriangleP1::isInside(double xyz[3], double tol) 
{
  printf("FIXME isInside");
  exit(-1);
  return true;
}

bool TriangleP1::xyz2uvw(double xyz[3], double uvw[3])
{
  const double O[3] = {_v[0]->x(), _v[0]->y(), _v[0]->z()};

  const double d[3] = {xyz[0] - O[0], xyz[1] - O[1], xyz[2] - O[2]};
  const double d1[3] = {_v[1]->x() - O[0], _v[1]->y() - O[1], _v[1]->z() - O[2]};
  const double d2[3] = {_v[2]->x() - O[0], _v[2]->y() - O[1], _v[2]->z() - O[2]};

  const double Jxy = d1[0] * d2[1] - d1[1] * d2[0];
  const double Jxz = d1[0] * d2[2] - d1[2] * d2[0];
  const double Jyz = d1[1] * d2[2] - d1[2] * d2[1];

  if((fabs(Jxy) > fabs(Jxz)) && (fabs(Jxy) > fabs(Jyz))) {
    uvw[0] = (d[0] * d2[1] - d[1] * d2[0]) / Jxy;
    uvw[1] = (d[1] * d1[0] - d[0] * d1[1]) / Jxy;
  } else if(fabs(Jxz) > fabs(Jyz)) {
    uvw[0] = (d[0] * d2[2] - d[2] * d2[0]) / Jxz;
    uvw[1] = (d[2] * d1[0] - d[0] * d1[2]) / Jxz;
  } else {
    uvw[0] = (d[1] * d2[2] - d[2] * d2[1]) / Jyz;
    uvw[1] = (d[2] * d1[1] - d[1] * d1[2]) / Jyz;
  }
  uvw[2] = 0.0;
  return true;
}

bool TriangleP2::isInside(double xyz[3], double tol) 
{
  bool isInside = true;

  for(int i = 0; i < 3; ++i)
  {
    double x0 = _v[i]->x();
    double x1 = _v[(i+1) % 3]->x();
    double xm = _v[i+3]->x();
    double y0 = _v[i]->y();
    double y1 = _v[(i+1) % 3]->y();
    double ym = _v[i+3]->y();
    // Get coefficients of parabola y = ax^2 + bx + c
    double det = (x0 - x1)*(x0*x1 - x0*xm - x1*xm + xm*xm);
    double a = -(x0*y1 - x1*y0 - x0*ym + xm*y0 + x1*ym - xm*y1) / det;
    double b = (x0*x0*y1 - x1*x1*y0 - x0*x0*ym + xm*xm*y0 + x1*x1*ym - xm*xm*y1) / det;
    double c = -(- ym*x0*x0*x1 + y1*x0*x0*xm + ym*x0*x1*x1 - y1*x0*xm*xm - y0*x1*x1*xm + y0*x1*xm*xm) / det;
    bool orient = xyz[1] - (a * xyz[0] * xyz[0] + b * xyz[0] + c) - tol > 0;
    isInside &= orient;
  }

  return isInside;
}

// This is defined in multiple places, this should be fixed
void phiTriP2(const double r, const double s, double phi[6], double dphidr[6], double dphids[6]){
  phi[0] = (1. - r - s) * (1. - 2.*r - 2.*s);
  phi[1] = r * (2. * r - 1.);
  phi[2] = s * (2. * s - 1.);
  phi[3] = 4. * r * (1. - r - s);
  phi[4] = 4. * r * s;
  phi[5] = 4. * s * (1. - r - s);

  dphidr[0] = 4. * (r+s) - 3.;
  dphidr[1] = 4. * r - 1.;
  dphidr[2] = 0.;
  dphidr[3] = 4.*(1.-2.*r-s);
  dphidr[4] = 4.*s;
  dphidr[5] = -4.*s;

  dphids[0] = 4.*(r+s)-3.;
  dphids[1] = 0.;
  dphids[2] = 4.*s-1.;
  dphids[3] = -4.*r;
  dphids[4] = 4.*r;
  dphids[5] = 4.*(1.-r-2.*s);
}

bool TriangleP2::xyz2uvw(double xyz[3], double uvw[3])
{
  // printf("Looking for point %f - %f\n", xyz[0], xyz[1]);
  // printf("in element (%f - %f) - (%f - %f) - (%f - %f)\n", xyz[0], xyz[1]);
  double phi[6], dphidr[6], dphids[6];

  // Newton-Raphson to find (r,s) s.t. F(r,s) = (x,y)
  double J[4], invJ[4], f[2], r[2] = {0., 0.};
  int iter = 0, maxiter = 20;
  bool success = false, stop = false;
  double tol = 1e-5;

  // Use linear solution as initial Newton guess
  const double O[2] = {_v[0]->x(), _v[0]->y()};
  const double d[2] = {xyz[0] - O[0], xyz[1] - O[1]};
  const double d1[2] = {_v[1]->x() - O[0], _v[1]->y() - O[1]};
  const double d2[2] = {_v[2]->x() - O[0], _v[2]->y() - O[1]};
  const double Jxy = d1[0] * d2[1] - d1[1] * d2[0];
  r[0] = (d[0] * d2[1] - d[1] * d2[0]) / Jxy;
  r[1] = (d[1] * d1[0] - d[0] * d1[1]) / Jxy;

  while(!stop) {

    phiTriP2(r[0], r[1], phi, dphidr, dphids);

    // Tangent matrix and residual
    f[0] = f[1] = 0.;
    J[0] = J[1] = J[2] = J[3] = 0.;
    for(int i = 0; i < 6; ++i) {
      f[0] += _v[i]->x() * phi[i];
      f[1] += _v[i]->y() * phi[i];
      J[0] += _v[i]->x() * dphidr[i];
      J[1] += _v[i]->x() * dphids[i];
      J[2] += _v[i]->y() * dphidr[i];
      J[3] += _v[i]->y() * dphids[i];
    }
    f[0] -= xyz[0];
    f[1] -= xyz[1];

    // Inverse
    double det = J[0]*J[3] - J[1]*J[2];
    invJ[0] =  J[3]/det;
    invJ[1] = -J[1]/det;
    invJ[2] = -J[2]/det;
    invJ[3] =  J[0]/det;

    // Solve
    r[0] -= invJ[0] * f[0] + invJ[1] * f[1];
    r[1] -= invJ[2] * f[0] + invJ[3] * f[1];

    // printf("NR iter %2d: %f - %f ---- %f - %f\n", iter, r[0], r[1], f[0], f[1]);

    if(sqrt(f[0]*f[0] + f[1]*f[1]) < tol || iter > maxiter) {
      stop = true;
      if(sqrt(f[0]*f[0] + f[1]*f[1]) < tol)
        success = true;
    }
    iter++;
  }

  uvw[0] = r[0];
  uvw[1] = r[1];
  uvw[2] = 0.0;

  return success;
}