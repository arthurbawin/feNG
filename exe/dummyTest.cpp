#include "SVector3.h"
#include "STensor3.h"
#include "ellipseToolbox.h"
#include <cmath>
#include <iostream>
#include "fullMatrix.h"

SVector3 E_X(1., 0., 0.);
SVector3 E_Y(0., 1., 0.);
SVector3 E_Z(0., 0., 1.);

static inline SVector3 interpolateDirTri2(double *dir, uint32_t n0, uint32_t n1, uint32_t n2,
                                          double *uvw) {
  double t0 = fmin(atan2(-dir[4 * n0 + 1], -dir[4 * n0]),
                   atan2(dir[4 * n0 + 1], dir[4 * n0])); // Angle associated to e1 at vertex 0
  double t1 = fmin(atan2(-dir[4 * n1 + 1], -dir[4 * n1]),
                   atan2(dir[4 * n1 + 1], dir[4 * n1])); // Angle associated to e1 at vertex 1
  double t2 = fmin(atan2(-dir[4 * n2 + 1], -dir[4 * n2]),
                   atan2(dir[4 * n2 + 1], dir[4 * n2])); // Angle associated to e1 at vertex 2

  // printf("Angles : %f deg et %f deg\n", atan2(-dir[4*n0+1],-dir[4*n0])/M_PI*180.0,
  // atan2(dir[4*n0+1],dir[4*n0])/M_PI*180.0);
  printf("Angle au noeud 0 = %f rad = %f deg\n", t0, t0 / M_PI * 180.0);

  // printf("Angles : %f deg et %f deg\n", atan2(-dir[4*n1+1],-dir[4*n1])/M_PI*180.0,
  // atan2(dir[4*n1+1],dir[4*n1])/M_PI*180.0);
  printf("Angle au noeud 0 = %f rad = %f deg\n", t1, t1 / M_PI * 180.0);

  // printf("Angles : %f deg et %f deg\n", atan2(-dir[4*n2+1],-dir[4*n2])/M_PI*180.0,
  // atan2(dir[4*n2+1],dir[4*n2])/M_PI*180.0);
  printf("Angle au noeud 0 = %f rad = %f deg\n", t2, t2 / M_PI * 180.0);

  double cI = (1. - uvw[0] - uvw[1]) * cos(t0) + uvw[0] * cos(t1) + uvw[1] * cos(t2);
  double sI = (1. - uvw[0] - uvw[1]) * sin(t0) + uvw[0] * sin(t1) + uvw[1] * sin(t2);
  double tI = atan2(sI, cI);

  printf("Angle interpolé = %f rad = %f deg\n", tI, tI / M_PI * 180.0);

  return SVector3(cos(tI), sin(tI), 0.);
}

bool isInside(double u, double v, double w) {
  double tol = 1e-12;
  if(u < (-tol) || v < (-tol) || u > ((1. + tol) - v) || fabs(w) > tol) return false;
  return true;
}

static inline SVector3 interpolateDirTri(double *dir, uint32_t n0, uint32_t n1, uint32_t n2,
                                         double *uvw) {
  // double t0 = (dir[9*n0]*dir[9*n0+1] < 0) ? atan2(-dir[9*n0+1],dir[9*n0]) :
  // atan2(dir[9*n0+1],dir[9*n0]); // Angle associated to e0 at vertex 0 double t1 =
  // (dir[9*n1]*dir[9*n1+1] < 0) ? atan2(-dir[9*n1+1],dir[9*n1]) : atan2(dir[9*n1+1],dir[9*n1]); //
  // Angle associated to e0 at vertex 1 double t2 = (dir[9*n2]*dir[9*n2+1] < 0) ?
  // atan2(-dir[9*n2+1],dir[9*n2]) : atan2(dir[9*n2+1],dir[9*n2]); // Angle associated to e0 at
  // vertex 2
  double t0 = fmin(atan2(-dir[4 * n0 + 1], -dir[4 * n0]),
                   atan2(dir[4 * n0 + 1], dir[4 * n0])); // Angle associated to e0 at vertex 0
  double t1 = fmin(atan2(-dir[4 * n1 + 1], -dir[4 * n1]),
                   atan2(dir[4 * n1 + 1], dir[4 * n1])); // Angle associated to e0 at vertex 1
  double t2 = fmin(atan2(-dir[4 * n2 + 1], -dir[4 * n2]),
                   atan2(dir[4 * n2 + 1], dir[4 * n2])); // Angle associated to e0 at vertex 2
  // double t0 = atan2(dir[4*n0+1],dir[4*n0]); // Angle associated to e0 at vertex 0
  // double t1 = atan2(dir[4*n1+1],dir[4*n1]); // Angle associated to e0 at vertex 1
  // double t2 = atan2(dir[4*n2+1],dir[4*n2]); // Angle associated to e0 at vertex 2

  // // Angle de e0 est dans [0,pi]
  // if(t0 < 0) t0 += M_PI/2.
  // if(t1 < 0) t1 += M_PI/2.
  // if(t2 < 0) t2 += M_PI/2.

  // Pour eviter les retournements : configuration avec seulement des variations entre 0 et pi/2
  // double d01 = fabs(t0-t1);
  // double d12 = fabs(t1-t2);
  // double d20 = fabs(t2-t0);

  // if(d01 > M_PI/2. && d12 > M_PI/2.) t1 += M_PI;
  // if(d01 > M_PI/2. && d20 > M_PI/2.) t0 += M_PI;
  // if(d12 > M_PI/2. && d20 > M_PI/2.) t2 += M_PI;

  double a = 2.0;
  double c = (1. - uvw[0] - uvw[1]) * cos(a * t0) + uvw[0] * cos(a * t1) + uvw[1] * cos(a * t2);
  double s = (1. - uvw[0] - uvw[1]) * sin(a * t0) + uvw[0] * sin(a * t1) + uvw[1] * sin(a * t2);
  double tI = 1.0 / a * atan2(s, c);

  return SVector3(cos(tI), sin(tI), 0.);
}

int mainTri(int argc, char **argv) {
  double theta0 = 45.0;
  double theta1 = 60.0;
  double theta2 = -45.0;

  theta0 *= M_PI / 180.0;
  // theta0 = atan2(-0.580021, 0.814601);
  double thetaOrth0 = theta0 + M_PI / 2.;

  theta1 *= M_PI / 180.0;
  // theta1 = atan2(-0.028528, 0.999593);
  double thetaOrth1 = theta1 + M_PI / 2.;

  theta2 *= M_PI / 180.0;
  // theta2 = atan2(-0.128028, -0.991771);
  double thetaOrth2 = theta2 + M_PI / 2.;

  double *dir = (double *)malloc(sizeof(double) * 4 * 3);

  dir[4 * 0 + 0] = cos(theta0);
  dir[4 * 0 + 1] = sin(theta0);
  dir[4 * 0 + 2] = cos(thetaOrth0);
  dir[4 * 0 + 3] = sin(thetaOrth0);

  dir[4 * 1 + 0] = cos(theta1);
  dir[4 * 1 + 1] = sin(theta1);
  dir[4 * 1 + 2] = cos(thetaOrth1);
  dir[4 * 1 + 3] = sin(thetaOrth1);

  dir[4 * 2 + 0] = cos(theta2);
  dir[4 * 2 + 1] = sin(theta2);
  dir[4 * 2 + 2] = cos(thetaOrth2);
  dir[4 * 2 + 3] = sin(thetaOrth2);

  for(int i = 0; i < 12; ++i) std::cout << dir[i] << std::endl;

  double uvw[3] = {1.0 / 3.0, 1.0 / 3.0, 0.};

  SVector3 res = interpolateDirTri(dir, 0, 1, 2, uvw);

  res.print("Vec interpolé");

  int np = 30;
  double *xp = (double *)malloc(sizeof(double) * np * np);
  double *yp = (double *)malloc(sizeof(double) * np * np);

  FILE *myfile = fopen("triTest4.pos", "w");
  fprintf(myfile, "View \"triTest\"{\n");

  fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", 0., 0., 0., 1., 0.,
          0., 1., 1.);
  fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", 1., 0., 0., 0., 1.,
          0., 1., 1.);
  fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", 0., 1., 0., 0., 0.,
          0., 1., 1.);

  for(int i = 0; i < np; ++i) {
    for(int j = 0; j < np; ++j) {
      xp[np * i + j] = i * 1.0 / np;
      yp[np * i + j] = j * 1.0 / np;

      uvw[0] = xp[np * i + j];
      uvw[1] = yp[np * i + j];
      uvw[2] = 0.0;

      if(isInside(uvw[0], uvw[1], uvw[2])) {
        res = interpolateDirTri(dir, 0, 1, 2, uvw);
        SVector3 v = crossprod(res, SVector3(0., 0., 1.));

        // std::cout << "norme "<< res.norm() << std::endl;

        fprintf(myfile, "VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n", xp[np * i + j],
                yp[np * i + j], 0., 0.25 * res[0], 0.25 * res[1], 0.);
        fprintf(myfile, "VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n", xp[np * i + j],
                yp[np * i + j], 0., 0.5 * v[0], 0.5 * v[1], 0.);
      }
    }
  }

  fprintf(myfile, "VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n", 0., 0., 0., 0.8 * dir[4 * 0 + 0],
          0.8 * dir[4 * 0 + 1], 0.);
  fprintf(myfile, "VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n", 0., 0., 0., dir[4 * 0 + 2],
          dir[4 * 0 + 3], 0.);
  fprintf(myfile, "VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n", 1., 0., 0., 0.8 * dir[4 * 1 + 0],
          0.8 * dir[4 * 1 + 1], 0.);
  fprintf(myfile, "VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n", 1., 0., 0., dir[4 * 1 + 2],
          dir[4 * 1 + 3], 0.);
  fprintf(myfile, "VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n", 0., 1., 0., 0.8 * dir[4 * 2 + 0],
          0.8 * dir[4 * 2 + 1], 0.);
  fprintf(myfile, "VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n", 0., 1., 0., dir[4 * 2 + 2],
          dir[4 * 2 + 3], 0.);

  // for(int i = 0; i < np*np; ++i){
  // fprintf(myfile,"VP(%.16g,%.16g,%.16g){%.16g, %.16g, %.16g};\n",xp[i],yp[i],0., 0.5, 0.5, 0.);
  // fprintf(myfile,"VP(%.16g,%.16g,%.16g){%.16g, %.16g, %.16g};\n",xp[i],yp[i],0., 0.5, 0.5, 0.);
  // }
  fprintf(myfile, "};");
  fclose(myfile);

  free(dir);
  free(xp);
  free(yp);

  return 0;
}

SMetric3 interpolateMetricSquare(SMetric3 m1, SMetric3 m2, SMetric3 m3, SMetric3 m4, double xsi,
                                 double eta) {
  std::cout << " Interpolating at (" << xsi << " - " << eta << ")" << std::endl;
  SMetric3 im1 = m1.invert(); // m1 : (-1,-1)
  SMetric3 im2 = m2.invert(); // m2 : ( 1,-1)
  SMetric3 im3 = m3.invert(); // m3 : ( 1, 1)
  SMetric3 im4 = m4.invert(); // m4 : (-1, 1)

  im1 *= (1 - xsi) * (1 - eta) / 4.;
  im2 *= (1 + xsi) * (1 - eta) / 4.;
  im3 *= (1 + xsi) * (1 + eta) / 4.;
  im4 *= (1 - xsi) * (1 + eta) / 4.;

  im1 += im2;
  im1 += im3;
  im1 += im4;
  return im1.invert();
}

SMetric3 interpolateMetricTriangle(SMetric3 m1, SMetric3 m2, SMetric3 m3, double xsi, double eta) {
  SMetric3 im1 = m1.invert(); // m1 : ( 0, 0)
  SMetric3 im2 = m2.invert(); // m2 : ( 0, 1)
  SMetric3 im3 = m3.invert(); // m3 : ( 1, 0)

  // if(isnan(im1(0,0)) || isnan(im1(0,1)) || isnan(im1(1,1)) || isnan(im1(0,2)) || isnan(im1(1,2))
  // || isnan(im1(2,2))){
  //   Msg::Warning("im1 is nan");
  // }
  // if(isnan(im2(0,0)) || isnan(im2(0,1)) || isnan(im2(1,1)) || isnan(im2(0,2)) || isnan(im2(1,2))
  // || isnan(im2(2,2))){
  //   Msg::Warning("im2 is nan");
  // }
  // if(isnan(im3(0,0)) || isnan(im3(0,1)) || isnan(im3(1,1)) || isnan(im3(0,2)) || isnan(im3(1,2))
  // || isnan(im3(2,2))){
  //   Msg::Warning("im3 is nan");
  // }

  im1 *= 1.0 - xsi - eta;
  im2 *= xsi;
  im3 *= eta;

  im1 += im2;
  im1 += im3;
  return im1.invert();
}

int mainOctreeInterpolation(int argc, char **argv) {
  int nt = 30;
  std::vector<double> x(nt, 0.);
  std::vector<double> y(nt, 0.);

  int nPoints = 8;
  std::vector<double> xP(nPoints * nPoints, 0.);
  std::vector<double> yP(nPoints * nPoints, 0.);

  for(int i = 0; i < nPoints; ++i) {
    for(int j = 0; j < nPoints; ++j) {
      xP[nPoints * i + j] = -1.0 + i * 2.0 / nPoints;
      yP[nPoints * i + j] = -1.0 + j * 2.0 / nPoints;
    }
  }

  double color = 10.0;

  FILE *myfile = fopen("octreeCellInterpolation.pos", "w");
  fprintf(myfile, "View \"octreeCellInterpolation\"{\n");

  fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", -1., -1., 0., 1., -1.,
          0., 1., 1.);
  fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", 1., -1., 0., 1., 1.,
          0., 1., 1.);
  fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", 1., 1., 0., -1., 1.,
          0., 1., 1.);
  fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", -1., 1., 0., -1., -1.,
          0., 1., 1.);

  double factor = 10.0;

  SVector3 v0(1., 0., 0.);
  v0.normalize();
  SVector3 n(0., 0., 1.);
  SVector3 v1 = crossprod(v0, n);
  v1.normalize();
  double h0 = 1.0 / factor;
  double h1 = 2.0 / factor;
  SMetric3 m1(1. / h0 / h0, 1. / h1 / h1, 1.0, v0, v1, n);
  SMetric3 m2 = m1;

  SVector3 v2(1., 1., 0.);
  v2.normalize();
  SVector3 v3 = crossprod(v0, n);
  v3.normalize();
  double h2 = 1.0 / factor;
  double h3 = 2.0 / factor;
  SMetric3 m3(1. / h2 / h2, 1. / h3 / h3, 1.0, v2, v3, n);
  SMetric3 m4 = m3;

  // Plot metrics at corners
  getEllipsePoints(m1(0, 0), 2.0 * m1(0, 1), m1(1, 1), -1., -1., x, y);
  for(int i = 0; i < nt; ++i) {
    if(i != nt - 1) {
      fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
              x[i + 1], y[i + 1], 0., color, color);
    } else {
      fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
              x[0], y[0], 0., color, color);
    }
  }
  getEllipsePoints(m2(0, 0), 2.0 * m2(0, 1), m2(1, 1), 1., -1., x, y);
  for(int i = 0; i < nt; ++i) {
    if(i != nt - 1) {
      fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
              x[i + 1], y[i + 1], 0., color, color);
    } else {
      fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
              x[0], y[0], 0., color, color);
    }
  }
  getEllipsePoints(m3(0, 0), 2.0 * m3(0, 1), m3(1, 1), 1., 1., x, y);
  for(int i = 0; i < nt; ++i) {
    if(i != nt - 1) {
      fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
              x[i + 1], y[i + 1], 0., color, color);
    } else {
      fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
              x[0], y[0], 0., color, color);
    }
  }
  getEllipsePoints(m4(0, 0), 2.0 * m4(0, 1), m4(1, 1), -1., 1., x, y);
  for(int i = 0; i < nt; ++i) {
    if(i != nt - 1) {
      fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
              x[i + 1], y[i + 1], 0., color, color);
    } else {
      fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
              x[0], y[0], 0., color, color);
    }
  }

  SMetric3 mI;

  for(int i = 0; i < nPoints; ++i) {
    for(int j = 0; j < nPoints; ++j) {
      fprintf(myfile, "SP(%.16g,%.16g,%.16g){%.16g};\n", xP[nPoints * i + j], yP[nPoints * i + j],
              0., (i + j) * color);

      // Interpolate ellipse
      mI = interpolateMetricSquare(m1, m2, m3, m4, xP[nPoints * i + j], yP[nPoints * i + j]);

      getEllipsePoints(mI(0, 0), 2.0 * mI(0, 1), mI(1, 1), xP[nPoints * i + j], yP[nPoints * i + j],
                       x, y);

      // Plot the ellipse
      for(int ii = 0; ii < nt; ++ii) {
        if(ii != nt - 1) {
          fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[ii], y[ii],
                  0., x[ii + 1], y[ii + 1], 0., color, color);
        } else {
          fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[ii], y[ii],
                  0., x[0], y[0], 0., color, color);
        }
      }
    }
  }

  fprintf(myfile, "};");
  fclose(myfile);

  return 0;
}

int mainTriangleInterpolation(int argc, char **argv) {
  int nt = 30;
  std::vector<double> x(nt, 0.);
  std::vector<double> y(nt, 0.);

  int nPoints = 14;
  std::vector<double> xP(nPoints * nPoints, 0.);
  std::vector<double> yP(nPoints * nPoints, 0.);

  for(int i = 0; i < nPoints; ++i) {
    for(int j = 0; j < nPoints; ++j) {
      xP[nPoints * i + j] = i * 2.0 / nPoints;
      yP[nPoints * i + j] = j * 2.0 / nPoints;
    }
  }

  double color = 1.0;

  FILE *myfile = fopen("triangleCellInterpolation.pos", "w");
  fprintf(myfile, "View \"triangleCellInterpolation\"{\n");

  fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", 0., 0., 0., 1., 0.,
          0., 1., 1.);
  fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", 1., 0., 0., 0., 1.,
          0., 1., 1.);
  fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", 0., 1., 0., 0., 0.,
          0., 1., 1.);

  double factor = 30.0;

  SVector3 v0(1., 2., 0.);
  v0.normalize();
  SVector3 n(0., 0., 1.);
  SVector3 v1 = crossprod(v0, n);
  v1.normalize();
  double h0 = 1.0 / factor;
  double h1 = 2.0 / factor;
  SMetric3 m1(1. / h0 / h0, 1. / h1 / h1, 1.0, v0, v1, n);

  SVector3 v4(2., 1., 0.);
  v4.normalize();
  SVector3 v5 = crossprod(v0, n);
  v5.normalize();
  double h4 = 1.0 / factor;
  double h5 = 2.0 / factor;
  SMetric3 m2(1. / h4 / h4, 1. / h5 / h5, 1.0, v4, v5, n);
  // SMetric3 m2 = m1;

  SVector3 v2(1., 1., 0.);
  v2.normalize();
  SVector3 v3 = crossprod(v0, n);
  v3.normalize();
  double h2 = 1.0 / factor;
  double h3 = 2.0 / factor;
  SMetric3 m3(1. / h2 / h2, 1. / h3 / h3, 1.0, v2, v3, n);
  // SMetric3 m4 = m3;

  // Plot metrics at corners
  getEllipsePoints(m1(0, 0), 2.0 * m1(0, 1), m1(1, 1), 0., 0., x, y);
  for(int i = 0; i < nt; ++i) {
    if(i != nt - 1) {
      fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
              x[i + 1], y[i + 1], 0., color, color);
    } else {
      fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
              x[0], y[0], 0., color, color);
    }
  }
  getEllipsePoints(m2(0, 0), 2.0 * m2(0, 1), m2(1, 1), 1., 0., x, y);
  for(int i = 0; i < nt; ++i) {
    if(i != nt - 1) {
      fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
              x[i + 1], y[i + 1], 0., color, color);
    } else {
      fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
              x[0], y[0], 0., color, color);
    }
  }
  getEllipsePoints(m3(0, 0), 2.0 * m3(0, 1), m3(1, 1), 0., 1., x, y);
  for(int i = 0; i < nt; ++i) {
    if(i != nt - 1) {
      fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
              x[i + 1], y[i + 1], 0., color, color);
    } else {
      fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
              x[0], y[0], 0., color, color);
    }
  }
  // SMetric3 mIntersection;
  // SMetric3 mIntersection = intersection_conserve_mostaniso(m3,m3);
  // getEllipsePoints(mIntersection(0,0), 2.0*mIntersection(0,1), mIntersection(1,1),  0.,  1., x,
  // y); for(int i = 0; i < nt; ++i){
  //   if(i != nt-1){
  //     fprintf(myfile,"SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
  //     x[i+1], y[i+1], 0., 2.0*color, 2.0*color);
  //   } else{
  //     fprintf(myfile,"SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
  //     x[0], y[0], 0., 2.0*color, 2.0*color);
  //   }
  // }
  // mIntersection = intersection_alauzet(m3,m3);
  // getEllipsePoints(mIntersection(0,0), 2.0*mIntersection(0,1), mIntersection(1,1),  0.,  1., x,
  // y); for(int i = 0; i < nt; ++i){
  //   if(i != nt-1){
  //     fprintf(myfile,"SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
  //     x[i+1], y[i+1], 0., 2.0*color, 2.0*color);
  //   } else{
  //     fprintf(myfile,"SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
  //     x[0], y[0], 0., 2.0*color, 2.0*color);
  //   }
  // }
  // mIntersection = intersection(m3,m3);
  // getEllipsePoints(mIntersection(0,0), 2.0*mIntersection(0,1), mIntersection(1,1),  0.,  1., x,
  // y); for(int i = 0; i < nt; ++i){
  //   if(i != nt-1){
  //     fprintf(myfile,"SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
  //     x[i+1], y[i+1], 0., 4.0*color, 4.0*color);
  //   } else{
  //     fprintf(myfile,"SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
  //     x[0], y[0], 0., 4.0*color,4.0*color);
  //   }
  // }
  // getEllipsePoints(m1(0,0), 2.0*m1(0,1), m1(1,1), 0., 1., x, y);
  // for(int i = 0; i < nt; ++i){
  //   if(i != nt-1){
  //     fprintf(myfile,"SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
  //     x[i+1], y[i+1], 0., 3.0*color, 3.0*color);
  //   } else{
  //     fprintf(myfile,"SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
  //     x[0], y[0], 0., 3.0*color, 3.0*color);
  //   }
  // }

  SMetric3 mii;

  for(int i = 0; i < nPoints; ++i) {
    for(int j = 0; j < nPoints; ++j) {
      if(xP[nPoints * i + j] + yP[nPoints * i + j] <= 1.0) {
        fprintf(myfile, "SP(%.16g,%.16g,%.16g){%.16g};\n", xP[nPoints * i + j], yP[nPoints * i + j],
                0., color);

        // Interpolate ellipse
        mii = interpolateMetricTriangle(m1, m2, m3, xP[nPoints * i + j], yP[nPoints * i + j]);

        getEllipsePoints(mii(0, 0), 2.0 * mii(0, 1), mii(1, 1), xP[nPoints * i + j],
                         yP[nPoints * i + j], x, y);

        // Plot the ellipse
        for(int ii = 0; ii < nt; ++ii) {
          if(ii != nt - 1) {
            fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[ii],
                    y[ii], 0., x[ii + 1], y[ii + 1], 0., color, color);
          } else {
            fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[ii],
                    y[ii], 0., x[0], y[0], 0., color, color);
          }
        }

        SMetric3 mfoo1 = mii;
        SMetric3 mfoo2 = mii;

        SMetric3 mfoo = intersection(mfoo1, mfoo2);
        getEllipsePoints(mfoo(0, 0), 2.0 * mfoo(0, 1), mfoo(1, 1), xP[nPoints * i + j],
                         yP[nPoints * i + j], x, y);

        // Plot the ellipse
        for(int ii = 0; ii < nt; ++ii) {
          if(ii != nt - 1) {
            fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[ii],
                    y[ii], 0., x[ii + 1], y[ii + 1], 0., color, color);
          } else {
            fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[ii],
                    y[ii], 0., x[0], y[0], 0., color, color);
          }
        }
      }
    }
  }

  fprintf(myfile, "};");
  fclose(myfile);

  return 0;
}

static SMetric3 gradationAlauzet(double grad, const SMetric3 &mp, const SMetric3 &mq, double xp,
                                 double yp, double xq, double yq) {
  double lpq = sqrt((xp - xq) * (xp - xq) + (yp - yq) * (yp - yq));
  fullMatrix<double> V(3, 3);
  fullVector<double> S(3);
  mq.eig(V, S, false);
  double eta1sq = (1. + sqrt(S(0)) * lpq * log(grad));
  eta1sq = 1. / (eta1sq * eta1sq);
  double eta2sq = (1. + sqrt(S(1)) * lpq * log(grad));
  eta2sq = 1. / (eta2sq * eta2sq);
  double eta3sq = (1. + sqrt(S(2)) * lpq * log(grad));
  eta3sq = 1. / (eta3sq * eta3sq);
  SMetric3 N(eta1sq, eta2sq, eta3sq, E_X, E_Y, E_Z);
  SMetric3 L(S(0), S(1), S(2), E_X, E_Y, E_Z);
  // Compute V^T*N*Diag(lambda)*V
  STensor3 R, Rt;
  R.setMat(V);
  Rt = R.transpose();
  Rt *= N;
  Rt *= L;
  Rt *= R;
  Rt.set_m33(1.0); // Enforce M_33 = 1
  // Check result is symmetric
  if(fabs(Rt(0, 1) - Rt(1, 0)) > 1e-6) {
    Rt.print("Rt");
    exit(-1);
  }
  // Transfer result from STensor3 to SMetric3
  Rt.getMat(V);
  SMetric3 tmp;
  tmp.setMat(V);

  return tmp;

  // return intersection_reductionSimultanee(tmp, mp);
}

static SMetric3 gradationAlauzetExplicite(double grad, const SMetric3 &mp, const SMetric3 &mq,
                                          double xp, double yp, double xq, double yq) {
  double lpq = sqrt((xp - xq) * (xp - xq) + (yp - yq) * (yp - yq));
  double lxyGrad = lpq * log(grad);
  double a1 = mq(0, 0);
  double b1 = mq(0, 1);
  double c1 = mq(1, 1);

  double a =
    (pow(
       (c1 / b1 - (a1 / 2 + c1 / 2 - sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2) / b1),
       2) *
     (a1 / 2 + c1 / 2 - sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2)) /
      (pow((lxyGrad *
              sqrt(a1 / 2 + c1 / 2 - sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2) +
            1),
           2) *
       (pow(fabs(c1 / b1 -
                 (a1 / 2 + c1 / 2 - sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2) / b1),
            2) +
        1)) +
    (pow(
       (c1 / b1 - (a1 / 2 + c1 / 2 + sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2) / b1),
       2) *
     (a1 / 2 + c1 / 2 + sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2)) /
      (pow((lxyGrad *
              sqrt(a1 / 2 + c1 / 2 + sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2) +
            1),
           2) *
       (pow(fabs(c1 / b1 -
                 (a1 / 2 + c1 / 2 + sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2) / b1),
            2) +
        1));

  double b =
    -((c1 / b1 - (a1 / 2 + c1 / 2 - sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2) / b1) *
      (a1 / 2 + c1 / 2 - sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2)) /
      (pow((lxyGrad *
              sqrt(a1 / 2 + c1 / 2 - sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2) +
            1),
           2) *
       (pow(fabs(c1 / b1 -
                 (a1 / 2 + c1 / 2 - sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2) / b1),
            2) +
        1)) -
    ((c1 / b1 - (a1 / 2 + c1 / 2 + sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2) / b1) *
     (a1 / 2 + c1 / 2 + sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2)) /
      (pow((lxyGrad *
              sqrt(a1 / 2 + c1 / 2 + sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2) +
            1),
           2) *
       (pow(fabs(c1 / b1 -
                 (a1 / 2 + c1 / 2 + sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2) / b1),
            2) +
        1));

  double c =
    (a1 / 2 + c1 / 2 - sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2) /
      (pow((lxyGrad *
              sqrt(a1 / 2 + c1 / 2 - sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2) +
            1),
           2) *
       (pow(fabs(c1 / b1 -
                 (a1 / 2 + c1 / 2 - sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2) / b1),
            2) +
        1)) +
    (a1 / 2 + c1 / 2 + sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2) /
      (pow((lxyGrad *
              sqrt(a1 / 2 + c1 / 2 + sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2) +
            1),
           2) *
       (pow(fabs(c1 / b1 -
                 (a1 / 2 + c1 / 2 + sqrt(a1 * a1 - 2 * a1 * c1 + 4 * b1 * b1 + c1 * c1) / 2) / b1),
            2) +
        1));

  SMetric3 m(1.0);
  fullMatrix<double> mat(3, 3);
  mat.set(0, 0, a);
  mat.set(0, 1, b);
  mat.set(1, 0, b);
  mat.set(1, 1, c);
  mat.set(0, 2, 0.0);
  mat.set(2, 0, 0.0);
  mat.set(1, 2, 0.0);
  mat.set(2, 1, 0.0);
  mat.set(2, 2, 1.0);
  m.setMat(mat);
  return m;
}

static SMetric3 gradationAlauzetExpliciteDiagonale(double grad, const SMetric3 &mp,
                                                   const SMetric3 &mq, double xp, double yp,
                                                   double xq, double yq) {
  double lpq = sqrt((xp - xq) * (xp - xq) + (yp - yq) * (yp - yq));
  double lxyGrad = lpq * log(grad);
  std::cout << lxyGrad << std::endl;
  double a1 = mq(0, 0);
  double b1 = mq(0, 1);
  double c1 = mq(1, 1);
  double a = a1 / pow((sqrt(a1) * lxyGrad + 1), 2);
  double b = 0.0;
  double c = c1 / pow((sqrt(c1) * lxyGrad + 1), 2);
  SMetric3 m(0.0);
  fullMatrix<double> mat(3, 3);
  mat.set(0, 0, a);
  mat.set(0, 1, b);
  mat.set(1, 0, b);
  mat.set(1, 1, c);
  mat.set(0, 2, 0.0);
  mat.set(2, 0, 0.0);
  mat.set(1, 2, 0.0);
  mat.set(2, 1, 0.0);
  mat.set(2, 2, 1.0);
  m.setMat(mat);
  return m;
}

static SMetric3 gradationAlauzet2(double grad, const SMetric3 &mp, const SMetric3 &mq, double xp,
                                  double yp, double xq, double yq) {
  SVector3 pq(xp - xq, yp - yq, 0.);
  double eta1sq = (1. + sqrt(dot(pq, mq, pq)) * log(grad));
  eta1sq = 1. / (eta1sq * eta1sq);

  SMetric3 tmp(mq);
  tmp *= eta1sq;

  return tmp;

  // return intersection_reductionSimultanee(tmp, mp);
}

static void drawEllipse(SMetric3 &m, double xC, double yC, FILE *file, double color) {
  int nt = 30;
  std::vector<double> x(nt, 0.);
  std::vector<double> y(nt, 0.);
  getEllipsePoints(m(0, 0), 2.0 * m(0, 1), m(1, 1), xC, yC, x, y);
  for(int i = 0; i < nt; ++i) {
    if(i != nt - 1) {
      fprintf(file, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
              x[i + 1], y[i + 1], 0., color, color);
    } else {
      fprintf(file, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i], 0.,
              x[0], y[0], 0., color, color);
    }
  }
}

int mainGradation(int argc, char **argv) {
  double color = 10.0;

  SVector3 v0(1., 1., 0.);
  v0.normalize();
  SVector3 v1(-1., 1., 0.);
  v1.normalize();
  SVector3 v2 = E_Z;
  SMetric3 m = SMetric3(0.004, 1.0, 1., v0, v1, E_Z);
  // SMetric3 m = SMetric3(0.004, 1.0, 1., E_X, E_Y, E_Z);

  // Draw the initial metric
  FILE *myfile = fopen("gradationTest.pos", "w");
  fprintf(myfile, "View \"gradationTest\"{\n");
  double factor = 50.0;
  double x1 = 0.0;
  double y1 = 0.0;
  drawEllipse(m, x1, y1, myfile, color);

  // int nPoints = 7;
  // double L = 10.0;
  // std::vector<double> xP(nPoints * nPoints, 0.);
  // std::vector<double> yP(nPoints * nPoints, 0.);

  // for(int i = 0; i < nPoints; ++i) {
  //   for(int j = 0; j < nPoints; ++j) {
  //     xP[nPoints * i + j] = i * L / nPoints;
  //     yP[nPoints * i + j] = j * L / nPoints;
  //   }
  // }

  // for(int i = 0; i < nPoints; ++i) {
  //   for(int j = 0; j < nPoints; ++j) {
  //     // Grow the metric to (xP,yP) and draw it
  //     double x2 = xP[nPoints * i + j], y2 = yP[nPoints * i + j];
  //     SMetric3 m2 = gradationAlauzet2(1.05, SMetric3(1.0), m, x1, y1, x2, y2);

  //     getEllipsePoints(m2(0, 0), 2.0 * m2(0, 1), m2(1, 1), x2, y2, x, y);
  //     for(int i = 0; i < nt; ++i) {
  //       if(i != nt - 1) {
  //         fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i],
  //         0.,
  //                 x[i + 1], y[i + 1], 0., color, color);
  //       } else {
  //         fprintf(myfile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g};\n", x[i], y[i],
  //         0.,
  //                 x[0], y[0], 0., color, color);
  //       }
  //     }
  //   }
  // }

  SMetric3 m2 = SMetric3(0.004, 0.004, 1., E_X, E_Y, E_Z);
  double x2 = 10.0, y2 = 0.0;
  // drawEllipse(m2, x2, y2, myfile, color);

  m.print("input");
  SMetric3 mFrom1At2 = gradationAlauzet(1.2, SMetric3(10000.0), m, x1, y1, x2, y2);
  drawEllipse(mFrom1At2, x2, y2, myfile, 2 * color);

  m.print("input");
  SMetric3 mFrom1At2Sym = gradationAlauzetExplicite(1.2, SMetric3(10000.0), m, x1, y1, x2, y2);
  drawEllipse(mFrom1At2Sym, x2, y2, myfile, 3 * color);

  mFrom1At2.print("reg");
  mFrom1At2Sym.print("sym");

  // SMetric3 m2Inter = intersection_reductionSimultanee(m2, mFrom1At2);
  // drawEllipse(m2Inter, x2, y2, myfile, color);

  SMetric3 mFrom2At1 = gradationAlauzet(1.2, m, m2, x1, y1, x2, y2);
  // drawEllipse(mFrom2At1, x1, y1, myfile, 3*color);

  // SMetric3 m1Inter = intersection_reductionSimultanee(m, mFrom2At1);
  // drawEllipse(m1Inter, x1, y1, myfile, color);


  // SMetric3 M1(1.0);
  // M1.set_m11(2.501000000000E+00);
  // M1.set_m22(2.501000000000E+00);
  // M1.set_m33(1.);
  // M1.set_m21(-2.499000000000E+00);
  // M1.set_m31(0.);
  // M1.set_m32(0.);
  // M1.print("M1");
  // SMetric3 M1(1.0);
  // M1.set_m11(1.464029283623E+00);
  // M1.set_m22(1.496992768857E+00);
  // M1.set_m33(1.);
  // M1.set_m21(-1.449113059614E+00);
  // M1.set_m31(0.);
  // M1.set_m32(0.);
  // M1.print("M1");
  SMetric3 M1(1.0);
  M1.set_m11(2);
  M1.set_m22(2);
  M1.set_m33(1.);
  M1.set_m21(1);
  M1.set_m31(0.);
  M1.set_m32(0.);
  M1.print("M1");

  drawEllipse(M1, 50., 10., myfile, 2*color);

  // SMetric3 M2(1.0);
  // M2.set_m11(5.002000000000E+00);
  // M2.set_m22(5.002000000000E+00);
  // M2.set_m33(1.);
  // M2.set_m21(-4.998000000000E+00);
  // M2.set_m31(0.);
  // M2.set_m32(0.);
  // M2.print("M2");
  // SMetric3 M2(1.0);
  // M2.set_m11(2.458930005425E+00);
  // M2.set_m22(2.523713604758E+00);
  // M2.set_m33(1.);
  // M2.set_m21(-2.455485150437E+00);
  // M2.set_m31(0.);
  // M2.set_m32(0.);
  // M2.print("M2");
  SMetric3 M2(1.0);
  M2.set_m11(4);
  M2.set_m22(4);
  M2.set_m33(1.);
  M2.set_m21(2);
  M2.set_m31(0.);
  M2.set_m32(0.);
  M2.print("M2");

  drawEllipse(M2, 50., 10., myfile, .5*color);

  // SMetric3 M3(1.0);
  // M3.set_m11(5.002000000000E+00);
  // M3.set_m22(5.002000000000E+00);
  // M3.set_m33(1.);
  // M3.set_m21(0.);
  // M3.set_m31(0.);
  // M3.set_m32(0.);
  // M3.print("M3");

  // drawEllipse(M3, 50., 10., myfile, .5*color);

  fprintf(myfile, "};");
  fclose(myfile);
}

int main(int argc, char **argv) {
  // mainOctreeInterpolation(argc, argv);
  // mainTriangleInterpolation(argc, argv);
  mainGradation(argc, argv);
  return 0;
}
