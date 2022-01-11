#ifndef _FEMETRICTOOLS_
#define _FEMETRICTOOLS_

#include "feRecovery.h"

double evaluateFieldFromRecovery(int indexDerivative, feRecovery *rec, double x, double y);

double f(feRecovery *rec, double x, double y);
double fx(feRecovery *rec, double x, double y);
double fy(feRecovery *rec, double x, double y);
double fxx(feRecovery *rec, double x, double y);
double fxy(feRecovery *rec, double x, double y);
double fyx(feRecovery *rec, double x, double y);
double fyy(feRecovery *rec, double x, double y);
double fxxx(feRecovery *rec, double x, double y);
double fxxy(feRecovery *rec, double x, double y);
double fxyx(feRecovery *rec, double x, double y);
double fxyy(feRecovery *rec, double x, double y);
double fyxx(feRecovery *rec, double x, double y);
double fyxy(feRecovery *rec, double x, double y);
double fyyx(feRecovery *rec, double x, double y);
double fyyy(feRecovery *rec, double x, double y);

void computeDirectionFieldFromGradient(double x, double y, double &C, double &S, double tol,
                                       feRecovery *rec, FILE *F = nullptr);
void computeDirectionFieldFromHessian(double x, double y, double &C, double &S, double tol,
                                      feRecovery *rec, FILE *F = nullptr);

void smoothDirections(std::map<size_t, double> &C, std::map<size_t, double> &S, FILE *F = nullptr,
                      int nIter = 100, double tol = 0.1);

double dtt(const double x, const double y, double C, double S, feRecovery *rec);
double dttt(const double x, const double y, double C, double S, feRecovery *rec, int direction);

void metricHechtKuate(int nbpoints, double *x, double *y, double &A, double &B, double &C,
                      double epsilon, double *xNew, double *yNew, int dimRecovery,
                      int degreeSolution);

int computePointsUsingScaledCrossFieldPlanarP2(
  const char *modelForMetric, const char *modelForMesh, int VIEW_TAG, int faceTag,
  std::vector<double> &pts, double er(double *, double *, double *, double *, double *, double *),
  bool inside(double *),
  void grad_er(double *, double *, double *, double *, double *, double *, std::vector<double> &));

#endif