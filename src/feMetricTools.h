#ifndef _FEMETRICTOOLS_
#define _FEMETRICTOOLS_

#include "feRecovery.h"
#include "STensor3.h"

struct gmshEdgeLessThan {
  bool operator()(const std::pair<int,int> &e1, const std::pair<int,int> &e2) const
  {
    int tag10 = e1.first, tag11 = e1.second;
    int tag20 = e2.first, tag21 = e2.second;
    int diffMin = fmin(tag10, tag11) - fmin(tag20, tag21);
    int diffMax = fmax(tag10, tag11) - fmax(tag20, tag21);
    if(diffMin < 0) return true;
    if(diffMin > 0) return false;
    if(diffMax < 0) return true;
    if(diffMax > 0) return false;
    return false;
  }
};

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

void smoothSizes(std::map<size_t, double> &L1, std::map<size_t, double> &L2, int nIter,
                      double tol);

double dtt(const double x, const double y, double C, double S, feRecovery *rec);
double dttt(const double x, const double y, double C, double S, feRecovery *rec, int direction);
double dtttAnalytical(const double x, const double y, double C, double S, feRecovery *rec, int direction);

void computeWorstMetric(int nTheta, int nIncr, double e, double x, double y, double C, double S, 
  double &h1Opt, double &h2Opt, feRecovery *rec, SMetric3 &M, double lMin, double lMax);

void gradationMetriques(double gradation, std::vector<double> &coord, std::map<int, SMetric3> &metricsOnGmshModel);

void metricHechtKuate(int nbpoints, double *x, double *y, double &A, double &B, double &C,
                      double epsilon, double *xNew, double *yNew, int dimRecovery,
                      int degreeSolution);

int computePointsUsingScaledCrossFieldPlanarP2(
  const char *modelForMetric, const char *modelForMesh, int VIEW_TAG, int faceTag,
  std::vector<double> &pts, double er(double *, double *, double *, double *, double *, double *),
  bool inside(double *),
  void grad_er(double *, double *, double *, double *, double *, double *, std::vector<double> &),
  double pointwiseError(double *),
  int onlyGenerateVertices);

#endif