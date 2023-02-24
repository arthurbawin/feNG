#ifndef _FEMETRICTOOLS_
#define _FEMETRICTOOLS_

#include "feRecovery.h"
#include "STensor3.h"

#include "../contrib/Eigen/Eigen"

#if defined(HAVE_SOPLEX)
#include "soplex.h"

typedef struct linearProblemStruct {
  soplex::SoPlex problem;
  soplex::SPxSolver::Status stat;
  soplex::DSVector row;
  soplex::LPRow lprow;
  soplex::LPRowSet lprowset;
  soplex::DVector prim;
  std::vector<double> lvl1;
  std::vector<double> constraints;
} linearProblem;
#endif

struct gmshEdgeLessThan {
  bool operator()(const std::pair<int, int> &e1, const std::pair<int, int> &e2) const
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

double evaluateFieldFromRecovery(int indexDerivative, feRecovery *rec, double *x);

double evaluateFieldFromRecoveryCallback(int indexDerivative, void *recUserPtr, double *x);

double f(feRecovery *rec, double *x);
double fx(feRecovery *rec, double *x);
double fy(feRecovery *rec, double *x);
double fxx(feRecovery *rec, double *x);
double fxy(feRecovery *rec, double *x);
double fyx(feRecovery *rec, double *x);
double fyy(feRecovery *rec, double *x);
double fxxx(feRecovery *rec, double *x);
double fxxy(feRecovery *rec, double *x);
double fxyx(feRecovery *rec, double *x);
double fxyy(feRecovery *rec, double *x);
double fyxx(feRecovery *rec, double *x);
double fyxy(feRecovery *rec, double *x);
double fyyx(feRecovery *rec, double *x);
double fyyy(feRecovery *rec, double *x);

void computeDirectionFieldFromGradient(double *x, double &C, double &S, double tol, feRecovery *rec,
                                       FILE *F_grad = nullptr, FILE *F_iso = nullptr);
void computeDirectionFieldFromHessian(double *x, double &C, double &S, double tol, feRecovery *rec,
                                      FILE *F = nullptr);

void smoothDirections(std::map<size_t, double> &C, std::map<size_t, double> &S, int nIter = 100,
                      double tol = 0.1);

void smoothSizes(std::map<size_t, double> &L1, std::map<size_t, double> &L2, int nIter, double tol);

double dtt(double *x, double C, double S, feRecovery *rec);
double dttt(double *x, double C, double S, feRecovery *rec, int direction);
double dtttAnalytical(double *x, double C, double S, feRecovery *rec, int direction);

void computeWorstMetric(int nTheta, int nIncr, double e, double *x, double C, double S,
                        double &h1Opt, double &h2Opt, feRecovery *rec, SMetric3 &M, double lMin,
                        double lMax);

void gradationMetriques(double gradation, int nmax, std::vector<double> &coord,
                        std::map<int, SMetric3> &metricsOnGmshModel);

void metricHechtKuate(int nbpoints, double *x, double *y, double &A, double &B, double &C,
                      double epsilon, double *xNew, double *yNew, int dimRecovery,
                      int degreeSolution);

#if defined(HAVE_SOPLEX)
bool computeMetricLogSimplexCurved(double *x, double cG, double sG, feRecovery *rec,
                                   Eigen::Matrix2d &Q, int maxIter, int nThetaPerQuadrant,
                                   double tol, int &numIter, linearProblem &myLP);
#endif

int computePointsUsingScaledCrossFieldPlanarP2(
  const char *modelForMetric, const char *modelForMesh, int VIEW_TAG, int faceTag,
  std::vector<double> &pts, double er(double *, double *, double *, double *, double *, double *),
  bool inside(double *, bool), double pointwiseError(double *), int onlyGenerateVertices,
  double evaluateFieldFromRecovery(int, void *, double *), void *recoveryUserPointer,

  void interpolateMetricP1WithDerivativesWrapper(void *, const double *, Eigen::Matrix2d &,
                                                 Eigen::Matrix2d &, Eigen::Matrix2d &),

  void interpolateMetricP1Wrapper(void *, const double *, Eigen::Matrix2d &, Eigen::Matrix2d &,
                                  Eigen::Matrix2d &),

  void interpolateMetricAndDerivativeOnP2EdgeWrapper(
    void *, const double, const Eigen::Matrix2d &, const Eigen::Matrix2d &, const Eigen::Matrix2d &,
    const Eigen::Matrix2d &, const Eigen::Matrix2d &, Eigen::Matrix2d &, Eigen::Matrix2d &,
    Eigen::Matrix2d &),

  void interpolateMetricP1Wrapper1D(void *, const double *, const double *, Eigen::Matrix2d &,
                                    Eigen::Matrix2d &),

  void interpolateMetricAndDerivativeOnP2EdgeWrapper1D(
    void *, const double, const Eigen::Matrix2d &, const Eigen::Matrix2d &, const Eigen::Matrix2d &,
    const Eigen::Matrix2d &, Eigen::Matrix2d &, Eigen::Matrix2d &),

  void *metric);

void fooTestStuff();

#endif