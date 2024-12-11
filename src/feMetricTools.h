#ifndef _FEMETRICTOOLS_
#define _FEMETRICTOOLS_

#include "feMatrixInterface.h"
#include "feRecovery.h"
#include "feNewRecovery.h"
#include "STensor3.h"

#include <Eigen/Eigen>

#if defined(HAVE_GMSH)
#include "gmsh.h"
#include "curvedMesh_structs.h"
#endif

#if defined(HAVE_SOPLEX)
#include "soplex.h"

typedef struct {
  soplex::SoPlex problem;
  soplex::SPxSolver::Status stat;
  soplex::DSVector row;
  // soplex::LPRow lprow;
  soplex::LPRowSet lprowset;
  soplex::DVector prim;
  size_t numConstraints;
  std::vector<double> lvl1;
  std::vector<double> constraints;
  MetricTensor Hij;
  bool uniformErrorCurve;
  int numLoopsUniformErrorCurve;
  Eigen::VectorXd lvl1PolynomialCoeffs_minus;
  Eigen::VectorXd lvl1PolynomialCoeffs_plus;
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

bool computeAnalyticMetricP2ForLpNorm(std::vector<double> &errorCoeff, MetricTensor &Q, const double maxDiameter);
bool computeAnalyticMetricP2ForH1semiNorm(std::vector<double> &errorCoeff, MetricTensor &Q);

void computeDirectionFieldFromGradient(const double x[2], const int vertex, const double firstDerivatives[2], const double tol, double grad[2]);

void solveSizePolynomialLinear(const double targetError, const double directionGrad[2], const double du[2], const double d2u[4], double &Liso, double &Lgrad);
void solveSizePolynomialQuadratic(const double targetError, const double directionGrad[2], const double du[2], const double d2u[4], const double d3u[8], double &Liso, double &Lgrad);

double evaluateFieldFromRecovery(int indexDerivative, feRecovery *rec, double *x);

double evaluateFieldFromRecoveryCallback(int indexDerivative, void *recUserPtr, double *x);

void computeDirectionFieldFromGradient(double *x, double &C, double &S, double tol, feRecovery *rec,
                                       FILE *F_grad = nullptr, FILE *F_iso = nullptr);
void computeDirectionFieldFromHessian(double *x, double &C, double &S, double tol, feRecovery *rec,
                                      FILE *F = nullptr);

void smoothDirections(std::map<size_t, double> &C, std::map<size_t, double> &S, int nIter = 100,
                      double tol = 0.1, bool plot = false);

void smoothDirections(const std::vector<std::size_t> &nodeTags,
                      const std::vector<double> &coord,
                      std::map<int, MetricTensor> &metrics, int nIter = 100, double tol = 0.1);

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
bool computeMetricLogSimplexStraight(const double *x, const std::vector<double> &errorCoefficients,
                                     const int degree, const int nThetaPerQuadrant,
                                     const int maxIter, const double tol, 
                                     const bool applyBinomialCoefficients, MetricTensor &Qres,
                                     int &numIter, linearProblem &myLP);

bool logSimplexOneIteration(const double *x, const std::vector<double> &errorCoefficients,
                            const int degree, const int nThetaPerQuadrant, const bool applyBinomialCoefficients,
                            MetricTensor &Qres, linearProblem &myLP);

bool computeMetricLogSimplexCurved(const int vertex, double *x, double directionGradient[2], feNewRecovery *rec,
                                   MetricTensor &Qres, int maxIter, int nThetaPerQuadrant,
                                   double tol, int &numIter, linearProblem &myLP);
#endif

void reconstructManifoldFromMetric(const std::vector<std::size_t> &nodeTags,
  const std::vector<double> &coord,
  const std::map<int, MetricTensor> &metricsAtNodeTags,
  const std::vector<int> sequentialTag2nodeTag,
  feNewRecovery *recoveredField);

void computeGraphSurfaceMetric(const double *x, const double *df, const double *d2f, const double *d3f, MetricTensor &M, FILE *file);

#endif