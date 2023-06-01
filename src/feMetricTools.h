#ifndef _FEMETRICTOOLS_
#define _FEMETRICTOOLS_

#include "feMatrixInterface.h"
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
  size_t numConstraints;
  std::vector<double> lvl1;
  std::vector<double> constraints;
  MetricTensor Hij;
  bool uniformErrorCurve;
} linearProblem;
#endif

typedef struct {
  int nElements;
  std::vector<int> whichElements;
  int iVertexToMove_globalTag;
  std::vector<int> iVertexToMove_localTag;
  double modifiedCoord[2];
  double gradient[2];
  double initialCost;
  double finalCost;
} edgeAndVertexData;

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

template <class MetricType>
void drawSingleEllipse(FILE *file, const double *x, MetricType &M, double sizeFactor, int nPoints);
template <class MetricType>
void drawEllipsoids(const std::string &posFile, std::map<int, MetricType> &metrics,
                    const std::vector<std::size_t> &nodeTags, const std::vector<double> &coord,
                    double sizeFactor, int nPoints);

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
bool computeMetricLogSimplexStraight(const double *x, const std::vector<double> &errorCoefficients,
                                     const int degree, const int nThetaPerQuadrant,
                                     const int maxIter, const double tol, MetricTensor &Qres,
                                     int &numIter, linearProblem &myLP);

bool computeMetricLogSimplexCurved(double *x, double cG, double sG, feRecovery *rec,
                                   Eigen::Matrix2d &Q, int maxIter, int nThetaPerQuadrant,
                                   double tol, int &numIter, linearProblem &myLP);
#endif

int computePointsUsingScaledCrossFieldPlanarP2(
  const char *modelForMetric,
  const char *modelForMesh,
  int VIEW_TAG,
  int faceTag,
  std::vector<double> &pts,
  // double computeInterpolationError(const int whichElements[2],
  //                                  const int iVertexToMove_globalTag,
  //                                  const int iVertexToMove_localTag[2],
  //                                  double *modifiedCoord),
  double computeInterpolationError(const edgeAndVertexData &data, double *modifiedCoord),
  void computeInterpolationErrorGradient(const edgeAndVertexData &data,
                                         const double *modifiedCoord,
                                         double gradient[2]),
  void applyCurvatureToFeMesh(const edgeAndVertexData &data, const double *modifiedCoord),
  void getMidnodeTags(const SPoint2 edge[2], const double tol, int &elementTag, int &localTag, int &globalTag),
  void getPolyMeshVertexTags(const SPoint2 &p, const double tol, 
                             std::vector<int> &elementTags, std::vector<int> &localTags, int &globalTag),
  bool inside(double *, bool), double pointwiseError(double *),
  int onlyGenerateVertices,
  double evaluateFieldFromRecovery(int, void *, double *),
  void *recoveryUserPointer,

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

#endif