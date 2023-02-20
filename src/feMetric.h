#ifndef _FEMETRIC_
#define _FEMETRIC_

#include "feRecovery.h"
// #include "eigen3/Eigen/Eigen"
#include "../contrib/Eigen/Eigen"

#include "STensor3.h"

typedef struct feMetricOptions {
  bool debug = true;

  int computationMethod;
  int polynomialDegree;
  int nTargetVertices;
  double eTargetError;
  double hMin;
  double hMax;
  double LpNorm;
  int nPhi;

  // The name of the Gmsh model on which the metric field will be computed
  std::string modelForMetric;

  std::string modelForMesh;

  bool isGmshModelReady = false;
  // Geometric data to create the model with Gmsh api
  double modelSizeX;
  double modelSizeY;
  int directionFieldFromDerivativesOfOrder;

  bool (*inside)(double *, bool);

  // MMG parameters
  std::string meshName;
  std::string adaptedMeshName;
  std::string metricMeshNameForMMG;
  std::string metricMeshNameForMMG_out;
  std::string recoveryName;

  bool enableGradation;
  double gradation;

} feMetricOptions;

class feMetric
{
protected:
  feRecovery *_recovery;
  feMetricOptions _options;

  std::map<int, SMetric3> _metrics;
  std::map<int, SMetric3> _metricsOnGmshModel;
  std::map<int, Eigen::Matrix2d> _metricsOnGmshModel_eigen;

  std::map<Vertex *, int> _v2n;

  // The tag of the gmsh view in which the metric field is stored (if using Gmsh)
  int _metricViewTag = -1;

public:
  feMetric(feRecovery *recovery, feMetricOptions metricOptions);
  ~feMetric() {}

  void setGmshMetricModel(std::string metricModel) { _options.modelForMetric = metricModel; }

  void computeMetrics();
  void computeMetricsHechtKuate();
  void computeMetricsWithDirectionField();
  void computeMetricsLogSimplex();
  void computeMetricsExtremeSizesOnly();

  void setMetricViewTag(int tag) { _metricViewTag = tag; }
  int getMetricViewTag() { return _metricViewTag; }

  void metricScaling();
  void metricScalingFromGmshSubstitute();

  // Intersecter avec un autre feMetric
  void writeSizeFieldSol2D(std::string solFileName);
  void writeSizeFieldSol3D(std::string solFileName);
  void writeSizeFieldGmsh(std::string meshName, std::string metricMeshName);
  void drawEllipsoids(std::string posFile);

  // With gradient wrt to physical coordinates x,y
  void interpolateMetricP1WithDerivatives(const double *x, Eigen::Matrix2d &M,
                                          Eigen::Matrix2d &dMdx, Eigen::Matrix2d &dMdy);

  // With gradient of M wrt 2d position alpha (a1,a2)
  void interpolateMetricP1(const double *x, Eigen::Matrix2d &M, Eigen::Matrix2d &sumduda1M,
                           Eigen::Matrix2d &sumduda2M);
  void interpolateMetricAndDerivativeOnP2Edge(const double t, const Eigen::Matrix2d &M11,
                                              const Eigen::Matrix2d &M20,
                                              const Eigen::Matrix2d &M02,
                                              const Eigen::Matrix2d &sumduda1M,
                                              const Eigen::Matrix2d &sumduda2M, Eigen::Matrix2d &M,
                                              Eigen::Matrix2d &dMda1, Eigen::Matrix2d &dMda2);

  // With derivative of M wrt 1D alpha (dMda)
  void interpolateMetricP1(const double *x, const double *gammaOrth, Eigen::Matrix2d &M,
                           Eigen::Matrix2d &sumdudaM);
  void interpolateMetricAndDerivativeOnP2Edge(double t, const Eigen::Matrix2d &M11,
                                              const Eigen::Matrix2d &M20,
                                              const Eigen::Matrix2d &M02,
                                              const Eigen::Matrix2d &sumdudaM, Eigen::Matrix2d &M,
                                              Eigen::Matrix2d &dMda);
};

void interpolateMetricP1WithDerivativesWrapper(void *metric, const double *x, Eigen::Matrix2d &M,
                                               Eigen::Matrix2d &dMdx, Eigen::Matrix2d &dMdy);

void interpolateMetricP1Wrapper(void *metric, const double *x, Eigen::Matrix2d &M,
                                Eigen::Matrix2d &sumduda1M, Eigen::Matrix2d &sumduda2M);

void interpolateMetricP1Wrapper1D(void *metric, const double *x, const double *gammaOrth,
                                  Eigen::Matrix2d &M, Eigen::Matrix2d &sumdudaM);

void interpolateMetricAndDerivativeOnP2EdgeWrapper(
  void *metric, const double t, const Eigen::Matrix2d &M11, const Eigen::Matrix2d &M20,
  const Eigen::Matrix2d &M02, const Eigen::Matrix2d &sumduda1M, const Eigen::Matrix2d &sumduda2M,
  Eigen::Matrix2d &M, Eigen::Matrix2d &dMda1, Eigen::Matrix2d &dMda2);

void interpolateMetricAndDerivativeOnP2EdgeWrapper1D(void *metric, const double t,
                                                     const Eigen::Matrix2d &M11,
                                                     const Eigen::Matrix2d &M20,
                                                     const Eigen::Matrix2d &M02,
                                                     const Eigen::Matrix2d &sumdudaM,
                                                     Eigen::Matrix2d &M, Eigen::Matrix2d &dMda);

#endif