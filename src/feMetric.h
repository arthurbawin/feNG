#ifndef _FEMETRIC_
#define _FEMETRIC_

#include "feRecovery.h"
#include "feNewRecovery.h"

#include "../contrib/Eigen/Eigen"

#include "STensor3.h"

enum class adaptationMethod { ANISO_P1, ANISO_PN, CURVED_LS, CURVED_EXTREME_SIZES };

class feMetricOptions
{
public:
  // These members are default-initialized but should be changed:
  //
  // Basic parameters:
  //
  // Min/max allowed mesh size
  double hMin = 1e-12;
  double hMax = 1.;
  // Mesh adaptation method
  adaptationMethod method = adaptationMethod::ANISO_P1;
  // Degree of the polynomial solution for the field driving
  // the mesh adaptation
  int polynomialDegree = 1;
  // Target number of vertices in adapted mesh
  // Used for anisotropic adaptation and curved adaptation
  // with the log-simplex method, ignored otherwise
  int nTargetVertices = 100;
  // Target interpolation error in adapted mesh
  // Used for curved adaptation along iso/grad, ignored otherwise
  double eTargetError = 1e-2;
  // Limit variations in the metric field for smoother mesh
  bool enableGradation = false;
  // Maximum growth in the metric field if gradation is enabled
  // See F. Alauzet, Size gradation control of anisotropic meshes
  double gradation = 1.5;
  // Target Lp norm in which the interpolation error is minimized
  // See e.g. Alauzet & Loseille, Multi-Dimensional Continuous Metric
  // for Mesh Adaptation
  double LpNorm = 2.;
  // Number of cycles to produce the final anisotropic mesh
  // One cycle = evaluate metric field, then remesh
  int nLoopsAnisoMesh = 1;

  //  to create the model with Gmsh api
  double modelSizeX = 1.;
  double modelSizeY = 1.;

  // MMG parameters
  std::string adaptedMeshName = "";
  std::string mmgInputMeshfile = "inputMMG.msh";
  std::string mmgOutputMeshfile = "outputMMG.msh";
  std::string recoveryName = "recoveredSolution.msh";

  // Print a bunch of ellipses and other visualisation files if set to true
  bool debug = true;

  //
  // More advanced parameters:
  //
  // Recover principal directions of metric field from
  // derivatives of this order (default is 1: iso/gradient)
  int directionFieldFromDerivativesOfOrder = 1;
  // Number of points for unit circle discretization
  int nPhi = 51;
  // Callback that returns true if a point is inside the
  // domain to be meshed
  bool (*inside)(double *, bool) = nullptr;

public:
  // These members are set automatically and should not be changed
  std::string backgroundMeshfile;
  // Name of the Gmsh model on which the metric field will be computed
  std::string modelForMetric = "";
  // Name of the Gmsh model on which the curved mesh is generated
  std::string modelForMesh = "modelForMesh";
  bool isGmshModelReady = false;

public:
  feMetricOptions(const std::string &meshFile) : backgroundMeshfile(meshFile){};
};

class feMetric
{
protected:
  feRecovery *_recovery;
  feNewRecovery *_newRecovery;
  feMetricOptions _options;

  std::map<int, SMetric3> _metrics;
  std::map<int, SMetric3> _metricsOnGmshModel;

  std::map<int, SMetric3> _metricsP1;
  std::map<int, Eigen::Matrix2d> _metricsP1AtNodeTags;

  // Ultimately, the "good" type: mapping nodeTags to "Eigen-agnostic" MetricTensors
  std::map<int, MetricTensor> _metricTensorAtNodetags;

  std::map<const Vertex *, Eigen::Matrix2d> _metricsOnGmshModelP1;

  std::map<int, Eigen::Matrix2d> _metricsOnGmshModel_eigen;

  // Vertex to nodeTag map (and inverse)
  std::map<const Vertex *, int> _v2n;
  std::map<int, const Vertex *> _n2v;

  // The tag of the gmsh view in which the metric field is stored (if using Gmsh)
  int _metricViewTag = -1;

public:
  feMetric(feRecovery *recovery, feMetricOptions metricOptions);
  feMetric(feNewRecovery *recovery, feMetricOptions metricOptions);
  ~feMetric() {}

  void setGmshMetricModel(std::string metricModel) { _options.modelForMetric = metricModel; }

  feStatus computeMetrics();
  feStatus computeMetricsP1();
  feStatus computeMetricsHechtKuate();
  feStatus computeMetricsLogSimplex();
  feStatus computeMetricsExtremeSizesOnly();

  void setMetricViewTag(int tag) { _metricViewTag = tag; }
  int getMetricViewTag() { return _metricViewTag; }

  void metricScaling();
  void metricScalingFromGmshSubstitute();

  template <class MetricType>
  void metricScalingFromGmshSubstitute(std::map<int, MetricType> &metrics,
                                       const std::vector<size_t> &nodeTags,
                                       double exponentInIntegral, double exponentForDeterminant);

  // Intersecter avec un autre feMetric
  void writeSizeFieldSol2D(std::string solFileName);
  void writeSizeFieldSol3D(std::string solFileName);
  void writeSizeFieldGmsh(std::string meshName, std::string metricMeshName);

  template <class MetricType>
  void drawEllipsoids(const std::string &posFile,
                      const std::map<const Vertex *, MetricType> &metrics, double sizeFactor,
                      int nPoints);

  template <class MetricType>
  void drawEllipsoids(const std::string &posFile, std::map<int, MetricType> &metrics,
                      const std::vector<std::size_t> &nodeTags, const std::vector<double> &coord,
                      double sizeFactor, int nPoints);

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
                                              const Eigen::Matrix2d &sumduda2M, Eigen::Matrix2d
                                              &M, Eigen::Matrix2d &dMda1, Eigen::Matrix2d
                                              &dMda2);

  // With derivative of M wrt 1D alpha (dMda)
  void interpolateMetricP1(const double *x, const double *gammaOrth, Eigen::Matrix2d &M,
                           Eigen::Matrix2d &sumdudaM);
  void interpolateMetricAndDerivativeOnP2Edge(double t, const Eigen::Matrix2d &M11,
                                              const Eigen::Matrix2d &M20,
                                              const Eigen::Matrix2d &M02,
                                              const Eigen::Matrix2d &sumdudaM, Eigen::Matrix2d
                                              &M, Eigen::Matrix2d &dMda);
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
