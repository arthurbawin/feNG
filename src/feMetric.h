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
  double hMin = 1e-3;
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
  
  // SHOULD NOT BE CHANGED CURRENTLY BECAUSE WE CANNOT LOOP
  // AND ADAPT SEVERAL TIMES WITH THE SAME METRIC FIELD
  // (en même temps est-ce que c'est utile de boucler pour un
  // même champ de métriques? Le maillage obtenu est censé être bon pour ce champ-là)
  int nLoopsAnisoMesh = 1;

  //  to create the model with Gmsh api
  double modelSizeX = 1.;
  double modelSizeY = 1.;

  // MMG parameters
  std::string adaptedMeshName = "adapted.msh";
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

  // Options for the log-simplex method
  struct { 
    // Number of points for the discretization of a quarter of the unit circle
    int nThetaPerQuadrant = 25;
    int maxIter = 50;
    // Tolerance on the difference Q-Qprev
    double tol = 1e-5;
  } logSimplexOptions;
  
  // Callback that returns true if a point is inside the
  // domain to be meshed
  bool (*insideCallback)(double *, bool) = nullptr;

  double userValue;

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
public:
  feMetricOptions _options;
  int _nVerticesPerElmOnBackmesh;

protected:
  feRecovery *_recovery;
  feNewRecovery *_newRecovery;

  std::map<int, SMetric3> _metrics;
  std::map<int, SMetric3> _metricsOnGmshModel;

  std::map<int, SMetric3> _metricsP1;
  std::map<int, Eigen::Matrix2d> _metricsP1AtNodeTags;

  // =========================================================
  // New interface uses only those types:
  // Ultimately, the "good" type: mapping nodeTags to "Eigen-agnostic" MetricTensors
  std::map<int, MetricTensor> _metricTensorAtNodetags;
  std::map<int, SMetric3> _smetric3AtNodetags;
  // =========================================================

  std::map<const Vertex *, Eigen::Matrix2d> _metricsOnGmshModelP1;

  std::map<int, Eigen::Matrix2d> _metricsOnGmshModel_eigen;

  // Vertex to nodeTag map (and inverse)
  // nodeTag is the Gmsh node tag obtained from gmsh::model::mesh::getNodes
  // _nodeTag2sequentialTag maps the nodeTags to the internal tag, which is higher 
  // than the total number of nodeTags when a P2 mesh is used.
  std::map<const Vertex *, int> _v2n;
  std::map<int, const Vertex *> _n2v;
  std::map<int, int> _nodeTag2sequentialTag;
  std::map<int, int> _sequentialTag2nodeTag;

  // The tag of the gmsh view in which the metric field is stored (if using Gmsh)
  int _metricViewTag = -1;

public:
  feMetric(feRecovery *recovery, feMetricOptions metricOptions);
  feMetric(feNewRecovery *recovery, feMetricOptions metricOptions);
  ~feMetric() {}

  void setGmshMetricModel(std::string metricModel) { _options.modelForMetric = metricModel; }
  void setRecovery(feRecovery *recovery){ _recovery = recovery; }

  // If computing metrics on a P1 mesh with derivatives recovered on a P1 mesh,
  // nodeTags is the sorted list of all vertices returned by gmsh (minus 1), 
  // hence a trivial vector from 0 to nVertices-1.
  //
  // If computing metrics on a P1 mesh with derivatives recovered on a P2 mesh,
  // nodeTags is the sorted list of only P1 vertices from the P2 mesh, to access the derivatives
  // at the right indices.
  //
  // If computing metrics on a P2 mesh, nodeTags is also the trivial list of vertices
  // computed by gmsh (minus 1)
  feStatus computeMetrics();

  feStatus computeMetricsP1();
  feStatus computeMetricsHechtKuate();
  feStatus computeMetricsLogSimplex();
  feStatus computeMetricsExtremeSizesOnly();

  // ====================
  // New interface
  feStatus createVertex2NodeMap(std::vector<std::size_t> &nodeTags, std::vector<double> &coord);
  feStatus computeMetricsP1(std::vector<std::size_t> &nodeTags, std::vector<double> &coord);
  feStatus computeMetricsPn(std::vector<std::size_t> &nodeTags, std::vector<double> &coord);
  void applyGradation(std::vector<std::size_t> &nodeTags, std::vector<double> &coord);
  void writeMetricField(std::vector<std::size_t> &nodeTags, std::vector<double> &coord);

  MetricTensor &getMetricAtSequentialTag(int tag){ return _metricTensorAtNodetags[_sequentialTag2nodeTag[tag]]; };
  // ====================

  void setMetricViewTag(int tag) { _metricViewTag = tag; }
  int getMetricViewTag() { return _metricViewTag; }

  void metricScaling();
  void metricScalingFromGmshSubstitute();

  template <class MetricType>
  void metricScalingFromGmshSubstitute(std::map<int, MetricType> &metrics,
                                       const std::vector<size_t> &nodeTags,
                                       double exponentInIntegral,
                                       double exponentForDeterminant);

  // Intersecter avec un autre feMetric
  void writeSizeFieldSol2D(std::string solFileName);
  void writeSizeFieldSol3D(std::string solFileName);
  void writeSizeFieldGmsh(std::string meshName, std::string metricMeshName);

  // template <class MetricType>
  // void drawEllipsoids(const std::string &posFile,
  //                     const std::map<const Vertex *, MetricType> &metrics, double sizeFactor,
  //                     int nPoints);

  // template <class MetricType>
  // void drawEllipsoids(const std::string &posFile, std::map<int, MetricType> &metrics,
  //                     const std::vector<std::size_t> &nodeTags, const std::vector<double> &coord,
  //                     double sizeFactor, int nPoints);

  template <class MetricType>
  void classicalP1Interpolation(const double *xsi,
                                const MetricType &M0,
                                const MetricType &M1,
                                const MetricType &M2,
                                MetricType &result)
  {
    // Hardcoded P1 Lagrange basis functions for simplicity
    double phi[3] = {1. - xsi[0] - xsi[1], xsi[0], xsi[1]};
    result = M0 * phi[0] + M1 * phi[1] + M2 * phi[2];
  }

  template <class MetricType>
  void logEuclidianP1Interpolation(const double *xsi,
                                   const MetricType &M0,
                                   const MetricType &M1,
                                   const MetricType &M2,
                                   MetricType &result)
  {
    // Hardcoded P1 Lagrange basis functions for simplicity
    double phi[3] = {1. - xsi[0] - xsi[1], xsi[0], xsi[1]};
    // Interpolate log(M) then take exponential
    MetricType tmp = M0.log() * phi[0] + M1.log() * phi[1] + M2.log() * phi[2];
    result = tmp.exp();
  }

  template <class MetricType>
  void logEuclidianP2Interpolation(const double *xsi,
                                   const MetricType &M0,
                                   const MetricType &M1,
                                   const MetricType &M2,
                                   const MetricType &M3,
                                   const MetricType &M4,
                                   const MetricType &M5,
                                   MetricType &result)
  {
    // Hardcoded P2 Lagrange basis functions for simplicity
    double phi[6] = {(1. - xsi[0] - xsi[1]) * (1. - 2. * xsi[0] - 2. * xsi[1]),
                      xsi[0] * (2. * xsi[0] - 1.),
                      xsi[1] * (2. * xsi[1] - 1.),
                      4. * xsi[0] * (1. - xsi[0] - xsi[1]),
                      4. * xsi[0] * xsi[1],
                      4. * xsi[1] * (1. - xsi[0] - xsi[1])};
    // Interpolate log(M) then take exponential
    MetricType tmp = M0.log() * phi[0] + M1.log() * phi[1] + M2.log() * phi[2] + 
                     M3.log() * phi[3] + M4.log() * phi[4] + M5.log() * phi[5];
    result = tmp.exp();
  }

  void interpolationTest(const std::vector<size_t> &nodeTags, std::vector<double> &coord);

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
