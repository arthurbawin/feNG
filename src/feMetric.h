#ifndef _FEMETRIC_
#define _FEMETRIC_

#include "feRecovery.h"
#include "feNewRecovery.h"

#include "../contrib/Eigen/Eigen"
#include "../contrib/unsupported/Eigen/MatrixFunctions"

#include "STensor3.h"

enum class Norm
{ 
  // Target norm or seminorm to minimize with mesh adaptation
  Lp, Wmp
};

enum class adaptationMethod
{ 
  // Isotropic adaptation minimizing linear interpolation error
  // The metric is simply the most critical size from the ANISO_P1
  // computation times the identity matrix.
  ISO_P1,
  // Anisotropic adaptation minimizing linear interpolation error
  ANISO_P1,
  // Same but for any order on straight meshes
  ISO_PN,
  // Compute the analytic solution for the upper bound Q on the 3rd order derivatives
  ANISO_P2,
  // Compute an approximate the upper bound Q on the n+1th order derivatives
  // using the logsimplex method from Coulaud & Loseille (2016)
  ANISO_PN,
  // Goal-oriented anisotropic adaptation for linear interpolant
  // Requires the computation of the adjoint solution p
  GOAL_ORIENTED_ANISO_P1,
  // Metric for curved adaptation minimizing P2 interpolation error on P2 meshes.
  // Extension of ANISO_PN to P2 meshes with the log-simplex method, may or may not be correct...
  CURVED_LS,
  // Metric for curved adaptation minimizing P2 interpolation error on P2 meshes.
  // Principal sizes are computed from the CURVED_LS metric
  // Principal directions are those of the induced metric (gradient/isolines)
  CURVED_LS_INDUCED_DIRECTIONS,
  // Metric for curved adaptation minimizing P2 interpolation error on P2 meshes.
  // Principal directions (eigenvectors) are the recovered gradient/isolines of the solution
  // Principal sizes are computed from a Taylor expansion on curves approximating grad/iso
  CURVED_EXTREME_SIZES
};

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
  // Target norm or seminorm to minimize - choice between Lp or Wmp
  Norm targetNorm = Norm::Lp;
  // Target Lp norm in which the interpolation error is minimized
  // See e.g. Alauzet & Loseille, Multi-Dimensional Continuous Metric
  // for Mesh Adaptation
  double LpNorm = 2.;
  double WmpNorm_mParam = 1.;
  double WmpNorm_pParam = 2.;

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
  double plotSizeFactor = 100.;

  // Use analytic derivatives? For debug only
  bool useAnalyticDerivatives = false;
  feVectorFunction *analyticDerivatives = nullptr;

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
    double tol = 5e-2;
    bool uniformErrorCurve = false;
    int numLoopsUniformErrorCurve = 1;
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
  int _backmeshOrder;

protected:
  feRecovery *_recovery;
  std::vector<feNewRecovery*> _recoveredFields;

  double _lambdaMin = 1e-14;
  double _lambdaMax = 1e22;

  std::map<int, SMetric3> _metrics;
  std::map<int, SMetric3> _metricsOnGmshModel;

  std::map<int, SMetric3> _metricsP1;
  std::map<int, Eigen::Matrix2d> _metricsP1AtNodeTags;

  // =========================================================
  // New interface uses only those types:
  // Ultimately, the "good" type: mapping nodeTags to "Eigen-agnostic" MetricTensors
  std::map<int, MetricTensor> _metricTensorAtNodetags;
  std::map<int, MetricTensor> _logMetricTensorAtNodetags;
  std::vector<Eigen::Matrix2d> _logMetricTensorAtNodetags_eigen;
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
  feMetric(std::vector<feNewRecovery*> &recoveredFields, feMetricOptions metricOptions);
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
  feStatus computeMetricsP1(std::vector<std::size_t> &nodeTags, std::vector<double> &coord, bool isotropic = false);
  void computeAnalyticMetricP2ForLpNorm(std::vector<double> &errorCoeff, MetricTensor &Q, const double maxDiameter);
  feStatus computeMetricsP2(std::vector<std::size_t> &nodeTags, std::vector<double> &coord);
  feStatus computeMetricsPn(std::vector<std::size_t> &nodeTags, std::vector<double> &coord, bool isotropic = false);
  feStatus computeMetricsGoalOrientedP1(std::vector<std::size_t> &nodeTags, std::vector<double> &coord);
  feStatus computeMetricsCurvedLogSimplex(std::vector<std::size_t> &nodeTags, std::vector<double> &coord, bool useInducedDirections = false);
  void computeDirectionFieldFromGradient(const int vertex, double directionV1[2], const double tol);
  double dttt(const int vertex, double directionV1[2], const int direction);
  feStatus computeMetricsExtremeSizesOnly(std::vector<std::size_t> &nodeTags, std::vector<double> &coord);

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
                                   const MetricType &logM0,
                                   const MetricType &logM1,
                                   const MetricType &logM2,
                                   MetricType &result)
  {
    // Hardcoded P1 Lagrange basis functions for simplicity
    double phi[3] = {1. - xsi[0] - xsi[1], xsi[0], xsi[1]};
    // Interpolate log(M) then take exponential
    result = logM0 * phi[0] + logM1 * phi[1] + logM2 * phi[2];
    result = result.exp();
  }

  template <class MetricType>
  void logEuclidianP2Interpolation(const double *xsi,
                                   const MetricType &logM0,
                                   const MetricType &logM1,
                                   const MetricType &logM2,
                                   const MetricType &logM3,
                                   const MetricType &logM4,
                                   const MetricType &logM5,
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
    result = (logM0 * phi[0] + logM1 * phi[1] + logM2 * phi[2] + 
              logM3 * phi[3] + logM4 * phi[4] + logM5 * phi[5]).exp();
  }

  void interpolationTest(const std::vector<size_t> &nodeTags, std::vector<double> &coord);


  ///////////////////////////////////////////////////////////
  // Metric interpolation on background mesh
  void interpolateMetricP1(const double *x, Eigen::Matrix2d &M,
                           Eigen::Matrix2d &dMdx, Eigen::Matrix2d &dMdy);
  void interpolateMetricP2(const double *x, Eigen::Matrix2d &M,
                           Eigen::Matrix2d &dMdx, Eigen::Matrix2d &dMdy);
  void interpolateMetricP2Log(const double *x, Eigen::Matrix2d &M,
                                   Eigen::Matrix2d &dMdx, Eigen::Matrix2d &dMdy,
                                   Eigen::Matrix2d &L,
                                   Eigen::Matrix2d &dLdx, Eigen::Matrix2d &dLdy,
                                   double &l1,
                                   double &dl1dx,
                                   double &dl1dy,
                                   double &l2,
                                   double &dl2dx,
                                   double &dl2dy,
                                   Eigen::Vector2d &u1,
                                   Eigen::Vector2d &du1dx,
                                   Eigen::Vector2d &du1dy,
                                   Eigen::Vector2d &u2,
                                   Eigen::Vector2d &du2dx,
                                   Eigen::Vector2d &du2dy);
  void gradLogEuclidianP1Interpolation(const double *xsi,
                                       const int element,
                                       const MetricTensor &logM0,
                                       const MetricTensor &logM1,
                                       const MetricTensor &logM2,
                                       Eigen::Matrix2d &dMdx,
                                       Eigen::Matrix2d &dMdy);
  void gradLogEuclidianP2Interpolation(const double xsi[2],
                                       const int element,
                                       const Eigen::Matrix2d &logM0,
                                       const Eigen::Matrix2d &logM1,
                                       const Eigen::Matrix2d &logM2,
                                       const Eigen::Matrix2d &logM3,
                                       const Eigen::Matrix2d &logM4,
                                       const Eigen::Matrix2d &logM5,
                                       Eigen::Matrix2d &dMdx,
                                       Eigen::Matrix2d &dMdy,
                                       Eigen::Matrix2d &LRES,
                                       Eigen::Matrix2d &DLDXRES,
                                       Eigen::Matrix2d &DLDYRES,
                                       double &l1res,
                                       double &dl1dxres,
                                       double &dl1dyres,
                                       double &l2res,
                                       double &dl2dxres,
                                       double &dl2dyres,
                                       Eigen::Vector2d &u1res,
                                       Eigen::Vector2d &du1dxres,
                                       Eigen::Vector2d &du1dyres,
                                       Eigen::Vector2d &u2res,
                                       Eigen::Vector2d &du2dxres,
                                       Eigen::Vector2d &du2dyres);

  // With gradient wrt to physical coordinates x,y
  // void interpolateMetricP1WithDerivatives(const double *x, Eigen::Matrix2d &M,
  //                                         Eigen::Matrix2d &dMdx, Eigen::Matrix2d &dMdy);

//   // With gradient of M wrt 2d position alpha (a1,a2)
//   void interpolateMetricP1(const double *x, Eigen::Matrix2d &M, Eigen::Matrix2d &sumduda1M,
//                            Eigen::Matrix2d &sumduda2M);
//   void interpolateMetricAndDerivativeOnP2Edge(const double t, const Eigen::Matrix2d &M11,
//                                               const Eigen::Matrix2d &M20,
//                                               const Eigen::Matrix2d &M02,
//                                               const Eigen::Matrix2d &sumduda1M,
//                                               const Eigen::Matrix2d &sumduda2M, Eigen::Matrix2d
//                                               &M, Eigen::Matrix2d &dMda1, Eigen::Matrix2d
//                                               &dMda2);

//   // With derivative of M wrt 1D alpha (dMda)
//   void interpolateMetricP1(const double *x, const double *gammaOrth, Eigen::Matrix2d &M,
//                            Eigen::Matrix2d &sumdudaM);
//   void interpolateMetricAndDerivativeOnP2Edge(double t, const Eigen::Matrix2d &M11,
//                                               const Eigen::Matrix2d &M20,
//                                               const Eigen::Matrix2d &M02,
//                                               const Eigen::Matrix2d &sumdudaM, Eigen::Matrix2d
//                                               &M, Eigen::Matrix2d &dMda);
};

// Non-class wrapper
void interpolateMetricP1Callback(void *metricPtr, const double *x, Eigen::Matrix2d &M,
                           Eigen::Matrix2d &dMdx, Eigen::Matrix2d &dMdy);
void interpolateMetricP2Callback(void *metricPtr, const double *x, Eigen::Matrix2d &M,
                           Eigen::Matrix2d &dMdx, Eigen::Matrix2d &dMdy);
void interpolateMetricP2CallbackLog(void *metricPtr, const double *x, Eigen::Matrix2d &M,
                           Eigen::Matrix2d &dMdx, Eigen::Matrix2d &dMdy,
                           Eigen::Matrix2d &L,
                           Eigen::Matrix2d &dLdx, Eigen::Matrix2d &dLdy,
                           double &l1,
                           double &dl1dx,
                           double &dl1dy,
                           double &l2,
                           double &dl2dx,
                           double &dl2dy,
                           Eigen::Vector2d &u1,
                           Eigen::Vector2d &du1dx,
                           Eigen::Vector2d &du1dy,
                           Eigen::Vector2d &u2,
                           Eigen::Vector2d &du2dx,
                           Eigen::Vector2d &du2dy);

// void interpolateMetricP1WithDerivativesWrapper(void *metric, const double *x, Eigen::Matrix2d &M,
//                                                Eigen::Matrix2d &dMdx, Eigen::Matrix2d &dMdy);

// void interpolateMetricP1Wrapper(void *metric, const double *x, Eigen::Matrix2d &M,
//                                 Eigen::Matrix2d &sumduda1M, Eigen::Matrix2d &sumduda2M);

// void interpolateMetricP1Wrapper1D(void *metric, const double *x, const double *gammaOrth,
//                                   Eigen::Matrix2d &M, Eigen::Matrix2d &sumdudaM);

// void interpolateMetricAndDerivativeOnP2EdgeWrapper(
//   void *metric, const double t, const Eigen::Matrix2d &M11, const Eigen::Matrix2d &M20,
//   const Eigen::Matrix2d &M02, const Eigen::Matrix2d &sumduda1M, const Eigen::Matrix2d &sumduda2M,
//   Eigen::Matrix2d &M, Eigen::Matrix2d &dMda1, Eigen::Matrix2d &dMda2);

// void interpolateMetricAndDerivativeOnP2EdgeWrapper1D(void *metric, const double t,
//                                                      const Eigen::Matrix2d &M11,
//                                                      const Eigen::Matrix2d &M20,
//                                                      const Eigen::Matrix2d &M02,
//                                                      const Eigen::Matrix2d &sumdudaM,
//                                                      Eigen::Matrix2d &M, Eigen::Matrix2d &dMda);

#endif
