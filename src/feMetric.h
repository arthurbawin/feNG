#ifndef _FEMETRIC_
#define _FEMETRIC_

#include "feRecovery.h"
#include "feNewRecovery.h"

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include "STensor3.h"

enum class Norm {
  // Target norm or seminorm to minimize with mesh adaptation
  Lp,
  Wmp
};

enum class GradationSpace {
  // In which space the graded metric field will be homoegeneous
  // Metric space preserves anisotropy
  // Physical space leads to more isotropic metrics
  Metric,
  Physical,
  Mixed,
  ExpMetric
};

enum class adaptationMethod {
  // Adapt with respect to an analytic Riemannian metric provided
  //
  ANALYTIC_METRIC,
  // Isotropic adaptation minimizing linear interpolation error
  // The metric is simply the most critical size from the ANISO_P1
  // computation times the identity matrix.
  ISO_P1,
  // Anisotropic adaptation minimizing linear interpolation error
  ANISO_P1,
  ANISO_P1_REFERENCE_SPACE,
  // Same but for any order on straight meshes
  ISO_PN,
  // Compute the analytic solution for the upper bound Q on the 3rd order derivatives
  ANISO_P2,
  ANISO_P2_REFERENCE_SPACE,
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
  CURVED_ISOLINES,
  //
  //
  CURVED_REFERENCE_SPACE,
  //
  //
  CURVED_GRAPH_SURFACE
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
  //
  double maxAnalyticEllipseDiameter = 1e2;
  // Mesh adaptation method
  adaptationMethod method = adaptationMethod::ANISO_P1;
  // Degree of the polynomial solution for the field driving
  // the mesh adaptation
  int polynomialDegree = 1;
  // Target number of vertices in adapted mesh
  // Used for anisotropic adaptation and curved adaptation
  // with the log-simplex method, ignored otherwise
  long int nTargetVertices = 100;
  // Target interpolation error in adapted mesh
  // Used for curved adaptation along iso/grad, ignored otherwise
  double eTargetError = 1e-2;
  // Limit variations in the metric field for smoother mesh
  bool enableGradation = false;
  // Maximum growth in the metric field if gradation is enabled
  // See F. Alauzet, Size gradation control of anisotropic meshes
  double gradation = 1.5;
  int gradationMaxIter = 100;
  GradationSpace gradationSpace = GradationSpace::Metric;
  // Smooth the directions field
  bool smoothDirectionField = true;
  // Scale the isocontour metric field as for the optimal metric for straight simplices
  bool scaleIsolineMetric = false;
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
  feVectorFunction *firstDerivatives = nullptr;
  feVectorFunction *secondDerivatives = nullptr;
  feVectorFunction *thirdDerivatives = nullptr;
  feVectorFunction *fourthDerivatives = nullptr;
  feVectorFunction *fifthDerivatives = nullptr;

  // Function to compute analytic metric instead of error-based metric
  feMetricTensorFunction *analyticMetric = nullptr;

  // Derivatives of M^(-1/2) (Mmud = M moins une demi (-: )
  feMetricTensorFunction *analytic_dMmud_dx = nullptr;
  feMetricTensorFunction *analytic_dMmud_dy = nullptr;

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

  struct {
    // Number of points for the discretization of a quarter of the unit circle
    bool enableCurvature = true;

    int maxIterMetrics = 10;
    double tolMetrics = 1e-10;
    bool printMetricConvergence = false;

    int maxIterDerivatives = 10;
    double tolDerivatives = 1e-10;

    double alpha = 0.;
    double alphaMax = 1.;

    // Number of iterations to go from alpha to alphaMax
    int numIterAlpha = 5;

    // If true: compute derivatives of M^(-1/2) by averaging adjacent FE gradients
    // If false: recover the gradient of each components of M^(-1/2) with Zhang and Naga's method
    bool averageMetricGradient = true;

    // Source term for verification
    bool addSourceTerm = false;
    feMetricTensorFunction *sourceTerm = nullptr;

    bool useAnalyticSolution = true;
    bool enableSmoothing = false;
    double toleranceDirectionsSmoothing = 0.02;
    bool enableGradation = false;
    double gradation = 2.;
    bool debugPlot = false;

    bool displayDerivatives = false;
    bool displayQtri = false;

    bool dampMetrics = true;
    double dampingMetricCoefficient = 0.5;

    bool dampDerivatives = true;
    double dampingCoefficient = 0.5;
  } curvedMetricReferenceSpaceOptions;

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

feStatus createMetricField(feMetric *&metricField,
                           feMetricOptions &options,
                           feMesh2DP1 *mesh,
                           std::vector<feNewRecovery*> recoveredFields = std::vector<feNewRecovery*>{});

class feMetric
{
public:
  feMetricOptions _options;
  int _nVerticesPerElmOnBackmesh;
  int _backmeshOrder;

protected:
  const feMesh2DP1 *_mesh;
  std::vector<feNewRecovery *> _recoveredFields;

  double _lambdaMin = 1e-14;
  double _lambdaMax = 1e22;

  // Node tags and coordinates of the Gmsh model
  // FIXME: Should eventually be removed to use only internal mesh structure
  std::vector<std::size_t> _nodeTags;
  std::vector<double> _coord;

  int _nVerticesPerElm;

  // =========================================================
  // New interface uses only those types:
  // Ultimately, the "good" type: mapping nodeTags to "Eigen-agnostic" MetricTensors
  std::map<int, MetricTensor> _metricTensorAtNodetags;
  std::map<int, MetricTensor> _logMetricTensorAtNodetags;
  std::vector<Eigen::Matrix2d> _logMetricTensorAtNodetags_eigen;
  std::map<int, SMetric3> _smetric3AtNodetags;

  // Stored in a vector, not in a map
  // With the internal structure node order, not Gmsh
  std::vector<MetricTensor> _metrics;
  std::vector<MetricTensor> _logMetrics;

  // M^(-1/2) and its derivatives
  std::map<int, MetricTensor> _mmud;
  std::map<int, MetricTensor> _dmmud_dy;
  std::map<int, MetricTensor> _dmmud_dx;
  // =========================================================

  // std::map<const Vertex *, Eigen::Matrix2d> _metricsOnGmshModelP1;
  // std::map<int, Eigen::Matrix2d> _metricsOnGmshModel_eigen;

  // Vertex to nodeTag map (and inverse)
  // nodeTag is the Gmsh node tag obtained from gmsh::model::mesh::getNodes
  // _nodeTag2sequentialTag maps the nodeTags to the internal tag, which is higher
  // than the total number of nodeTags when a P2 mesh is used.
  std::map<const Vertex *, int> _v2n;
  std::map<int, const Vertex *> _n2v;
  std::map<int, int> _nodeTag2sequentialTag;
  std::map<int, int> _sequentialTag2nodeTag;
  std::vector<int> _sequentialTag2nodeTagVec;

  // The tag of the gmsh view in which the metric field is stored (if using Gmsh)
  int _metricViewTag = -1;

  double _currentTime = 0.;

  bool _scaleMetricField = true;

  std::string _debugPrefix = "";

public:
  feMetric(const feMetricOptions &options,
           feMesh2DP1 *mesh,
           std::vector<feNewRecovery*> recoveredFields = std::vector<feNewRecovery*>{});
  ~feMetric()
  {
  #if defined(HAVE_GMSH)
    finalizeGmsh();
  #endif
  };

  void setOptions(feMetricOptions &other) { _options = other; }
  const feMetricOptions &getOptions() const { return _options; }
  void setGmshMetricModel(std::string metricModel) { _options.modelForMetric = metricModel; }
  // void setRecovery(feRecovery *other) { _recovery = other; }
  void setRecoveredFields(std::vector<feNewRecovery*> other) { _recoveredFields = other; }
  void setCurrentTime(const double time) { _currentTime = time; }
  void setMetricScaling(const bool flag) { _scaleMetricField = flag; }
  void setDebugPrefix(const std::string prefix) { _debugPrefix = prefix; }

  const std::vector<std::size_t> &getGmshNodeTags() const { return _nodeTags; }
  const std::vector<double> &getGmshNodeCoord() const { return _coord; }

  const std::map<int, MetricTensor> &getMetrics() const { return _metricTensorAtNodetags; }
  void getMetrics(std::map<int, MetricTensor> &other) const 
  {
    for(auto n : _nodeTags) {
      other[n].assignMatrixFrom(_metricTensorAtNodetags.at(n));
    }
  }
  void setMetrics(const std::map<int, MetricTensor> &other)
  {
    for(auto n : _nodeTags) {
      _metricTensorAtNodetags.at(n).assignMatrixFrom(other.at(n));
      _metrics[_nodeTag2sequentialTag[n]].assignMatrixFrom(other.at(n));
    }
  }
  void setMetricsToZero()
  {
    for(auto n : _nodeTags) {
      _metricTensorAtNodetags.at(n).setToZero();
      _metrics[_nodeTag2sequentialTag[n]].setToZero();
    }
  }
  void setMetricsToIdentity()
  {
    for(auto n : _nodeTags) {
      _metricTensorAtNodetags.at(n).setToIdentity();
      _metrics[_nodeTag2sequentialTag[n]].setToIdentity();
    }
  }

  // Apply other += alpha * metrics
  void addMetricsToOther(const double alpha, std::map<int, MetricTensor> &other)
  {
    for(const auto &tag : _nodeTags) {
      other.at(tag).incrementMatrix(alpha, _metricTensorAtNodetags.at(tag));
    }
  }

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
  #if defined(HAVE_GMSH)
  feStatus computeMetrics();
  #endif

  feStatus computeMetricsP1();
  feStatus computeMetricsHechtKuate();
  feStatus computeMetricsLogSimplex();
  feStatus computeMetricsExtremeSizesOnly();

  // ====================
  // New interface
  #if defined(HAVE_GMSH)
  feStatus createVertex2NodeMap(std::vector<std::size_t> &nodeTags, std::vector<double> &coord);
  void metricScalingFromGmshSubstitute(std::map<int, MetricTensor> &metrics,
                                       double exponentInIntegral,
                                       double exponentForDeterminant);
  void writeMetricField();
  #endif

  void computeLogMetrics();
  double computeIntegralOfDeterminant(const double exponentInIntegral);

  // Simply scale the metric field according to :
  //
  // M = scaling * M,
  //
  void scaleMetricsByConstant(const double scaling)
  {
    for(auto &pair : _metricTensorAtNodetags) { pair.second *= scaling; }
  }
  // Scale the metric field according to :
  //
  // M = C * det(M)^a * M,
  //
  // with C = globalScalingFactor and a = exponentForDeterminant
  void scaleMetricsByDeterminant(const double globalScalingFactor,
                                 const double exponentForDeterminant);

  feStatus setAnalyticMetric(const std::vector<std::size_t> &nodeTags,
                             const std::vector<double> &coord);
  feStatus computeMetricsP1(std::vector<std::size_t> &nodeTags, std::vector<double> &coord,
                            bool isotropic = false);
  feStatus computeMetricsP2(std::vector<std::size_t> &nodeTags, std::vector<double> &coord);
  feStatus computeMetricsP1_referenceSpace(std::vector<std::size_t> &nodeTags,
                                           std::vector<double> &coord, int setToIdentity = false);
  feStatus computeMetricsP2_referenceSpace(std::vector<std::size_t> &nodeTags,
                                           std::vector<double> &coord, int setToIdentity = false);
  feStatus computeMetricsP2_forGraphSurface(std::vector<std::size_t> &nodeTags,
                                            std::vector<double> &coord);
  feStatus computeMetricsPn(std::vector<std::size_t> &nodeTags, std::vector<double> &coord,
                            bool isotropic = false);
  feStatus computeMetricsGoalOrientedP1(std::vector<std::size_t> &nodeTags,
                                        std::vector<double> &coord);
  feStatus computeMetricsCurvedLogSimplex(std::vector<std::size_t> &nodeTags,
                                          std::vector<double> &coord,
                                          bool useInducedDirections = false);

  feStatus computeMetricsCurvedReferenceSpace(std::vector<std::size_t> &nodeTags,
                                              std::vector<double> &coord);
  feStatus computeMetricsCurvedReferenceSpace2(std::vector<std::size_t> &nodeTags,
                                               std::vector<double> &coord);

  void getDerivativesAtVertex(const int derivativesOrder, const double *x, const size_t vertex,
                              feNewRecovery *recoveredField, double *du, double *d2u = nullptr,
                              double *d3u = nullptr, double *d4u = nullptr, double *d5u = nullptr);

  void computeSizeField(const double pos[2], const int vertex, const double directionV1[2],
                        double &hGrad, double &hIso, const double *du, const double *d2u,
                        const double *d3u);

  feStatus computeMetricsCurvedIsolines(const std::vector<std::size_t> &nodeTags,
                                        const std::vector<double> &coord);

  void applyGradation();
  feStatus newGradation(const std::vector<std::size_t> &nodeTags, const std::vector<double> &coord,
                        const double gradation, std::map<int, MetricTensor> &metricsAtNodeTags);

  MetricTensor &getMetricAtSequentialTag(int tag)
  {
    return _metricTensorAtNodetags[_sequentialTag2nodeTag[tag]];
  };

  #if defined(HAVE_GMSH)
  void averageMetricGradientAtVertices(const std::vector<std::size_t> &nodeTags,
                                       const std::vector<double> &coord,
                                       const std::map<int, MetricTensor> &inputMetrics,
                                       std::map<int, MetricTensor> &dMdx,
                                       std::map<int, MetricTensor> &dMdy);
  void computeErrorOnMetricDerivatives(double &errorOnMmud_dx, double &errorOnMmud_dy);
  void computeContinuousErrorModel(const bool withAnalyticMetric, bool computeError, bool plotError,
                                   const std::string &errorFileName, double &errorStraight,
                                   double &errorCurved, double &complexity);
  #endif
  // ====================

  void setMetricViewTag(int tag) { _metricViewTag = tag; }
  int getMetricViewTag() { return _metricViewTag; }

  // Intersecter avec un autre feMetric
  void writeSizeFieldSol2D(std::string solFileName);
  void writeSizeFieldSol3D(std::string solFileName);
  void writeSizeFieldGmsh(std::string meshName, std::string metricMeshName);

  void interpolationTest(const std::vector<size_t> &nodeTags, std::vector<double> &coord);

  ///////////////////////////////////////////////////////////
  // Metric interpolation on background mesh
  void interpolateMetricP1(const double *x, Eigen::Matrix2d &M, Eigen::Matrix2d &dMdx,
                           Eigen::Matrix2d &dMdy);
  void interpolateMetricP2(const double *x, Eigen::Matrix2d &M);
  void interpolateMetricP2Explicit(const double *x, Eigen::Matrix2d &M);
  void interpolateMetricP2(const double *x, Eigen::Matrix2d &M, Eigen::Matrix2d &dMdx,
                           Eigen::Matrix2d &dMdy);
  void interpolateMetricP2Log(const double *x, Eigen::Matrix2d &M, Eigen::Matrix2d &dMdx,
                              Eigen::Matrix2d &dMdy, Eigen::Matrix2d &L, Eigen::Matrix2d &dLdx,
                              Eigen::Matrix2d &dLdy, double &l1, double &dl1dx, double &dl1dy,
                              double &l2, double &dl2dx, double &dl2dy, Eigen::Vector2d &u1,
                              Eigen::Vector2d &du1dx, Eigen::Vector2d &du1dy, Eigen::Vector2d &u2,
                              Eigen::Vector2d &du2dx, Eigen::Vector2d &du2dy);
  void gradLogEuclidianP1Interpolation(const double *xsi, const int element,
                                       const MetricTensor &logM0, const MetricTensor &logM1,
                                       const MetricTensor &logM2, Eigen::Matrix2d &dMdx,
                                       Eigen::Matrix2d &dMdy);
  void gradLogEuclidianP2Interpolation(
    const double xsi[2], const int element, const Eigen::Matrix2d &logM0,
    const Eigen::Matrix2d &logM1, const Eigen::Matrix2d &logM2, const Eigen::Matrix2d &logM3,
    const Eigen::Matrix2d &logM4, const Eigen::Matrix2d &logM5, Eigen::Matrix2d &dMdx,
    Eigen::Matrix2d &dMdy, Eigen::Matrix2d &LRES, Eigen::Matrix2d &DLDXRES,
    Eigen::Matrix2d &DLDYRES, double &l1res, double &dl1dxres, double &dl1dyres, double &l2res,
    double &dl2dxres, double &dl2dyres, Eigen::Vector2d &u1res, Eigen::Vector2d &du1dxres,
    Eigen::Vector2d &du1dyres, Eigen::Vector2d &u2res, Eigen::Vector2d &du2dxres,
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

SMetric3 intersectionReductionSimultaneeExplicite(const SMetric3 &m1, const SMetric3 &m2);

template <class MetricType>
void drawSingleEllipse(FILE *file, const double *x, MetricType &M, double sizeFactor = 1.,
                       int nPoints = 30);
template <class MetricType>
void drawEllipsoids(const std::string &posFile, std::map<int, MetricType> &metrics,
                    const std::vector<std::size_t> &nodeTags, const std::vector<double> &coord,
                    double sizeFactor, int nPoints);

// Non-class wrapper
void interpolateMetricP1Callback(void *metricPtr, const double *x, Eigen::Matrix2d &M,
                                 Eigen::Matrix2d &dMdx, Eigen::Matrix2d &dMdy);
void interpolateMetricP2CallbackWithoutDerivatives(void *metricPtr, const double *x,
                                                   Eigen::Matrix2d &M);
void interpolateMetricP2CallbackWithoutDerivativesExplicit(void *metricPtr, const double *x,
                                                           Eigen::Matrix2d &M);
void interpolateMetricP2CallbackWithDerivatives(void *metricPtr, const double *x,
                                                Eigen::Matrix2d &M, Eigen::Matrix2d &dMdx,
                                                Eigen::Matrix2d &dMdy);
void interpolateMetricP2CallbackLog(void *metricPtr, const double *x, Eigen::Matrix2d &M,
                                    Eigen::Matrix2d &dMdx, Eigen::Matrix2d &dMdy,
                                    Eigen::Matrix2d &L, Eigen::Matrix2d &dLdx,
                                    Eigen::Matrix2d &dLdy, double &l1, double &dl1dx, double &dl1dy,
                                    double &l2, double &dl2dx, double &dl2dy, Eigen::Vector2d &u1,
                                    Eigen::Vector2d &du1dx, Eigen::Vector2d &du1dy,
                                    Eigen::Vector2d &u2, Eigen::Vector2d &du2dx,
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
