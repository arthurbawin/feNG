
#include "feMinimizeInterpolationError.h"

feNewRecovery *activeRecovery;
feSpace *activeIntSpace;
feNumber *activeNumbering;
feSolution *activeSolution;
feFunction *activeExactSolution;
feVectorFunction *activeExactSolutionGradient;

feMesh *activeMesh;
feCncGeo *activeConnectivity;
feNorm *activeNorm;

// #define ERROR_FROM_ESTIMATOR
// #define GRAD_FROM_ESTIMATOR

#define USE_AVERAGED_EVALUATIONS false

static double INITIAL_COORD[3];
static double LARGE_VALUE_WHEN_FAILED = std::numeric_limits<double>::infinity();;
static std::vector<double> ELEM_ERROR_SQUARED;
static double INITIAL_ERROR_SQUARED;
static std::vector<double> GEO_COORD(18, 0.);

#if defined(HAVE_GMSH)
feStatus moveVertexAndRecomputeJacobians(const std::vector<int> &whichElements,
                                         const int iVertexToMove_globalTag,
                                         const double *modifiedCoord)
{
  activeMesh->setVertexCoord(iVertexToMove_globalTag, modifiedCoord);
  bool failed = false;
  for(size_t i = 0; i < whichElements.size(); ++i) {
    feStatus s = activeConnectivity->recomputeElementJacobian(whichElements[i]);
    failed = (s != FE_STATUS_OK);
  }
  if(failed) return FE_STATUS_FAILED;
  return FE_STATUS_OK;
}

void computeInterpolationErrorOnEachElement()
{
  // Store the initial interpolation error on each element
  ELEM_ERROR_SQUARED.resize(activeConnectivity->getNumElements());
  INITIAL_ERROR_SQUARED = 0.;
  for(size_t i = 0; i < activeConnectivity->getNumElements(); ++i) {
#if defined(ERROR_FROM_ESTIMATOR)
    ELEM_ERROR_SQUARED[i] =
      activeNorm->computeSquaredErrorFromEstimatorOnElement(i, USE_AVERAGED_EVALUATIONS);
#else
    ELEM_ERROR_SQUARED[i] = activeNorm->computeSquaredErrorOnElement(i);
#endif
    INITIAL_ERROR_SQUARED += ELEM_ERROR_SQUARED[i];
  }
}

void setNewError(const int nElements, const std::vector<int> &whichElements)
{
  for(size_t i = 0; i < nElements; ++i) {
    INITIAL_ERROR_SQUARED -= ELEM_ERROR_SQUARED[whichElements[i]];
  }
  for(size_t i = 0; i < nElements; ++i) {
#if defined(ERROR_FROM_ESTIMATOR)
    ELEM_ERROR_SQUARED[whichElements[i]] = activeNorm->computeSquaredErrorFromEstimatorOnElement(
      whichElements[i], USE_AVERAGED_EVALUATIONS);
#else
    ELEM_ERROR_SQUARED[whichElements[i]] =
      activeNorm->computeSquaredErrorOnElement(whichElements[i]);
#endif
  }
  for(size_t i = 0; i < nElements; ++i) {
    INITIAL_ERROR_SQUARED += ELEM_ERROR_SQUARED[whichElements[i]];
  }
}

double getNewError(const int nElements, const std::vector<int> &whichElements)
{
  double res = INITIAL_ERROR_SQUARED;
  for(int i = 0; i < nElements; ++i) {
    res -= ELEM_ERROR_SQUARED[whichElements[i]];
  }
  for(int i = 0; i < nElements; ++i) {
#if defined(ERROR_FROM_ESTIMATOR)
    res += activeNorm->computeSquaredErrorFromEstimatorOnElement(whichElements[i],
                                                                 USE_AVERAGED_EVALUATIONS);
#else
    res += activeNorm->computeSquaredErrorOnElement(whichElements[i]);
#endif
  }
  return sqrt(res);
}

static std::vector<feInt> ADR(6);
static std::vector<double> POS(3);

void applyCurvatureToFeMesh(const edgeAndVertexData &data, const double *modifiedCoord)
{
  // Get interpolated solution at initial vertex
  activeIntSpace->initializeAddressingVector(data.whichElements[0], ADR);
  int vertexDOF = ADR[data.iVertexToMove_localTag[0]];

  feStatus s = moveVertexAndRecomputeJacobians(data.whichElements, data.iVertexToMove_globalTag,
                                               modifiedCoord);

  // Compute solution at new location and assign it
  POS[0] = modifiedCoord[0];
  POS[1] = modifiedCoord[1];
  POS[2] = 0.;

#if defined(ERROR_FROM_ESTIMATOR)
  // Option 1: New solution will be the evaluation of the current recovery
  double newUh =
    activeRecovery->evaluateRecovery(PPR::RECOVERY, 0, POS.data(), USE_AVERAGED_EVALUATIONS);

// Option 2: Compute a new recovery at the moved vertex (expensive?)
// activeRecovery->recomputeLeastSquareMatrix(data.iVertexToMove_globalTag);
// activeRecovery->recomputeRHSAndSolve(data.iVertexToMove_globalTag);
// activeRecovery->recomputeFirstDerivative(data.iVertexToMove_globalTag);
// double newUh = activeRecovery->evaluateRecoveryAtVertex(PPR::RECOVERY, 0,
// data.iVertexToMove_globalTag);
#else
  double newUh = activeExactSolution->eval(0., POS);
#endif

  activeSolution->setSolAtDOF(vertexDOF, newUh);

  // Update error vector and total error (substract old error and add new error)
  setNewError(data.nElements, data.whichElements);

  if(s != FE_STATUS_OK) {
    feInfo("Modified vertex caused at least one negative jacobian :/");
    exit(-1);
  }
}

void feMesh2DP1::drawConnectivityToPOSfile(const std::string &cncName, const std::string &fileName)
{
  feCncGeo *cnc = this->getCncGeoByName(cncName);
  if(cnc) {
    FILE *myfile = fopen(fileName.data(), "w");
    fprintf(myfile, "View\"%s\"{\n", cncName.data());
    std::vector<double> geoCoord(cnc->getNumVerticesPerElem() * 3, 0.);
    for(int iElm = 0; iElm < cnc->getNumElements(); ++iElm) {
      this->getCoord(cnc, iElm, geoCoord);
      cnc->writeElementToPOS(myfile, geoCoord, 1.);
    }
    fprintf(myfile, "};\n");
    fclose(myfile);
  }
}

static void getP2BernsteinBasis(const double uvw[3], double basis[6])
{
  double u = uvw[0];
  double v = uvw[1];
  double w = uvw[2];
  basis[0] = u * u;
  basis[1] = 2. * u * v;
  basis[2] = v * v;
  basis[3] = 2. * v * w;
  basis[4] = w * w;
  basis[5] = 2. * w * u;
}

// Lagrange points are given as {L200, L020, L002, L110, L011, L101}
static void getBezierControlPoints(const SPoint2 L[6], SPoint2 P[6])
{
  P[0] = L[0]; // P_200
  P[1] = L[1]; // P_020
  P[2] = L[2]; // P_002
  P[3] = (L[0] * (-1.) + L[3] * 4. - L[1]) * 0.5; // P_110
  P[4] = (L[1] * (-1.) + L[4] * 4. - L[2]) * 0.5; // P_011
  P[5] = (L[2] * (-1.) + L[5] * 4. - L[0]) * 0.5; // P_101

  // for(int i = 0; i < 6; ++i){
  //   feInfo("Lagrange point %d : (%f, %f) - Bézier point : (%f, %f)",
  //     i, L[i].x(), L[i].y(), P[i].x(), P[i].y());
  // }
}

// Compute determinant of the 2x2 matrix |P0 P1|,
// points are column vectors.
inline double computeDeterminant(const SPoint2 &P0, const SPoint2 &P1)
{
  return (P0[0] * P1[1]) - (P1[0] * P0[1]);
}

static double controlCoefficientN200(const SPoint2 &P200, const SPoint2 &P011, const SPoint2 &P020,
                                     const SPoint2 &P101, const SPoint2 &P002, const SPoint2 &P110)
{
  return 4. * computeDeterminant(P200 - P110, P200 - P101);
}

static double controlCoefficientN110(const SPoint2 &P200, const SPoint2 &P011, const SPoint2 &P020,
                                     const SPoint2 &P101, const SPoint2 &P002, const SPoint2 &P110)
{
  return 2. * computeDeterminant(P200 - P101, P020 - P011) +
         2. * computeDeterminant(P110 - P011, P110 - P101);
}

// Bezier points are given as {P200, P020, P002, P110, P011, P101}
static void getJacobianControlCoefficients(const SPoint2 P[6], double N[6])
{
  const SPoint2 &P200 = P[0];
  const SPoint2 &P020 = P[1];
  const SPoint2 &P002 = P[2];
  const SPoint2 &P110 = P[3];
  const SPoint2 &P011 = P[4];
  const SPoint2 &P101 = P[5];
  // N[0] = 4. * computeDeterminant(P200 - P110, P200 - P101); // N200
  // N[1] = 4. * computeDeterminant(P020 - P110, P011 - P110); // N020
  // N[2] = 4. * computeDeterminant(P002 - P011, P101 - P011); // N002
  // N[3] = 2. * computeDeterminant(P200 - P101, P020 - P011)
  //      + 2. * computeDeterminant(P110 - P011, P110 - P101); // N110 
  // N[4] = 2. * computeDeterminant(P020 - P110, P002 - P101)
  //      + 2. * computeDeterminant(P011 - P101, P011 - P110); // N011
  // N[5] = 2. * computeDeterminant(P002 - P011, P200 - P110)
  //      + 2. * computeDeterminant(P101 - P110, P101 - P011); // N101

  // Or using permutations:
  N[0] = controlCoefficientN200(P200, P011, P020, P101, P002, P110);
  N[1] = controlCoefficientN200(P020, P101, P002, P110, P200, P011);
  N[2] = controlCoefficientN200(P002, P110, P200, P011, P020, P101);
  N[3] = controlCoefficientN110(P200, P011, P020, P101, P002, P110);
  N[4] = controlCoefficientN110(P020, P101, P002, P110, P200, P011);
  N[5] = controlCoefficientN110(P002, P110, P200, P011, P020, P101);
}

inline SPoint2 orthogonal(const SPoint2 &P) { return SPoint2(P[1], -P[0]); }

// The three non-trivial gradients with respect to vertex control point P200
static SPoint2 gradControlCoefficientP200N200(const SPoint2 &P200, const SPoint2 &P011,
                                              const SPoint2 &P020, const SPoint2 &P101,
                                              const SPoint2 &P002, const SPoint2 &P110)
{
  return orthogonal(P110 - P101) * 4.;
}
static SPoint2 gradControlCoefficientP200N110(const SPoint2 &P200, const SPoint2 &P011,
                                              const SPoint2 &P020, const SPoint2 &P101,
                                              const SPoint2 &P002, const SPoint2 &P110)
{
  return orthogonal(P020 - P011) * 2.;
}
static SPoint2 gradControlCoefficientP200N101(const SPoint2 &P200, const SPoint2 &P011,
                                              const SPoint2 &P020, const SPoint2 &P101,
                                              const SPoint2 &P002, const SPoint2 &P110)
{
  return orthogonal(P011 - P002) * 2.;
}

// The five non-trivial gradients with respect to edge control point P110
static SPoint2 gradControlCoefficientP110N200(const SPoint2 &P200, const SPoint2 &P011,
                                              const SPoint2 &P020, const SPoint2 &P101,
                                              const SPoint2 &P002, const SPoint2 &P110)
{
  return orthogonal(P101 - P200) * 4.;
}
static SPoint2 gradControlCoefficientP110N110(const SPoint2 &P200, const SPoint2 &P011,
                                              const SPoint2 &P020, const SPoint2 &P101,
                                              const SPoint2 &P002, const SPoint2 &P110)
{
  return orthogonal(P011 - P101) * 2.;
}
static SPoint2 gradControlCoefficientP110N011(const SPoint2 &P200, const SPoint2 &P011,
                                              const SPoint2 &P020, const SPoint2 &P101,
                                              const SPoint2 &P002, const SPoint2 &P110)
{
  return orthogonal(P011 - P002) * 2.;
}
static SPoint2 gradControlCoefficientP110N020(const SPoint2 &P200, const SPoint2 &P011,
                                              const SPoint2 &P020, const SPoint2 &P101,
                                              const SPoint2 &P002, const SPoint2 &P110)
{
  return orthogonal(P011 - P020) * (-4.);
}
static SPoint2 gradControlCoefficientP110N101(const SPoint2 &P200, const SPoint2 &P011,
                                              const SPoint2 &P020, const SPoint2 &P101,
                                              const SPoint2 &P002, const SPoint2 &P110)
{
  return orthogonal(P101 - P002) * (-2.);
}

static void getGradientJacobianControlCoefficients(const SPoint2 P[6], SPoint2 gradN[6][6])
{
  const SPoint2 &P200 = P[0];
  const SPoint2 &P020 = P[1];
  const SPoint2 &P002 = P[2];
  const SPoint2 &P110 = P[3];
  const SPoint2 &P011 = P[4];
  const SPoint2 &P101 = P[5];

  // grad_P200 of N200, N020, N002, N110, N011, N101
  gradN[0][0] = gradControlCoefficientP200N200(P200, P011, P020, P101, P002, P110); // N200
  gradN[0][1] = 0.; // N020
  gradN[0][2] = 0.; // N002
  gradN[0][3] = gradControlCoefficientP200N110(P200, P011, P020, P101, P002, P110); // N110
  gradN[0][4] = 0.; // N011
  gradN[0][5] = gradControlCoefficientP200N101(P200, P011, P020, P101, P002, P110); // N101

  // grad_P020 of N200, N020, N002, N110, N011, N101
  gradN[1][0] = 0.; // N200
  gradN[1][1] = gradControlCoefficientP200N200(P020, P101, P002, P110, P200, P011); // N020
  gradN[1][2] = 0.; // N002
  gradN[1][3] = gradControlCoefficientP200N101(P020, P101, P002, P110, P200, P011); // N110
  gradN[1][4] = gradControlCoefficientP200N110(P020, P101, P002, P110, P200, P011); // N011
  gradN[1][5] = 0.; // N101

  // grad_P002 of N200, N020, N002, N110, N011, N101
  gradN[2][0] = 0.; // N200
  gradN[2][1] = 0.; // N020
  gradN[2][2] = gradControlCoefficientP200N200(P002, P110, P200, P011, P020, P101); // N002
  gradN[2][3] = 0.; // N110
  gradN[2][4] = gradControlCoefficientP200N101(P002, P110, P200, P011, P020, P101); // N011
  gradN[2][5] = gradControlCoefficientP200N110(P002, P110, P200, P011, P020, P101); // N101

  // grad_P110 of N200, N020, N002, N110, N011, N101
  gradN[3][0] = gradControlCoefficientP110N200(P200, P011, P020, P101, P002, P110);
  gradN[3][1] = gradControlCoefficientP110N020(P200, P011, P020, P101, P002, P110);
  gradN[3][2] = 0.;
  gradN[3][3] = gradControlCoefficientP110N110(P200, P011, P020, P101, P002, P110);
  gradN[3][4] = gradControlCoefficientP110N011(P200, P011, P020, P101, P002, P110);
  gradN[3][5] = gradControlCoefficientP110N101(P200, P011, P020, P101, P002, P110);

  // grad_P011 of N200, N020, N002, N110, N011, N101
  gradN[4][0] = 0.;
  gradN[4][1] = gradControlCoefficientP110N200(P020, P101, P002, P110, P200, P011);
  gradN[4][2] = gradControlCoefficientP110N020(P020, P101, P002, P110, P200, P011);
  gradN[4][3] = gradControlCoefficientP110N101(P020, P101, P002, P110, P200, P011);
  gradN[4][4] = gradControlCoefficientP110N110(P020, P101, P002, P110, P200, P011);
  gradN[4][5] = gradControlCoefficientP110N011(P020, P101, P002, P110, P200, P011);

  // grad_P101 of N200, N020, N002, N110, N011, N101
  gradN[5][0] = gradControlCoefficientP110N020(P002, P110, P200, P011, P020, P101);
  gradN[5][1] = 0.;
  gradN[5][2] = gradControlCoefficientP110N200(P002, P110, P200, P011, P020, P101);
  gradN[5][3] = gradControlCoefficientP110N011(P002, P110, P200, P011, P020, P101);
  gradN[5][4] = gradControlCoefficientP110N101(P002, P110, P200, P011, P020, P101);
  gradN[5][5] = gradControlCoefficientP110N110(P002, P110, P200, P011, P020, P101);
}

// static void getDBezierdLagrange(double grad[6][6])
// {
//   // d_P200/dL_j : P200 = L200
//   grad[0][0] = 1.;
//   grad[0][1] = 0.;
//   grad[0][2] = 0.;
//   grad[0][3] = 0.;
//   grad[0][4] = 0.;
//   grad[0][5] = 0.;

//   // d_P020/dL_j : P020 = L020
//   grad[1][0] = 0.;
//   grad[1][1] = 1.;
//   grad[1][2] = 0.;
//   grad[1][3] = 0.;
//   grad[1][4] = 0.;
//   grad[1][5] = 0.;

//   // d_P020/dL_j : P002 = L002
//   grad[2][0] = 0.;
//   grad[2][1] = 0.;
//   grad[2][2] = 1.;
//   grad[2][3] = 0.;
//   grad[2][4] = 0.;
//   grad[2][5] = 0.;

//   // d_P110/dL_j : P110 = 1/2 * (4*L110 - L200 - L020)
//   grad[3][0] = -0.5;
//   grad[3][1] = -0.5;
//   grad[3][2] = 0.;
//   grad[3][3] = 2.;
//   grad[3][4] = 0.;
//   grad[3][5] = 0.;

//   // d_P011/dL_j : P011 = 1/2 * (4*L011 - L020 - L002)
//   grad[4][0] = 0.;
//   grad[4][1] = -0.5;
//   grad[4][2] = -0.5;
//   grad[4][3] = 0.;
//   grad[4][4] = 2.;
//   grad[4][5] = 0.;

//   // d_P101/dL_j : P101 = 1/2 * (4*L101 - L002 - L200)
//   grad[5][0] = -0.5;
//   grad[5][1] = 0.;
//   grad[5][2] = -0.5;
//   grad[5][3] = 0.;
//   grad[5][4] = 0.;
//   grad[5][5] = 2.;
// }

static bool checkValidity(const SPoint2 lagrangePoints[6], SPoint2 bezierPoints[6], double N[6])
{
  getBezierControlPoints(lagrangePoints, bezierPoints);
  getJacobianControlCoefficients(bezierPoints, N);
  for(int i = 0; i < 6; ++i) {
    if(N[i] <= 0.) return false;
  }
  return true;
}

static bool checkValidityOnModifiedElements(const std::vector<int> &whichElements)
{
  SPoint2 lagrangePoints[6], bezierPoints[6];
  double N[6];
  bool valid = true;
  for(size_t iElm = 0; iElm < whichElements.size(); ++iElm) {
    activeMesh->getCoord(activeConnectivity, whichElements[iElm], GEO_COORD);
    for(int i = 0; i < 6; ++i) {
      lagrangePoints[i] = SPoint2(GEO_COORD[i * 3 + 0], GEO_COORD[i * 3 + 1]);
    }
    valid = valid && checkValidity(lagrangePoints, bezierPoints, N);
  }
  return valid;
}

static double computeJacobian(const double uvw[3], const SPoint2 lagrangePoints[6],
                              SPoint2 bezierPoints[6], double N[6], double basis[6])
{
  getBezierControlPoints(lagrangePoints, bezierPoints);
  getJacobianControlCoefficients(bezierPoints, N);
  getP2BernsteinBasis(uvw, basis);
  double J = 0.;
  for(int i = 0; i < 6; ++i) {
    J += N[i] * basis[i];
  }
  return J;
}

static void computeGradientJacobian(const double uvw[3], const SPoint2 lagrangePoints[6],
                                    const int whichControlPoint, SPoint2 bezierPoints[6],
                                    SPoint2 gradN[6][6], double basis[6], double gradient[2])
{
  getBezierControlPoints(lagrangePoints, bezierPoints);
  getGradientJacobianControlCoefficients(bezierPoints, gradN);
  getP2BernsteinBasis(uvw, basis);
  SPoint2 gradJ(0., 0.);

  // The chain rule coefficient dLdP = 1 for the vertices control points
  // (they are identical) and 1/2 for edge control points
  // since they are given by permutations of
  // L110 = 1/4 * (P200 + 2*P110 + P020)
  // It's the inverse for dPdL.
  // double dBezierdLagrange = (whichControlPoint < 3) ? 1. : 2.; //  CHECK COEFFICIENT
  for(int i = 0; i < 6; ++i) {
    // gradJ += gradN[whichControlPoint][i] * basis[i] * dBezierdLagrange;
    gradJ += gradN[whichControlPoint][i] * basis[i];
  }
  gradient[0] = gradJ[0];
  gradient[1] = gradJ[1];
}

static void getGradientInterpolationErrorAtReferenceNode(
  const double xsi[2], const std::vector<double> &xPhysQuadratureNode,
  const std::vector<double> &xPhysControlPoint, const int whichControlPoint,
  const std::vector<double> &gradUAtQuadratureNode, const std::vector<double> &gradUAtControlPoint,
  double gradient[2])
{
  double phi[6] = {(1. - xsi[0] - xsi[1]) * (1. - 2. * xsi[0] - 2. * xsi[1]),
                   xsi[0] * (2. * xsi[0] - 1.),
                   xsi[1] * (2. * xsi[1] - 1.),
                   4. * xsi[0] * (1. - xsi[0] - xsi[1]),
                   4. * xsi[0] * xsi[1],
                   4. * xsi[1] * (1. - xsi[0] - xsi[1])};
  double phi0 = phi[whichControlPoint];
  gradient[0] = phi0 * (gradUAtQuadratureNode[0] - gradUAtControlPoint[0]);
  gradient[1] = phi0 * (gradUAtQuadratureNode[1] - gradUAtControlPoint[1]);
}

// void feNorm::computeInterpolationErrorGradientRochery(const int whichElements[2],
//                                                       const int whichControlPoint_localTag[2],
//                                                       double gradient[2])
// {
//   double uh, u, t = _solution->getCurrentTime();

//   double gradJ[2], gradE[2], uvw[3], basis[6], xsi[2], N[6];
//   std::vector<double> &xsi_quad = _spaces[0]->getRQuadraturePoints();
//   std::vector<double> &eta_quad = _spaces[0]->getSQuadraturePoints();
//   SPoint2 lagrangePoints[6], bezierPoints[6];
//   SPoint2 gradN[6][6];

//   gradient[0] = 0.;
//   gradient[1] = 0.;

//   std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
//   std::vector<double> dxds(3, 0.0);

//   // Only two elements are affected when moving an edge
//   for(int iElm = 0; iElm < 2; ++iElm) {
//     int elem = whichElements[iElm];
//     this->initializeLocalSolutionOnSpace(0, elem);

//     _spaces[0]->_mesh->getCoord(_cnc, elem, _geoCoord);

//     for(int i = 0; i < 6; ++i){
//       lagrangePoints[i] = SPoint2(_geoCoord[i * 3 + 0], _geoCoord[i * 3 + 1]);
//     }
//     std::vector<double> posControlPoint = {lagrangePoints[whichControlPoint_localTag[iElm]].x(),
//                                            lagrangePoints[whichControlPoint_localTag[iElm]].y(),
//                                            0.0};

//     for(int k = 0; k < _nQuad; ++k)
//     {
//       // Compute jacobian and its gradient
//       uvw[0] = 1. - xsi_quad[k] - eta_quad[k];
//       uvw[1] = xsi_quad[k];
//       uvw[2] = eta_quad[k];

//       double jac = computeJacobian(uvw, lagrangePoints, bezierPoints, N, basis);
//       computeGradientJacobian(uvw, lagrangePoints, whichControlPoint_localTag[iElm],
//       bezierPoints, gradN, basis, gradJ);

//       // Verification: compare with precomputed jacobian
//       // if(fabs(jac - _J[_nQuad * elem + k]) > 1e-13){
//       //   feInfo("Computed jacobians %+-1.6e and %+-1.6e", jac, _J[_nQuad * elem + k]);
//       //   exit(-1);
//       // }

//       // Verification: compare with derivatives of shape functions
//       double gradJ_Xj;
//       double gradJ_Yj;
//       for(int iii = 0; iii < 6; ++iii){
//         double dphi_jdxsi = _geoSpace->getdFunctiondrAtQuadNode(iii, k);
//         double dphi_jdeta = _geoSpace->getdFunctiondsAtQuadNode(iii, k);
//         _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(_geoCoord, k, dxdr);
//         _geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(_geoCoord, k, dxds);
//         double dxdxsi = dxdr[0];
//         double dydxsi = dxdr[1];
//         double dxdeta = dxds[0];
//         double dydeta = dxds[1];

//         gradJ_Xj = dydeta * dphi_jdxsi - dydxsi * dphi_jdeta;
//         gradJ_Yj = dxdxsi * dphi_jdeta - dxdeta * dphi_jdxsi;
//         // feInfo("Computed jacobian moving %d grad_x %+-1.6e and %+-1.6e", iii, gradJ[0],
//         gradJ_Xj);
//         // feInfo("Computed jacobian moving %d grad_y %+-1.6e and %+-1.6e", iii, gradJ[1],
//         gradJ_Yj);
//       }

//       xsi[0] = xsi_quad[k];
//       xsi[1] = eta_quad[k];
//       _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);
//       getGradientInterpolationErrorAtReferenceNode(xsi, _pos, posControlPoint,
//       whichControlPoint_localTag[iElm],
//         gradUAtQuadratureNode, gradUAtControlPoint, gradE);

//       uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
//       u = _scalarSolution->eval(t, _pos);
//       double error = u - uh;

//       // gradient[0] += _w[k] * (2. * error * gradE[0] * jac + error * gradJ[0]);
//       // gradient[1] += _w[k] * (2. * error * gradE[1] * jac + error * gradJ[1]);
//       gradient[0] += _w[k] * (2. * error * gradE[0] * jac + error * error * gradJ_Xj);
//       gradient[1] += _w[k] * (2. * error * gradE[1] * jac + error * error * gradJ_Yj);
//     }
//   }

//   // feInfo("Computed interpolation error gradient_x = %+-1.6e", gradient[0]);
//   // feInfo("Computed interpolation error gradient_y = %+-1.6e", gradient[1]);
// }

static std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
static std::vector<double> dxds(3, 0.0);
static std::vector<double> posControlPoint(3, 0.);
static std::vector<double> gradUAtQuadratureNode(2, 0.);
static std::vector<double> gradUAtControlPoint(2, 0.);
static std::vector<double> PLACEHOLDER(2, 0.);

void feNorm::computeInterpolationErrorGradient(const std::vector<int> &whichElements,
                                               const std::vector<int> &whichControlPoint_localTag,
                                               double gradient[2])
{
  double uh, u, uRec, t = _solution->getCurrentTime();

  double gradJ[2], gradESquared[2], gradE[2], xsi[2];
  std::vector<double> &xsi_quad = _spaces[0]->getRQuadraturePoints();
  std::vector<double> &eta_quad = _spaces[0]->getSQuadraturePoints();

  gradient[0] = 0.;
  gradient[1] = 0.;

  // Only two elements are affected when moving an edge
  for(int iElm = 0; iElm < 2; ++iElm) {
    int elem = whichElements[iElm];
    this->initializeLocalSolutionOnSpace(0, elem);

    _spaces[0]->_mesh->getCoord(_cnc, elem, _geoCoord);

    // Coordinates of the active Lagrange control point
    double X_j = _geoCoord[whichControlPoint_localTag[iElm] * 3 + 0];
    double Y_j = _geoCoord[whichControlPoint_localTag[iElm] * 3 + 1];
    posControlPoint[0] = X_j;
    posControlPoint[1] = Y_j;

// Evaluate gradient at control point
#if defined(GRAD_FROM_ESTIMATOR)
    gradUAtControlPoint[0] =
      _rec->evaluateRecovery(PPR::DERIVATIVE, 0, posControlPoint.data(), true);
    gradUAtControlPoint[1] =
      _rec->evaluateRecovery(PPR::DERIVATIVE, 1, posControlPoint.data(), true);
#else
    // Exact gradient
    activeExactSolutionGradient->eval(0., posControlPoint, gradUAtControlPoint);
#endif

    for(int k = 0; k < _nQuad; ++k) {
      // Compute jacobian and its gradient (jacobian was updated when moving the control point)
      double jac = _J[_nQuad * elem + k];

      double dphi_jdxsi = _geoSpace->getdFunctiondrAtQuadNode(whichControlPoint_localTag[iElm], k);
      double dphi_jdeta = _geoSpace->getdFunctiondsAtQuadNode(whichControlPoint_localTag[iElm], k);
      _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(_geoCoord, k, dxdr);
      _geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(_geoCoord, k, dxds);
      double dxdxsi = dxdr[0];
      double dydxsi = dxdr[1];
      double dxdeta = dxds[0];
      double dydeta = dxds[1];

      double gradJ_Xj = dydeta * dphi_jdxsi - dydxsi * dphi_jdeta;
      double gradJ_Yj = dxdxsi * dphi_jdeta - dxdeta * dphi_jdxsi;
      gradJ[0] = gradJ_Xj;
      gradJ[1] = gradJ_Yj;

      // Compute interpolation error and its gradient
      xsi[0] = xsi_quad[k];
      xsi[1] = eta_quad[k];
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);

      uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);

#if defined(ERROR_FROM_ESTIMATOR)
      uRec = _rec->evaluateRecovery(PPR::RECOVERY, 0, _pos.data(), USE_AVERAGED_EVALUATIONS);
      double error = uRec - uh;
#else
      u = _scalarSolution->eval(t, _pos);
      double error = u - uh;
#endif

// Evaluate the gradient at quad node
#if defined(GRAD_FROM_ESTIMATOR)
      gradUAtQuadratureNode[0] =
        _rec->evaluateRecovery(PPR::DERIVATIVE, 0, _pos.data(), USE_AVERAGED_EVALUATIONS);
      gradUAtQuadratureNode[1] =
        _rec->evaluateRecovery(PPR::DERIVATIVE, 1, _pos.data(), USE_AVERAGED_EVALUATIONS);
#else
      // Exact gradient
      activeExactSolutionGradient->eval(0., _pos, gradUAtQuadratureNode);
#endif

      getGradientInterpolationErrorAtReferenceNode(
        xsi, _pos, posControlPoint, whichControlPoint_localTag[iElm], gradUAtQuadratureNode,
        gradUAtControlPoint, gradE);

      gradESquared[0] = 2. * error * gradE[0];
      gradESquared[1] = 2. * error * gradE[1];

      gradient[0] += _w[k] * (gradESquared[0] * jac + error * error * gradJ[0]);
      gradient[1] += _w[k] * (gradESquared[1] * jac + error * error * gradJ[1]);
    }
  }
}

void computeInterpolationErrorGradient(const edgeAndVertexData &data, const double *modifiedCoord,
                                       double gradient[2])
{
  activeMesh->getVertexCoord(data.iVertexToMove_globalTag, INITIAL_COORD);

  // Get interpolated solution at initial vertex
  activeIntSpace->initializeAddressingVector(data.whichElements[0], ADR);
  int vertexDOF = ADR[data.iVertexToMove_localTag[0]];
  double uh = activeSolution->getSolAtDOF(vertexDOF);

  feStatus s = moveVertexAndRecomputeJacobians(data.whichElements, data.iVertexToMove_globalTag,
                                               modifiedCoord);

  // Compute solution at new location and assign it
  POS[0] = modifiedCoord[0];
  POS[1] = modifiedCoord[1];
  POS[2] = 0.;
#if defined(ERROR_FROM_ESTIMATOR)
  // Option 1
  double newUh =
    activeRecovery->evaluateRecovery(PPR::RECOVERY, 0, POS.data(), USE_AVERAGED_EVALUATIONS);
// Option 2
// activeRecovery->recomputeLeastSquareMatrix(data.iVertexToMove_globalTag);
// activeRecovery->recomputeRHSAndSolve(data.iVertexToMove_globalTag);
// activeRecovery->recomputeFirstDerivative(data.iVertexToMove_globalTag);
// double newUh = activeRecovery->evaluateRecoveryAtVertex(PPR::RECOVERY, 0,
// data.iVertexToMove_globalTag);
#else
  double newUh = activeExactSolution->eval(0., POS);
#endif
  activeSolution->setSolAtDOF(vertexDOF, newUh);

  // Compute gradient at modified vertex
  activeNorm->computeInterpolationErrorGradient(data.whichElements, data.iVertexToMove_localTag,
                                                gradient);

  double res = getNewError(data.nElements, data.whichElements);

  gradient[0] *= 1. / (2. * res);
  gradient[1] *= 1. / (2. * res);

  // Restore position and solution at initial vertex
  s = moveVertexAndRecomputeJacobians(data.whichElements, data.iVertexToMove_globalTag,
                                      INITIAL_COORD);
  activeSolution->setSolAtDOF(vertexDOF, uh);
}

// Callback to compute interpolation error during optimization
// The mesh is modified by moving vertex numbered "iVertexToMove"
// to the new position "modifiedCoord". The vertex position is restored
// before exiting.
// The mesh pointer in the norm has to match this mesh pointer.
double computeInterpolationErrorCallback(const edgeAndVertexData &data, const double *modifiedCoord,
                                         bool &valid)
{
  activeMesh->getVertexCoord(data.iVertexToMove_globalTag, INITIAL_COORD);

  // Get interpolated solution at initial vertex
  activeIntSpace->initializeAddressingVector(data.whichElements[0], ADR);
  int vertexDOF = ADR[data.iVertexToMove_localTag[0]];
  double uh = activeSolution->getSolAtDOF(vertexDOF);

  feStatus s = moveVertexAndRecomputeJacobians(data.whichElements, data.iVertexToMove_globalTag,
                                               modifiedCoord);

  valid = checkValidityOnModifiedElements(data.whichElements);
  if(!valid) return LARGE_VALUE_WHEN_FAILED;

  // Compute solution at new location and assign it
  POS[0] = modifiedCoord[0];
  POS[1] = modifiedCoord[1];
  POS[2] = 0.;
#if defined(ERROR_FROM_ESTIMATOR)
  double newUh =
    activeRecovery->evaluateRecovery(PPR::RECOVERY, 0, POS.data(), USE_AVERAGED_EVALUATIONS);
// activeRecovery->recomputeLeastSquareMatrix(data.iVertexToMove_globalTag);
// activeRecovery->recomputeRHSAndSolve(data.iVertexToMove_globalTag);
// activeRecovery->recomputeFirstDerivative(data.iVertexToMove_globalTag);
// double newUh = activeRecovery->evaluateRecoveryAtVertex(PPR::RECOVERY, 0,
// data.iVertexToMove_globalTag);
#else
  double newUh = activeExactSolution->eval(0., POS);
#endif
  activeSolution->setSolAtDOF(vertexDOF, newUh);

  double res = getNewError(data.nElements, data.whichElements);

  valid = checkValidityOnModifiedElements(data.whichElements);

  s = moveVertexAndRecomputeJacobians(data.whichElements, data.iVertexToMove_globalTag,
                                      INITIAL_COORD);

  // Restore solution at initial vertex
  activeSolution->setSolAtDOF(vertexDOF, uh);

  return res;
}

double computeInterpolationError(const edgeAndVertexData &data, double *modifiedCoord)
{
  bool valid = false;
  double err = computeInterpolationErrorCallback(data, modifiedCoord, valid);
  return valid ? err : LARGE_VALUE_WHEN_FAILED;
}

// Square [-1,1]
static double BBOX = 1.;
static double XMIN = -BBOX;
static double YMIN = -BBOX;
static double XMAX = BBOX;
static double YMAX = BBOX;

static int whichBoundary(PolyMesh::Vertex *v)
{
  double x = v->position[0];
  double y = v->position[1];
  double tol = 1e-4;
  bool isB = fabs(y - YMIN) < tol;
  bool isR = fabs(x - XMAX) < tol;
  bool isT = fabs(y - YMAX) < tol;
  bool isL = fabs(x - XMIN) < tol;
  bool isBBoxCorner = fabs(x - 1.1 * XMIN) < tol || fabs(x - 1.1 * XMAX) < tol ||
                      fabs(x - 1.1 * YMIN) < tol || fabs(x - 1.1 * YMAX) < tol;
  if(isB) return 0;
  if(isR) return 1;
  if(isT) return 2;
  if(isL) return 3;
  if(isBBoxCorner) return -2;
  return -1;
}

std::vector<double> MESH_INITIAL_COORD;
std::vector<double> MESH_NEW_COORD;
std::vector<double> MESH_INITIAL_SOL;
std::vector<double> MESH_NEW_SOL;
bool SHOULD_RESIZE = true;

// Same as callback above, but all the edges and vertices are modified
// The mesh pointer in the norm has to match this mesh pointer.
double computeInterpolationErrorCallback_EdgesAndVertices(
  const std::vector<double> &newPositions, const curvedAdaptData &adaptData,
  const std::vector<PolyMesh::Vertex *> &verticesToMove,
  const std::vector<PolyMesh::HalfEdge *> &edgesToMove)
{
  bool valid = true;

  if(SHOULD_RESIZE) {
    // Resize the big vector
    int totalsize = verticesToMove.size() + edgesToMove.size();
    MESH_INITIAL_COORD.resize(2 * totalsize); // 2 components
    MESH_NEW_COORD.resize(2 * totalsize); // 2 components
    MESH_INITIAL_SOL.resize(totalsize);
    MESH_NEW_SOL.resize(totalsize);
    SHOULD_RESIZE = false;
  }

  // Save current positions and solution of the feMesh (vertices and edges)
  int cnt = 0, cnt_u = 0;
  for(auto v : verticesToMove) {
    // Save position
    activeMesh->getVertexCoord(adaptData.v2globalTag.at(v), INITIAL_COORD);
    MESH_INITIAL_COORD[cnt] = INITIAL_COORD[0];
    cnt++;
    MESH_INITIAL_COORD[cnt] = INITIAL_COORD[1];
    cnt++;

    // Get interpolated solution at initial vertex
    activeIntSpace->initializeAddressingVector(adaptData.v2elements.at(v)[0], ADR);
    int vertexDOF = ADR[adaptData.v2localTags.at(v)[0]];
    MESH_INITIAL_SOL[cnt_u] = activeSolution->getSolAtDOF(vertexDOF);
    cnt_u++;
  }
  for(auto he : edgesToMove) {
    // Save position
    activeMesh->getVertexCoord(adaptData.he2midnodeGlobalTag.at(he), INITIAL_COORD);
    MESH_INITIAL_COORD[cnt] = INITIAL_COORD[0];
    cnt++;
    MESH_INITIAL_COORD[cnt] = INITIAL_COORD[1];
    cnt++;

    // Get interpolated solution at initial vertex
    activeIntSpace->initializeAddressingVector(adaptData.he2element.at(he), ADR);
    int vertexDOF = ADR[adaptData.he2midnodeLocalTag.at(he)];
    MESH_INITIAL_SOL[cnt_u] = activeSolution->getSolAtDOF(vertexDOF);
    cnt_u++;
  }

  // Move the DOFs at their new position
  feStatus s;
  cnt = 0;
  int cnt2 = 0;
  double modifiedCoord[2];
  std::vector<int> edgeElements(2);
  for(auto v : verticesToMove) {
    if(v->data == 1) {
      modifiedCoord[0] = newPositions[cnt];
      MESH_NEW_COORD[cnt2] = newPositions[cnt];
      cnt++;
      cnt2++;
      modifiedCoord[1] = newPositions[cnt];
      MESH_NEW_COORD[cnt2] = newPositions[cnt];
      cnt++;
      cnt2++;

    } else {
      int bdr = whichBoundary(v);
      if(bdr == 0 || bdr == 2) {
        modifiedCoord[0] = newPositions[cnt];
        MESH_NEW_COORD[cnt2] = newPositions[cnt];
        cnt++;
        cnt2++;
        modifiedCoord[1] = v->position[1];
        MESH_NEW_COORD[cnt2] = v->position[1];
        cnt2++;
      }
      if(bdr == 1 || bdr == 3) {
        modifiedCoord[0] = v->position[0];
        MESH_NEW_COORD[cnt2] = v->position[0];
        cnt2++;
        modifiedCoord[1] = newPositions[cnt];
        MESH_NEW_COORD[cnt2] = newPositions[cnt];
        cnt++;
        cnt2++;
      }
      if(bdr == -1) {
        feErrorMsg(FE_STATUS_ERROR, "Unexpected boundary %f - %f", v->position[0], v->position[1]);
      }
    }

    s = moveVertexAndRecomputeJacobians(adaptData.v2elements.at(v), adaptData.v2globalTag.at(v),
                                        modifiedCoord);
    if(s != FE_STATUS_OK) {
      feInfo("Could not move vertex during computeInterpolationError");
      exit(-1);
    }
    if(s != FE_STATUS_OK) {
      return std::numeric_limits<double>::infinity();
    }
    valid &= checkValidityOnModifiedElements(adaptData.v2elements.at(v));
  }
  for(auto he : edgesToMove) {
    modifiedCoord[0] = newPositions[cnt];
    MESH_NEW_COORD[cnt2] = newPositions[cnt];
    cnt++;
    cnt2++;
    modifiedCoord[1] = newPositions[cnt];
    MESH_NEW_COORD[cnt2] = newPositions[cnt];
    cnt++;
    cnt2++;
    edgeElements[0] = adaptData.he2element.at(he);
    edgeElements[1] = adaptData.he2element.at(he->opposite);
    s = moveVertexAndRecomputeJacobians(edgeElements, adaptData.he2midnodeGlobalTag.at(he),
                                        modifiedCoord);
    if(s != FE_STATUS_OK) {
      feInfo("Could not move edge during computeInterpolationError");
      exit(-1);
    }
    if(s != FE_STATUS_OK) {
      return std::numeric_limits<double>::infinity();
    }
    valid &= checkValidityOnModifiedElements(edgeElements);
  }

  // for(int i = 0; i < MESH_NEW_COORD.size(); ++i) {
  //   if(fabs(MESH_NEW_COORD[i] - MESH_INITIAL_COORD[i]) > 1e-10) {
  //     feInfo("Not identical : %+-1.12e vs %+-1.12e", MESH_NEW_COORD[i], MESH_INITIAL_COORD[i]);
  //     exit(-1);
  //   }
  //   else { feInfo("Identical"); };
  // }

  double res;
  if(valid) {
    // Compute solution at new location and assign it
    cnt = 0, cnt_u = 0;
    POS[2] = 0.;
    for(auto v : verticesToMove) {
      if(v->data == 1) {
        POS[0] = newPositions[cnt];
        cnt++;
        POS[1] = newPositions[cnt];
        cnt++;
      } else {
        int bdr = whichBoundary(v);
        if(bdr == 0 || bdr == 2) {
          POS[0] = newPositions[cnt];
          cnt++;
          POS[1] = v->position[1];
        }
        if(bdr == 1 || bdr == 3) {
          POS[0] = v->position[0];
          POS[1] = newPositions[cnt];
          cnt++;
        }
        if(bdr == -1) {
          feErrorMsg(FE_STATUS_ERROR, "Unexpected boundary %f - %f", v->position[0],
                     v->position[1]);
        }
      }

      double newUh = activeExactSolution->eval(0., POS);
      MESH_NEW_SOL[cnt_u] = newUh;
      cnt_u++;
      activeIntSpace->initializeAddressingVector(adaptData.v2elements.at(v)[0], ADR);
      int vertexDOF = ADR[adaptData.v2localTags.at(v)[0]];
      activeSolution->setSolAtDOF(vertexDOF, newUh);
    }
    for(auto he : edgesToMove) {
      POS[0] = newPositions[cnt];
      cnt++;
      POS[1] = newPositions[cnt];
      cnt++;
      double newUh = activeExactSolution->eval(0., POS);
      MESH_NEW_SOL[cnt_u] = newUh;
      cnt_u++;
      activeIntSpace->initializeAddressingVector(adaptData.he2element.at(he), ADR);
      int vertexDOF = ADR[adaptData.he2midnodeLocalTag.at(he)];
      activeSolution->setSolAtDOF(vertexDOF, newUh);
    }

    // double res = getNewError(data.nElements, data.whichElements);
    res = activeNorm->compute();
  } else {
    res = std::numeric_limits<double>::infinity();
  }

  // Move the DOFs back to initial position
  cnt = 0, cnt_u = 0;
  for(auto v : verticesToMove) {
    INITIAL_COORD[0] = MESH_INITIAL_COORD[cnt];
    cnt++;
    INITIAL_COORD[1] = MESH_INITIAL_COORD[cnt];
    cnt++;
    s = moveVertexAndRecomputeJacobians(adaptData.v2elements.at(v), adaptData.v2globalTag.at(v),
                                        INITIAL_COORD);
    if(s != FE_STATUS_OK) {
      feInfo("Could not move vertex BACK during computeInterpolationError");
      exit(-1);
    }
    valid &= checkValidityOnModifiedElements(edgeElements);
    // if(!valid) return std::numeric_limits<double>::infinity();

    // Restore solution
    activeIntSpace->initializeAddressingVector(adaptData.v2elements.at(v)[0], ADR);
    int vertexDOF = ADR[adaptData.v2localTags.at(v)[0]];
    activeSolution->setSolAtDOF(vertexDOF, MESH_INITIAL_SOL[cnt_u]);
    cnt_u++;
  }
  for(auto he : edgesToMove) {
    INITIAL_COORD[0] = MESH_INITIAL_COORD[cnt];
    cnt++;
    INITIAL_COORD[1] = MESH_INITIAL_COORD[cnt];
    cnt++;
    edgeElements[0] = adaptData.he2element.at(he);
    edgeElements[1] = adaptData.he2element.at(he->opposite);
    s = moveVertexAndRecomputeJacobians(edgeElements, adaptData.he2midnodeGlobalTag.at(he),
                                        INITIAL_COORD);
    if(s != FE_STATUS_OK) {
      feInfo("Could not move edge BACK during computeInterpolationError");
      exit(-1);
    }
    valid &= checkValidityOnModifiedElements(edgeElements);
    // if(!valid) return std::numeric_limits<double>::infinity();

    // Restore solution
    activeIntSpace->initializeAddressingVector(adaptData.he2element.at(he), ADR);
    int vertexDOF = ADR[adaptData.he2midnodeLocalTag.at(he)];
    activeSolution->setSolAtDOF(vertexDOF, MESH_INITIAL_SOL[cnt_u]);
    cnt_u++;
  }

  if(!valid) {
    feInfo("DOFs relocation did not cause a negative jacobian, but there is at least one invalid "
           "element. Should not happen?");
    exit(-1);
  }

  feInfo("ERROR IS %+-1.8e", res);
  return res;
}

// Return the tag of the midnode whose coordinates match to the prescribed tolerance
// Used in gmsh to identify a HalfEdge to a vertex tag in feNG.
void getMidnodeTags(const SPoint2 edge[2], const double tol, int &elementTag, int &localTag,
                    int &globalTag)
{
  SPoint2 midnode = (edge[0] + edge[1]) * 0.5;
  std::vector<Vertex> &meshVertices = activeMesh->getVertices();
  bool wasFound = false;
  for(auto v : meshVertices) {
    if(fabs(midnode[0] - v.x()) < tol && fabs(midnode[1] - v.y()) < tol) {
      globalTag = activeMesh->_verticesMap[v.getTag()];
      wasFound = true;
      break;
    }
  }
  if(!wasFound) {
    feErrorMsg(FE_STATUS_ERROR, "Midnode (%f, %f) was not found in the mesh!", midnode[0],
               midnode[1]);
    exit(-1);
  }

  // Get the element on which the edge lies and the local tag of the midnode
  feCncGeo *cnc = activeMesh->getCncGeoByName("Domaine");
  std::vector<double> geoCoord(cnc->getNumVerticesPerElem() * 3, 0.);
  for(int iElm = 0; iElm < cnc->getNumElements(); ++iElm) {
    activeMesh->getCoord(cnc, iElm, geoCoord);
    for(int iEdge = 0; iEdge < 3; ++iEdge) {
      double x0 = geoCoord[3 * iEdge + 0];
      double y0 = geoCoord[3 * iEdge + 1];
      double x1 = geoCoord[3 * ((iEdge + 1) % 3) + 0];
      double y1 = geoCoord[3 * ((iEdge + 1) % 3) + 1];
      if(fabs(edge[0].x() - x0) < tol && fabs(edge[0].y() - y0) < tol &&
         fabs(edge[1].x() - x1) < tol && fabs(edge[1].y() - y1) < tol) {
        elementTag = iElm;
        localTag = iEdge + 3;
        return;
      }
    }
  }
  feErrorMsg(FE_STATUS_ERROR, "Edge (%f, %f) - (%f, %f) was not found in the mesh!", edge[0].x(),
             edge[0].y(), edge[1].x(), edge[1].y());
  exit(-1);
}

void getPolyMeshVertexTags(const SPoint2 &p, const double tol, std::vector<int> &elementTags,
                           std::vector<int> &localTags, int &globalTag)
{
  // Check for vertex global tag
  std::vector<Vertex> &meshVertices = activeMesh->getVertices();
  bool wasFound = false;
  for(auto v : meshVertices) {
    if(fabs(p[0] - v.x()) < tol && fabs(p[1] - v.y()) < tol) {
      globalTag = activeMesh->_verticesMap[v.getTag()];
      wasFound = true;
      break;
    }
  }
  if(!wasFound) {
    feErrorMsg(FE_STATUS_ERROR, "PolyMesh::Vertex (%f, %f) was not found in the mesh!", p[0], p[1]);
    exit(-1);
  }

  // Check for elements containing p
  feCncGeo *cnc = activeMesh->getCncGeoByName("Domaine");
  std::vector<double> geoCoord(cnc->getNumVerticesPerElem() * 3, 0.);
  for(int iElm = 0; iElm < cnc->getNumElements(); ++iElm) {
    activeMesh->getCoord(cnc, iElm, geoCoord);
    for(int iVertex = 0; iVertex < 3; ++iVertex) {
      double x = geoCoord[3 * iVertex + 0];
      double y = geoCoord[3 * iVertex + 1];
      if(fabs(p[0] - x) < tol && fabs(p[1] - y) < tol) {
        elementTags.push_back(iElm);
        localTags.push_back(iVertex);
      }
    }
  }
}
#endif