#include "feMesh.h"
#include "feNorm.h"
#include "feAdaptMesh.h"
#include "feMetricTools.h"

#if defined(HAVE_GMSH)
#include "gmsh.h"
#endif

extern int FE_VERBOSE;

int gmshWasInitialized = false;

feRecovery *activeRecovery;
feSpace *activeIntSpace;
feNumber *activeNumbering;
feSolution *activeSolution;
feFunction *activeExactSolution;
feVectorFunction *activeExactSolutionGradient; 

feMesh *activeMesh;
feCncGeo *activeConnectivity;
feNorm *activeNorm;

static double INITIAL_COORD[3];
static double LARGE_VALUE_WHEN_FAILED = 1e10;
static std::vector<double> ELEM_ERROR_SQUARED;
static double INITIAL_ERROR_SQUARED;
static std::vector<double> GEO_COORD(18, 0.);

feStatus moveVertexAndRecomputeJacobians(const std::vector<int> &whichElements, const int iVertexToMove_globalTag, const double *modifiedCoord)
{
  activeMesh->setVertexCoord(iVertexToMove_globalTag, modifiedCoord);
  bool failed = false;
  for(size_t i = 0; i < whichElements.size(); ++i) {
    feStatus s = activeConnectivity->recomputeElementJacobian(whichElements[i]);
    failed = (s != FE_STATUS_OK);
  }
  if(failed)
    return FE_STATUS_FAILED;
  return FE_STATUS_OK;
}

void computeInterpolationErrorOnEachElement()
{
  // Store the initial interpolation error on each element
  ELEM_ERROR_SQUARED.resize(activeConnectivity->getNumElements());
  INITIAL_ERROR_SQUARED = 0.;
  for(size_t i = 0; i < activeConnectivity->getNumElements(); ++i) {
    ELEM_ERROR_SQUARED[i] = activeNorm->computeSquaredErrorOnElement(i);
    INITIAL_ERROR_SQUARED += ELEM_ERROR_SQUARED[i];
  }
}

void setNewError(const int nElements, const std::vector<int> &whichElements)
{
  for(size_t i = 0; i < nElements; ++i) {
    INITIAL_ERROR_SQUARED -= ELEM_ERROR_SQUARED[whichElements[i]];
  }
  for(size_t i = 0; i < nElements; ++i) {
    ELEM_ERROR_SQUARED[whichElements[i]] = activeNorm->computeSquaredErrorOnElement(whichElements[i]);
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
    res += activeNorm->computeSquaredErrorOnElement(whichElements[i]);
  }
  return sqrt(res);
}

static std::vector<feInt> ADR(6);
static std::vector<double> POS(3);

void applyCurvatureToFeMesh(const edgeAndVertexData &data, const double *modifiedCoord)
{
  // std::string posfile = "DomaineBefore.pos";
  // static_cast<feMesh2DP1 *>(activeMesh)->drawConnectivityToPOSfile("Domaine", posfile);

  // Get interpolated solution at initial vertex
  activeIntSpace->initializeAddressingVector(data.whichElements[0], ADR);
  int vertexDOF = ADR[data.iVertexToMove_localTag[0]];
  double uh = activeSolution->getSolAtDOF(vertexDOF);

  feStatus s = moveVertexAndRecomputeJacobians(data.whichElements, data.iVertexToMove_globalTag, modifiedCoord);

  // Compute solution at new location and assign it
  POS[0] = modifiedCoord[0];
  POS[1] = modifiedCoord[1];
  POS[2] = 0.;
  double newUh = activeExactSolution->eval(0., POS);
  activeSolution->setSolAtDOF(vertexDOF, newUh);

  // feStatus s = moveVertexAndRecomputeJacobians(data.whichElements, data.iVertexToMove_globalTag, modifiedCoord);
  // activeSolution->initializeEssentialBC(activeMesh);

  // Update error vector and total error (substract old error and add new error)
  setNewError(data.nElements, data.whichElements);
  // for(size_t i = 0; i < data.nElements; ++i) {
  //   INITIAL_ERROR_SQUARED -= ELEM_ERROR_SQUARED[data.whichElements[i]];
  //   ELEM_ERROR_SQUARED[data.whichElements[i]] = activeNorm->computeLpNormOnElement(2, true, data.whichElements[i]);
  //   INITIAL_ERROR_SQUARED += ELEM_ERROR_SQUARED[data.whichElements[i]];
  // }

  // feInfo("NEW ERROR IS %+-1.8e for modifiedCoord %f - %f", sqrt(INITIAL_ERROR_SQUARED), modifiedCoord[0], modifiedCoord[1]);
  // INITIAL_ERROR_SQUARED = INITIAL_ERROR_SQUARED - ELEM_ERROR_SQUARED[data.whichElements[0]] - ELEM_ERROR_SQUARED[data.whichElements[1]];
  // ELEM_ERROR_SQUARED[data.whichElements[0]] = activeNorm->computeLpNormOnElement(2, true, data.whichElements[0]);
  // ELEM_ERROR_SQUARED[data.whichElements[1]] = activeNorm->computeLpNormOnElement(2, true, data.whichElements[1]);
  // INITIAL_ERROR_SQUARED = INITIAL_ERROR_SQUARED + ELEM_ERROR_SQUARED[data.whichElements[0]] + ELEM_ERROR_SQUARED[data.whichElements[1]];

  // posfile = "DomaineAfter.pos";
  // static_cast<feMesh2DP1 *>(activeMesh)->drawConnectivityToPOSfile("Domaine", posfile);
  if(s != FE_STATUS_OK){
    feInfo("Modified vertex caused at least one negative jacobian :/");
    exit(-1);
  }
}

void feMesh2DP1::drawConnectivityToPOSfile(const std::string &cncName, const std::string &fileName)
{
  feCncGeo *cnc = this->getCncGeoByName(cncName);
  if(cnc){
    FILE* myfile = fopen(fileName.data(), "w");
    fprintf(myfile, "View\"%s\"{\n", cncName.data());
    std::vector<double> geoCoord(cnc->getNumVerticesPerElem() * 3, 0.);
    for(int iElm = 0; iElm < cnc->getNumElements(); ++iElm)
    {
      this->getCoord(cnc, iElm, geoCoord);
      cnc->writeElementToPOS(myfile, geoCoord, 1.);
    }
    fprintf(myfile, "};\n"); fclose(myfile);
  }
}

static void getP2BernsteinBasis(const double uvw[3], double basis[6])
{
  double u = uvw[0];
  double v = uvw[1];
  double w = uvw[2];
  basis[0] = u*u;
  basis[1] = 2.*u*v;
  basis[2] = v*v;
  basis[3] = 2.*v*w;
  basis[4] = w*w;
  basis[5] = 2.*w*u;
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

static double controlCoefficientN200(const SPoint2 &P200, const SPoint2 &P011,
                                     const SPoint2 &P020, const SPoint2 &P101,
                                     const SPoint2 &P002, const SPoint2 &P110)
{
  return 4. * computeDeterminant(P200 - P110, P200 - P101);
}

static double controlCoefficientN110(const SPoint2 &P200, const SPoint2 &P011,
                                     const SPoint2 &P020, const SPoint2 &P101,
                                     const SPoint2 &P002, const SPoint2 &P110)
{
  return 2. * computeDeterminant(P200 - P101, P020 - P011) + 2. * computeDeterminant(P110 - P011, P110 - P101);
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
  // N[3] = 2. * computeDeterminant(P200 - P101, P020 - P011) + 2. * computeDeterminant(P110 - P011, P110 - P101); // N110
  // N[4] = 2. * computeDeterminant(P020 - P110, P002 - P101) + 2. * computeDeterminant(P011 - P101, P011 - P110); // N011
  // N[5] = 2. * computeDeterminant(P002 - P011, P200 - P110) + 2. * computeDeterminant(P101 - P110, P101 - P011); // N101

  // Or using permutations:
  N[0] = controlCoefficientN200(P200, P011, P020, P101, P002, P110);
  N[1] = controlCoefficientN200(P020, P101, P002, P110, P200, P011);
  N[2] = controlCoefficientN200(P002, P110, P200, P011, P020, P101);
  N[3] = controlCoefficientN110(P200, P011, P020, P101, P002, P110);
  N[4] = controlCoefficientN110(P020, P101, P002, P110, P200, P011);
  N[5] = controlCoefficientN110(P002, P110, P200, P011, P020, P101);
}

inline SPoint2 orthogonal(const SPoint2 &P)
{
  return SPoint2(P[1], -P[0]);
}

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

static bool checkValidity(const SPoint2 lagrangePoints[6],
                          SPoint2 bezierPoints[6],
                          double N[6])
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

static double computeJacobian(const double uvw[3],
                              const SPoint2 lagrangePoints[6],
                              SPoint2 bezierPoints[6],
                              double N[6],
                              double basis[6])
{
  getBezierControlPoints(lagrangePoints, bezierPoints);
  getJacobianControlCoefficients(bezierPoints, N);
  getP2BernsteinBasis(uvw, basis);
  double J = 0.;
  for(int i = 0; i < 6; ++i){
    J += N[i] * basis[i];
  }
  return J;
}

static void computeGradientJacobian(const double uvw[3],
                                    const SPoint2 lagrangePoints[6],
                                    const int whichControlPoint,
                                    SPoint2 bezierPoints[6],
                                    SPoint2 gradN[6][6],
                                    double basis[6],
                                    double gradient[2])
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
  double dBezierdLagrange = (whichControlPoint < 3) ? 1. : 2.; //  CHECK COEFFICIENT
  for(int i = 0; i < 6; ++i){
    // gradJ += gradN[whichControlPoint][i] * basis[i] * dBezierdLagrange;
    gradJ += gradN[whichControlPoint][i] * basis[i];
  }
  gradient[0] = gradJ[0];
  gradient[1] = gradJ[1];
}

static std::vector<double> gradUAtQuadratureNode(2, 0.);
static std::vector<double> gradUAtControlPoint(2, 0.);

static void getGradientInterpolationErrorAtReferenceNode(const double xsi[2],
                                                         const std::vector<double> &xPhysQuadratureNode,
                                                         const std::vector<double> &xPhysControlPoint,
                                                         const int whichControlPoint,
                                                         double gradient[2])
{
  double phi[6] = {(1. - xsi[0] - xsi[1]) * (1. - 2. * xsi[0] - 2. * xsi[1]),
                   xsi[0] * (2. * xsi[0] - 1.),
                   xsi[1] * (2. * xsi[1] - 1.),
                   4. * xsi[0] * (1. - xsi[0] - xsi[1]),
                   4. * xsi[0] * xsi[1],
                   4. * xsi[1] * (1. - xsi[0] - xsi[1])};
  double phi0 = phi[whichControlPoint];
  activeExactSolutionGradient->eval(0., xPhysQuadratureNode, gradUAtQuadratureNode);
  activeExactSolutionGradient->eval(0., xPhysControlPoint, gradUAtControlPoint);

  gradient[0] = phi0 * (gradUAtQuadratureNode[0] - gradUAtControlPoint[0]);
  gradient[1] = phi0 * (gradUAtQuadratureNode[1] - gradUAtControlPoint[1]);
}

void feNorm::computeInterpolationErrorGradientRochery(const int whichElements[2],
                                                      const int whichControlPoint_localTag[2],
                                                      double gradient[2])
{
  double uh, u, t = _solution->getCurrentTime();

  double gradJ[2], gradE[2], uvw[3], basis[6], xsi[2], N[6];
  std::vector<double> &xsi_quad = _spaces[0]->getRQuadraturePoints();
  std::vector<double> &eta_quad = _spaces[0]->getSQuadraturePoints();
  SPoint2 lagrangePoints[6], bezierPoints[6];
  SPoint2 gradN[6][6];

  gradient[0] = 0.;
  gradient[1] = 0.;

  std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
  std::vector<double> dxds(3, 0.0);

  // Only two elements are affected when moving an edge
  for(int iElm = 0; iElm < 2; ++iElm) {
    int elem = whichElements[iElm];
    this->initializeLocalSolutionOnSpace(0, elem);

    _spaces[0]->_mesh->getCoord(_cnc, elem, _geoCoord);

    for(int i = 0; i < 6; ++i){
      lagrangePoints[i] = SPoint2(_geoCoord[i * 3 + 0], _geoCoord[i * 3 + 1]);
    }
    std::vector<double> posControlPoint = {lagrangePoints[whichControlPoint_localTag[iElm]].x(),
                                           lagrangePoints[whichControlPoint_localTag[iElm]].y(),
                                           0.0};

    for(int k = 0; k < _nQuad; ++k)
    {
      // Compute jacobian and its gradient
      uvw[0] = 1. - xsi_quad[k] - eta_quad[k];
      uvw[1] = xsi_quad[k];
      uvw[2] = eta_quad[k];

      double jac = computeJacobian(uvw, lagrangePoints, bezierPoints, N, basis);
      computeGradientJacobian(uvw, lagrangePoints, whichControlPoint_localTag[iElm], bezierPoints, gradN, basis, gradJ);

      // Verification: compare with precomputed jacobian
      // if(fabs(jac - _J[_nQuad * elem + k]) > 1e-13){
      //   feInfo("Computed jacobians %+-1.6e and %+-1.6e", jac, _J[_nQuad * elem + k]);
      //   exit(-1);
      // }

      // Verification: compare with derivatives of shape functions
      double gradJ_Xj;
      double gradJ_Yj;
      for(int iii = 0; iii < 6; ++iii){
        double dphi_jdxsi = _geoSpace->getdFunctiondrAtQuadNode(iii, k);
        double dphi_jdeta = _geoSpace->getdFunctiondsAtQuadNode(iii, k);
        _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(_geoCoord, k, dxdr);
        _geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(_geoCoord, k, dxds);
        double dxdxsi = dxdr[0];
        double dydxsi = dxdr[1];
        double dxdeta = dxds[0];
        double dydeta = dxds[1];

        gradJ_Xj = dydeta * dphi_jdxsi - dydxsi * dphi_jdeta;
        gradJ_Yj = dxdxsi * dphi_jdeta - dxdeta * dphi_jdxsi;
        // feInfo("Computed jacobian moving %d grad_x %+-1.6e and %+-1.6e", iii, gradJ[0], gradJ_Xj);
        // feInfo("Computed jacobian moving %d grad_y %+-1.6e and %+-1.6e", iii, gradJ[1], gradJ_Yj);
      }

      xsi[0] = xsi_quad[k];
      xsi[1] = eta_quad[k];
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);
      getGradientInterpolationErrorAtReferenceNode(xsi, _pos, posControlPoint, whichControlPoint_localTag[iElm], gradE);
      
      uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
      u = _scalarSolution->eval(t, _pos);
      double error = u - uh;

      // gradient[0] += _w[k] * (2. * error * gradE[0] * jac + error * gradJ[0]);
      // gradient[1] += _w[k] * (2. * error * gradE[1] * jac + error * gradJ[1]);
      gradient[0] += _w[k] * (2. * error * gradE[0] * jac + error * error * gradJ_Xj);
      gradient[1] += _w[k] * (2. * error * gradE[1] * jac + error * error * gradJ_Yj);
    }
  }

  // feInfo("Computed interpolation error gradient_x = %+-1.6e", gradient[0]);
  // feInfo("Computed interpolation error gradient_y = %+-1.6e", gradient[1]);
}

static std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
static std::vector<double> dxds(3, 0.0);
static std::vector<double> posControlPoint(3, 0.);

void feNorm::computeInterpolationErrorGradient(const std::vector<int> &whichElements,
                                               const std::vector<int> &whichControlPoint_localTag,
                                               double gradient[2])
{
  double uh, u, t = _solution->getCurrentTime();

  double gradJ[2], gradESquared[2], gradE[2], uvw[3], basis[6], xsi[2];
  std::vector<double> &xsi_quad = _spaces[0]->getRQuadraturePoints();
  std::vector<double> &eta_quad = _spaces[0]->getSQuadraturePoints();

  gradient[0] = 0.;
  gradient[1] = 0.;

  // Only two elements are affected when moving an edge
  for(int iElm = 0; iElm < 2; ++iElm)
  {
    int elem = whichElements[iElm];
    this->initializeLocalSolutionOnSpace(0, elem);

    _spaces[0]->_mesh->getCoord(_cnc, elem, _geoCoord);

    // Coordinates of the active Lagrange control point
    double X_j = _geoCoord[whichControlPoint_localTag[iElm] * 3 + 0];
    double Y_j = _geoCoord[whichControlPoint_localTag[iElm] * 3 + 1];
    posControlPoint[0] = X_j;
    posControlPoint[1] = Y_j;

    for(int k = 0; k < _nQuad; ++k)
    {
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
      u = _scalarSolution->eval(t, _pos);
      double error = u - uh;

      getGradientInterpolationErrorAtReferenceNode(xsi, _pos, posControlPoint, whichControlPoint_localTag[iElm], gradE);
      gradESquared[0] = 2. * error * gradE[0];
      gradESquared[1] = 2. * error * gradE[1];

      gradient[0] += _w[k] * (gradESquared[0] * jac + error * error * gradJ[0]);
      gradient[1] += _w[k] * (gradESquared[1] * jac + error * error * gradJ[1]);
    }
  }
}

void computeInterpolationErrorGradient(const edgeAndVertexData &data,
                                       const double *modifiedCoord,
                                       double gradient[2])
{
  activeMesh->getVertexCoord(data.iVertexToMove_globalTag, INITIAL_COORD);

  // Get interpolated solution at initial vertex
  activeIntSpace->initializeAddressingVector(data.whichElements[0], ADR);
  int vertexDOF = ADR[data.iVertexToMove_localTag[0]];
  double uh = activeSolution->getSolAtDOF(vertexDOF);

  feStatus s = moveVertexAndRecomputeJacobians(data.whichElements, data.iVertexToMove_globalTag, modifiedCoord);

  // Compute solution at new location and assign it
  POS[0] = modifiedCoord[0];
  POS[1] = modifiedCoord[1];
  POS[2] = 0.;
  double newUh = activeExactSolution->eval(0., POS);
  activeSolution->setSolAtDOF(vertexDOF, newUh);

  /////////////////////////////////////////////////////////
  // Compute gradient at modified vertex
  activeNorm->computeInterpolationErrorGradient(data.whichElements, data.iVertexToMove_localTag, gradient);

  // Recompute error only on the N modified elements
  // FIXME
  // double res = INITIAL_ERROR_SQUARED;
  // for(size_t i = 0; i < data.nElements; ++i) {
  //   res -= ELEM_ERROR_SQUARED[data.whichElements[i]];
  // // }
  // // for(size_t i = 0; i < data.nElements; ++i) {
  //   res += activeNorm->computeLpNormOnElement(2, true, data.whichElements[i]);
  // }
  // res = sqrt(res);
  double res = getNewError(data.nElements, data.whichElements);

  // double e0_squared = activeNorm->computeLpNormOnElement(2, true, data.whichElements[0]);
  // double e1_squared = activeNorm->computeLpNormOnElement(2, true, data.whichElements[1]);
  // double res2 = INITIAL_ERROR_SQUARED - ELEM_ERROR_SQUARED[data.whichElements[0]] - ELEM_ERROR_SQUARED[data.whichElements[1]] + e0_squared + e1_squared;
  // res2 = sqrt(res2);
  // if(fabs(res - res2) > 1e-14) {
  //   feInfo("res = %+-1.8e - res2 = %+-1.8e - diff = %+-1.8e", res, res2, fabs(res-res2));
  //   exit(-1);
  // }

  gradient[0] *= 1. / (2.*res);
  gradient[1] *= 1. / (2.*res);
  // feInfo("Computed gradient %+-1.16e - %+-1.16e", gradient[0], gradient[1]);
  /////////////////////////////////////////////////////////

  // Restore position and solution at initial vertex
  s = moveVertexAndRecomputeJacobians(data.whichElements, data.iVertexToMove_globalTag, INITIAL_COORD);
  activeSolution->setSolAtDOF(vertexDOF, uh);
}

// Callback to compute interpolation error during optimization
// The mesh is modified by moving vertex numbered "iVertexToMove"
// to the new position "modifiedCoord". The vertex position is restored
// before exiting.
// The mesh pointer in the norm has to match this mesh pointer.
double computeInterpolationErrorCallback(const edgeAndVertexData &data,
                                         const double *modifiedCoord,
                                         bool &valid)
{
  activeMesh->getVertexCoord(data.iVertexToMove_globalTag, INITIAL_COORD);

  // Get interpolated solution at initial vertex
  activeIntSpace->initializeAddressingVector(data.whichElements[0], ADR);
  int vertexDOF = ADR[data.iVertexToMove_localTag[0]];
  double uh = activeSolution->getSolAtDOF(vertexDOF);

  feStatus s = moveVertexAndRecomputeJacobians(data.whichElements, data.iVertexToMove_globalTag, modifiedCoord);

  // Compute solution at new location and assign it
  POS[0] = modifiedCoord[0];
  POS[1] = modifiedCoord[1];
  POS[2] = 0.;
  double newUh = activeExactSolution->eval(0., POS);
  activeSolution->setSolAtDOF(vertexDOF, newUh);

  // double res = pow(activeNorm->compute(), 2);
  // double res = activeNorm->compute();

  // Recompute error only on the N modified elements
  // double res = INITIAL_ERROR_SQUARED;
  // for(size_t i = 0; i < data.nElements; ++i) {
  //   res -= ELEM_ERROR_SQUARED[data.whichElements[i]];
  //   res += activeNorm->computeLpNormOnElement(2, true, data.whichElements[i]);
  // }
  // res = sqrt(res);
  double res = getNewError(data.nElements, data.whichElements);

  // // // Compute error only on the two modified elements
  // double e0_squared = activeNorm->computeLpNormOnElement(2, true, data.whichElements[0]);
  // double e1_squared = activeNorm->computeLpNormOnElement(2, true, data.whichElements[1]);
  // double res2 = INITIAL_ERROR_SQUARED - ELEM_ERROR_SQUARED[data.whichElements[0]] - ELEM_ERROR_SQUARED[data.whichElements[1]] + e0_squared + e1_squared;
  // res2 = sqrt(res2);

  // if(fabs(res - res2) > 1e-14) {
  //   feInfo("In Error callback : res = %+-1.8e - res2 = %+-1.8e - diff = %+-1.8e", res, res2, fabs(res-res2));
  //   exit(-1);
  // }

  valid = checkValidityOnModifiedElements(data.whichElements);

  s = moveVertexAndRecomputeJacobians(data.whichElements, data.iVertexToMove_globalTag, INITIAL_COORD);

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

// Return the tag of the midnode whose coordinates match to the prescribed tolerance
// Used in gmsh to identify a HalfEdge to a vertex tag in feNG.
void getMidnodeTags(const SPoint2 edge[2], const double tol, int &elementTag, int &localTag, int &globalTag)
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
  if(!wasFound){
    feErrorMsg(FE_STATUS_ERROR, "Midnode (%f, %f) was not found in the mesh!",
      midnode[0], midnode[1]);
    exit(-1);
  }

  // Get the element on which the edge lies and the local tag of the midnode
  feCncGeo *cnc = activeMesh->getCncGeoByName("Domaine");
  std::vector<double> geoCoord(cnc->getNumVerticesPerElem() * 3, 0.);
  for(int iElm = 0; iElm < cnc->getNumElements(); ++iElm){
    activeMesh->getCoord(cnc, iElm, geoCoord);
    for(int iEdge = 0; iEdge < 3; ++iEdge){
      double x0 = geoCoord[3 * iEdge + 0];
      double y0 = geoCoord[3 * iEdge + 1];
      double x1 = geoCoord[3 * ((iEdge + 1) % 3) + 0];
      double y1 = geoCoord[3 * ((iEdge + 1) % 3) + 1];
      if  (fabs(edge[0].x() - x0) < tol && fabs(edge[0].y() - y0) < tol
        && fabs(edge[1].x() - x1) < tol && fabs(edge[1].y() - y1) < tol)
      {
        elementTag = iElm;
        localTag = iEdge + 3;
        return;
      }
    }
  }
  feErrorMsg(FE_STATUS_ERROR, "Edge (%f, %f) - (%f, %f) was not found in the mesh!",
    edge[0].x(), edge[0].y(), edge[1].x(), edge[1].y());
  exit(-1);
}

void getPolyMeshVertexTags(const SPoint2 &p, const double tol, 
  std::vector<int> &elementTags, std::vector<int> &localTags, int &globalTag)
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
  if(!wasFound){
    feErrorMsg(FE_STATUS_ERROR, "PolyMesh::Vertex (%f, %f) was not found in the mesh!",
      p[0], p[1]);
    exit(-1);
  }

  // Check for elements containing p
  feCncGeo *cnc = activeMesh->getCncGeoByName("Domaine");
  std::vector<double> geoCoord(cnc->getNumVerticesPerElem() * 3, 0.);
  for(int iElm = 0; iElm < cnc->getNumElements(); ++iElm){
    activeMesh->getCoord(cnc, iElm, geoCoord);
    for(int iVertex = 0; iVertex < 3; ++iVertex){
      double x = geoCoord[3 * iVertex + 0];
      double y = geoCoord[3 * iVertex + 1];
      if(fabs(p[0] - x) < tol && fabs(p[1] - y) < tol) {
        elementTags.push_back(iElm);
        localTags.push_back(iVertex);
      }
    }
  }
}

inline feStatus checkMMGcall(std::string &command)
{
  // From mmg/src/mmg2d/libmmgtypes.h
  int MMG_SUCCESS = 0;

  auto status = system(command.data());
  if(status < 0) {
    return feErrorMsg(FE_STATUS_ERROR, "System call to MMG2D failed with system status %d", status);
  } else {
    if(WIFEXITED(status)) {
      if(WEXITSTATUS(status) == MMG_SUCCESS)
        // Process ran normally and MMG exited normally
        return FE_STATUS_OK;
      else
        // Process ran normally, but MMG failed to adapt for some reason
        return feErrorMsg(FE_STATUS_ERROR, "Adaptation with MMG2D failed with exit code %d",
                          WIFEXITED(status));
    } else
      return feErrorMsg(FE_STATUS_ERROR, "System call to MMG2D failed");
  }
}

feStatus feMesh2DP1::adapt(feNewRecovery *recoveredField, feMetricOptions &options,
                           const std::vector<feSpace *> &spaces,
                           const std::vector<feSpace *> &essentialSpaces,
                           feSpace *spaceForAdaptation,
                           feSolution *discreteSolution,
                           feFunction *exactSolution,
                           feVectorFunction *exactGradient,
                           bool curve, bool isBackmeshP2, feRecovery *oldRecovery)
{

#if defined(HAVE_GMSH)

  activeMesh = this;

  // Step 1: Set the gmsh background mesh on which the metric tensors are computed
  if(!gmshWasInitialized) {

    // Get max threads *before* initializing gmsh
    #if defined(HAVE_OMP)
    int maxNumThreads = omp_get_max_threads();
    #endif

    gmsh::initialize();

    #if defined(HAVE_OMP)
    gmsh::option::setNumber("General.NumThreads", maxNumThreads);
    #endif

    gmshWasInitialized = true;
  }

  if(FE_VERBOSE == VERBOSE_NONE) gmsh::option::setNumber("General.Verbosity", 2);

  gmsh::open(options.backgroundMeshfile);

  int elementOrder = getGeometricInterpolantDegree(recoveredField->_cnc->getInterpolant());
  if(!curve && elementOrder > 1)
  {
    // Set the mesh to P1.
    // This is because MMG only takes a P1 mesh for aniso adaptation.
    // The map from gmsh P1 nodeTags to the feMesh P2 vertices is 
    // built when computing the metrics in feMetric::createVertex2NodeMap
    // If we curve, the background mesh should be P2, so dont change order.
    gmsh::model::mesh::setOrder(1);
  }

  gmsh::model::getCurrent(options.modelForMetric);
  options.isGmshModelReady = true;

  feMetric metricField(recoveredField, options);

  // Step 2: Create aniso mesh
  // Save directly to .msh to preserve Physical Entities
  std::string cmd1 = "mmg2d " + options.mmgInputMeshfile + " -hgrad -1 -o " +
  options.mmgOutputMeshfile;
  if(FE_VERBOSE == VERBOSE_NONE) {
    // Write MMG console outputs to logMMG.txt
    cmd1 += " > logMMG.txt";
  }

  // Assign old feRecovery structure (for curved adaptation, temporary)
  metricField.setRecovery(oldRecovery);

  // The back feMesh can be P2, but we need to write a P1 view in the gmsh model
  // to give mmg for aniso adaptation. So what decides on the number of vertices
  // when metrics are computed is whether we curve (and dont have to give a mesh to mmg)
  // and not whether isP2backmesh is true
  metricField._nVerticesPerElmOnBackmesh = curve ? 6 : 3;

  if(curve && !isBackmeshP2) {
    return feErrorMsg(FE_STATUS_ERROR, "Backmesh must be P2 when curving the mesh.");
  }

  // Compute the metric tensor field and write metrics to options.mmgInputMeshfile
  tic();
  feStatus s = metricField.computeMetrics();
  options.userValue = metricField._options.userValue;
  feInfoCond(FE_VERBOSE > 0, "\t\tComputed metric tensors in %f s", toc());
  if(s != FE_STATUS_OK) {
    gmsh::finalize();
    return s;
  }

  if(!curve){
    // Create aniso mesh
    s = checkMMGcall(cmd1);
    if(s != FE_STATUS_OK) {
      gmsh::finalize();
      return s;
    }

    // Open adapted mesh (DEPRECATED: compute the next metric field on this mesh)
    gmsh::open(options.mmgOutputMeshfile);
    gmsh::model::getCurrent(options.modelForMetric);
    metricField.setGmshMetricModel(options.modelForMetric);

    // Step 3: Check physical entities
    gmsh::vectorpair physicalGroups;
    gmsh::model::getPhysicalGroups(physicalGroups);

    std::map<std::pair<int, int>, std::vector<int> > entitiesForPhysical;

    // Get the entities for existing physical groups
    for(auto pair : physicalGroups) {
      std::vector<int> entities;
      gmsh::model::getEntitiesForPhysicalGroup(pair.first, pair.second, entities);
      entitiesForPhysical[pair] = entities;
    }

    // Remove all the physical groups
    gmsh::model::removePhysicalGroups();

    // Re-add stored physical groups
    for(auto pair : _physicalEntitiesDescription) {
      int dim = pair.first.first;
      int tag = pair.first.second;
      std::string name = pair.second;

      bool OK = false;
      for(auto p : physicalGroups) {
        if(p.first == dim && p.second == tag) {
          gmsh::model::addPhysicalGroup(dim, entitiesForPhysical[p], tag);
          gmsh::model::setPhysicalName(dim, tag, name);
          OK = true;
        }
      }
      if(!OK) {
        return feErrorMsg(FE_STATUS_ERROR,
                          "Physical Entity \"%s\" with (dim,tag) = (%d,%d)"
                          " could not be reassigned after mesh adaptation :/",
                          name.data(), dim, tag);
      }
    }

    // gmsh::model::mesh::reverse();

    if(elementOrder > 1)
    {
      // Reset the mesh to Pn
      gmsh::model::mesh::setOrder(elementOrder);
    }

    // Write adapted anisotropic straight mesh with correct Physical Entities
    gmsh::option::setNumber("Mesh.MshFileVersion", 4.1);
    gmsh::write("aniso" + options.adaptedMeshName);
  }

  // Step 4: Interpolate solution on new mesh (TODO)
  // feMesh2DP1 newMesh("afterPhysical.msh");
  // feMetaNumber numbering(&newMesh, spaces, essentialSpaces);

  // Curve after a few aniso adaptations
  if(curve) {
    // gmsh::open(options.mmgInputMeshfile);
    // gmsh::model::getCurrent(options.modelForMetric);

    // gmsh::model::setCurrent(options.modelForMetric);
    // gmsh::model::setCurrent(options.mmgOutputMeshfile);

    // Tools to compute interpolation error while curving
    feNorm *norm;
    createNorm(norm, L2_ERROR, {spaceForAdaptation}, discreteSolution, exactSolution, exactGradient);
    activeNorm = norm;
    activeConnectivity = activeMesh->getCncGeoByName("Domaine");
    activeIntSpace = spaceForAdaptation;
    activeSolution = discreteSolution;
    activeExactSolution = exactSolution;
    activeExactSolutionGradient = exactGradient;

    computeInterpolationErrorOnEachElement();

    if(elementOrder > 1)
    {
      // GFace2PolyMesh only takes P1 meshes
      gmsh::model::mesh::setOrder(1);
    }

    int faceTag;
    gmsh::vectorpair dimTags;
    gmsh::model::getEntities(dimTags, 2);
    if(dimTags.size() > 1) {
      return feErrorMsg(FE_STATUS_ERROR, "Gmsh model has more than one surface");
    } else {
      faceTag = dimTags[0].second;
    }

    std::vector<int> viewTags;
    gmsh::view::getTags(viewTags);
    feInfo("There are %d views in the gmsh model : ", viewTags.size());
    for(auto val : viewTags) {
      feInfo("View %d with index %d", val, gmsh::view::getIndex(val));
    }

    metricField.setMetricViewTag(viewTags[0]);

    std::vector<double> pts;
    int onlyGenerateVertices = false;

    // Check there is an "inside" callback
    if(options.insideCallback == nullptr){
      return feErrorMsg(FE_STATUS_ERROR,
        "Provide an \"inside\" callback to peform curved mesh adaptation");
    }

    tic();
    computePointsUsingScaledCrossFieldPlanarP2(
      options.modelForMetric.c_str(),
      options.modelForMetric.c_str(),
      metricField.getMetricViewTag(),
      faceTag,
      pts,
      computeInterpolationError,
      computeInterpolationErrorGradient,
      applyCurvatureToFeMesh,
      getMidnodeTags,
      getPolyMeshVertexTags,
      options.insideCallback,
      nullptr,
      onlyGenerateVertices, 
      evaluateFieldFromRecoveryCallback, 
      (void *) oldRecovery,
      interpolateMetricP1WithDerivativesWrapper,
      interpolateMetricP1Wrapper,
      interpolateMetricAndDerivativeOnP2EdgeWrapper,
      interpolateMetricP1Wrapper1D,
      interpolateMetricAndDerivativeOnP2EdgeWrapper1D,
      (void *) &metricField);
    feInfoCond(FE_VERBOSE > 0, "\t\tGenerated curved mesh in %f s", toc());

    gmsh::model::setCurrent(options.modelForMetric.c_str());
    gmsh::option::setNumber("Mesh.MshFileVersion", 4.1);
    gmsh::write("curvedMesh.msh");
  }

#endif
  return FE_STATUS_OK;
}

/* Creates a curved mesh based on
    - the geometry stored in the active Gmsh model with name metricOptions.gmshModel
    - the metric field stored as a view in this model
    - the "inside" callback metricOptions.insideCallback, returning true if a point is inside the geometry
*/
// void createCurvedMesh(feFunction *solExact, feMetaNumber *metaNumber, feSolution *sol,
//                       feSpace *intSpace, feRecovery *recovery, feMetric *metric,
//                       feMetricOptions &metricOptions, int onlyGenerateVertices, int nLoopsAnisoMesh,
//                       bool curve)
// {
// #if defined(HAVE_GMSH)
//   // if(metricOptions.isGmshModelReady){
//   std::vector<double> pts;
//   int faceTag = 0;
//   activeRecovery = recovery;
//   activeIntSpace = intSpace;
//   activeNumbering = metaNumber->getNumbering(intSpace->getFieldID());
//   activeSolution = sol;
//   exactSolution = solExact;

//   gmsh::model::add(metricOptions.modelForMesh);
//   gmsh::model::setCurrent(metricOptions.modelForMesh);

//   // Aniso mesh with MMG (used to get the boundary vertices only)
//   // std::string cmd = "mmg2d " + metricOptions.mmgInputMeshfile + " -hgrad 3 -o " +
//   // metricOptions.mmgOutputMeshfile; std::string cmd1 = "mmg2d " +
//   // metricOptions.mmgInputMeshfile + " -hgrad 10 -o tmp.mesh";
//   std::string cmd1 = "mmg2d " + metricOptions.mmgInputMeshfile + " -hgrad -1 -o tmp.mesh";
//   system(cmd1.c_str());
//   std::string cmd2 = "gmsh tmp.mesh -o " + metricOptions.mmgOutputMeshfile + " -0";
//   system(cmd2.c_str());
//   std::string cmd3 = "gmsh " + metricOptions.mmgOutputMeshfile + " &";
//   // system(cmd3.c_str());

//   // Loop a few times
//   for(int i = 0; i < nLoopsAnisoMesh; ++i) {
//     gmsh::open(metricOptions.mmgOutputMeshfile);
//     gmsh::model::getCurrent(metricOptions.modelForMetric);
//     metric->setGmshMetricModel(metricOptions.modelForMetric);
//     // feInfo("Showing gmsh models at iter %d", i);
//     // gmsh::fltk::run();
//     metric->computeMetrics();
//     system(cmd1.c_str());
//     system(cmd2.c_str());
//     // system(cmd3.c_str());
//   };

//   // system(cmd3.c_str());

//   gmsh::clear();
//   gmsh::open(metricOptions.mmgOutputMeshfile);
//   // gmsh::open("thegmshModel.msh");
//   gmsh::model::getCurrent(metricOptions.modelForMesh);

//   gmsh::write("beforeVersion.msh");
//   gmsh::option::setNumber("Mesh.MshFileVersion", 4.1);
//   gmsh::write("afterVersion.msh");

//   // Determine if the mesh is reversed (MMG seems to reverse the mesh sometimes, unless I'm doing
//   // something wrong) Get quadrature rule and interpolation functions on the adapted mesh
//   int triP1 = gmsh::model::mesh::getElementType("Triangle", 1);
//   std::vector<double> localCoord;
//   std::vector<double> weights;
//   gmsh::model::mesh::getIntegrationPoints(triP1, "Gauss4", localCoord, weights);
//   // Get the jacobians
//   std::vector<double> jac, det, points;
//   gmsh::model::mesh::getJacobians(triP1, localCoord, jac, det, points);

//   // for(int i = 0; i < 20; ++i)
//   //   feInfo("det1 = %f", det[i]);

//   // Check the first determinant
//   if(det[0] < 0) {
//     feInfo("Adapted aniso mesh is reversed : reversing the numbering");
//     gmsh::model::mesh::reverse();
//   }

//   gmsh::model::mesh::getJacobians(triP1, localCoord, jac, det, points);
//   // for(int i = 0; i < 20; ++i)
//   //   feInfo("det2 = %f", det[i]);

//   // Get geometric entities and assign physical entities (should be improved)
//   gmsh::vectorpair dimTags;
//   gmsh::model::getEntities(dimTags);
//   std::vector<int> entities1D, entities2D;
//   for(auto p : dimTags) {
//     std::cout << p.first << " - " << p.second << std::endl;
//     if(p.first == 1) entities1D.push_back(p.second);
//     if(p.first == 2) entities2D.push_back(p.second);
//   }

//   gmsh::write("beforePhysical.msh");
//   // Add physical entities
//   gmsh::model::addPhysicalGroup(1, entities1D, 1);
//   gmsh::model::setPhysicalName(1, 1, "Bord");
//   gmsh::model::addPhysicalGroup(2, entities2D, 2);
//   gmsh::model::setPhysicalName(2, 2, "Domaine");
//   gmsh::write("afterPhysical.msh");

//   // For some reason (?) the physical tags added above can be negative, although
//   // the mesh is numbered in counterclockwise orientation.
//   // If it's the case, reverse the mesh and keep the negative physical tags.
//   std::vector<int> physicalTags;
//   bool atLeastOnePositive = false, atLeastOneNegative = false;
//   for(auto e : entities1D) {
//     gmsh::model::getPhysicalGroupsForEntity(1, e, physicalTags);
//     if(physicalTags[0] >= 0) {
//       atLeastOnePositive = true;
//     } else {
//       atLeastOneNegative = true;
//     }
//   }
//   for(auto e : entities2D) {
//     gmsh::model::getPhysicalGroupsForEntity(2, e, physicalTags);
//     if(physicalTags[0] >= 0) {
//       atLeastOnePositive = true;
//     } else {
//       atLeastOneNegative = true;
//     }
//   }
//   if(atLeastOnePositive && atLeastOneNegative) {
//     // If this happens then I'm very confused
//     feWarning(
//       "Some physical tags added to the MMG mesh are positive, whereas some are negative...");
//     exit(-1);
//   }

//   if(atLeastOneNegative) {
//     // All are negative : reverse the mesh
//     feInfo("Physical tags added to the MMG mesh are negative but numbering is positive :"
//            " reversing the numbering to match the signs of the physical tags");
//     gmsh::model::mesh::reverse();
//     exit(-1);
//   }
//   gmsh::write("afterPhysicalCheck.msh");

//   gmsh::write(metricOptions.adaptedMeshName);

//   // Curve after a few aniso adaptations
//   if(curve) {
//     gmsh::open(metricOptions.mmgInputMeshfile);
//     gmsh::model::getCurrent(metricOptions.modelForMetric);

//     gmsh::model::setCurrent(metricOptions.modelForMetric);
//     gmsh::model::setCurrent(metricOptions.modelForMesh);

//     gmsh::model::getEntities(dimTags, 2);
//     if(dimTags.size() > 1) {
//       feWarning("Gmsh model has more than one surface");
//     } else {
//       faceTag = dimTags[0].second;
//     }

//     std::vector<int> viewTags;
//     gmsh::view::getTags(viewTags);
//     feInfo("There are %d views in the gmsh model : ", viewTags.size());
//     for(auto val : viewTags) {
//       feInfo("View %d with index %d", val, gmsh::view::getIndex(val));
//     }

//     metric->setMetricViewTag(viewTags[0]);

//     // computePointsUsingScaledCrossFieldPlanarP2(
//     //   metricOptions.modelForMetric.c_str(), metricOptions.modelForMesh.c_str(),
//     //   metric->getMetricViewTag(), faceTag, pts, errorSquaredCallback, metricOptions.insideCallback,
//     //   nullptr, onlyGenerateVertices, evaluateFieldFromRecoveryCallback, (void *)recovery,
//     // interpolateMetricP1WithDerivativesWrapper, interpolateMetricP1Wrapper,
//     // interpolateMetricAndDerivativeOnP2EdgeWrapper, interpolateMetricP1Wrapper1D,
//     // interpolateMetricAndDerivativeOnP2EdgeWrapper1D,
//     // (void *)metric);

//     gmsh::write(metricOptions.adaptedMeshName);
//   }

// #else
//   printf("In feAdaptMesh : Error - Gmsh is required to generate curved meshes.\n");
// #endif
// }