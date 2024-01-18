#include "feMetric.h"
#include "../contrib/Eigen/QR"

// Wrapper to give gmsh
void interpolateMetricP1Callback(void *metricPtr, const double *x, Eigen::Matrix2d &M,
                           Eigen::Matrix2d &dMdx, Eigen::Matrix2d &dMdy)
{
  static_cast<feMetric*>(metricPtr)->interpolateMetricP1(x, M, dMdx, dMdy);
}

// Interpolate metric on P1 background mesh
void feMetric::interpolateMetricP1(const double *x, Eigen::Matrix2d &M,
                                   Eigen::Matrix2d &dMdx, Eigen::Matrix2d &dMdy)
{
  feInfo("Please use P2 interpolation!");
  exit(-1);
  // Locate point in the feMesh
  int elm;
  double UVW[3];
  bool isFound = static_cast<feMesh2DP1 *>(_recoveredFields[0]->_mesh)->locateVertex(x, elm, UVW);
  if(!isFound) {
    feWarning("In interpolateMetricP2 : Point (%f, %f) was not found in the mesh.", x[0], x[1]);
    M = Eigen::Matrix2d::Identity();
    dMdx = Eigen::Matrix2d::Zero();
    dMdy = Eigen::Matrix2d::Zero();
  } else {

    double xsi[2] = {UVW[0], UVW[1]};

    // Get the vertex tags from the feMesh
    int gmshNodeTags[3];
    for(int i = 0; i < 3; ++i) {
      int vertexTag = _recoveredFields[0]->_cnc->getVertexConnectivity(elm, i);
      gmshNodeTags[i] = _sequentialTag2nodeTag[vertexTag];
    }

    MetricTensor res;
    this->logEuclidianP1Interpolation(xsi, _logMetricTensorAtNodetags.at(gmshNodeTags[0]),
                                           _logMetricTensorAtNodetags.at(gmshNodeTags[1]),
                                           _logMetricTensorAtNodetags.at(gmshNodeTags[2]), res);
    M(0,0) = res(0,0);
    M(0,1) = res(0,1);
    M(1,0) = res(1,0);
    M(1,1) = res(1,1);

    gradLogEuclidianP1Interpolation(xsi, elm, _logMetricTensorAtNodetags.at(gmshNodeTags[0]),
                                              _logMetricTensorAtNodetags.at(gmshNodeTags[1]),
                                              _logMetricTensorAtNodetags.at(gmshNodeTags[2]), dMdx, dMdy);
  }
}

// Gradient of the interpolated log-metric on P1 background mesh
// with respect to front mesh coordinates
void feMetric::gradLogEuclidianP1Interpolation(const double xsi[2],
                                               const int element,
                                               const MetricTensor &logM0,
                                               const MetricTensor &logM1,
                                               const MetricTensor &logM2,
                                               Eigen::Matrix2d &dMdx,
                                               Eigen::Matrix2d &dMdy)
{
  // Step 1: compute inverse of grad_(xsiHat) of phi_pHat
  // Inverse of jacobian of P1 transformation in the background mesh
  Triangle *t = _recoveredFields[0]->_mesh->_elements[element];
  double x0 = t->getVertex(0)->x();
  double y0 = t->getVertex(0)->y();
  double x1 = t->getVertex(1)->x();
  double y1 = t->getVertex(1)->y();
  double x2 = t->getVertex(2)->x();
  double y2 = t->getVertex(2)->y();
  Eigen::Matrix2d Jac;
  Jac(0,0) = x1 - x0;
  Jac(0,1) = x2 - x0;
  Jac(1,0) = y1 - y0;
  Jac(1,1) = y2 - y0;

  // Step 2: compute gradient of BGM shape functions with respect to 
  // the coordinates in the front mesh. This is dphiHat/dxsiHat * dxsiHat/dx,
  // where hat are the background mesh quantities.

  // Hardcoded gradient of the 3 P1 shape functions with respect to xj, yj in front mesh
  double phi[3] = {1 - xsi[0] - xsi[1], xsi[0], xsi[1]};
  double dphidr[3] = {-1., 1., 0.};
  double dphids[3] = {-1., 0., 1.};
  Eigen::Matrix2d invJ = Jac.inverse();
  Eigen::RowVector2d dN0dx, dN1dx, dN2dx;
  dN0dx << dphidr[0], dphids[0]; dN0dx = dN0dx * invJ;
  dN1dx << dphidr[1], dphids[1]; dN1dx = dN1dx * invJ;
  dN2dx << dphidr[2], dphids[2]; dN2dx = dN2dx * invJ;

  MetricTensor L    = logM0 * phi[0] + logM1 * phi[1] + logM2 * phi[2];
  MetricTensor dLdx = logM0 * dN0dx(0) + logM1 * dN1dx(0) + logM2 * dN2dx(0);
  MetricTensor dLdy = logM0 * dN0dx(1) + logM1 * dN1dx(1) + logM2 * dN2dx(1);

  // Step 3: compute gradient of eigenvalues and eigenvectors of L　from Aparicio-Estrems et al, High order
  // metric interpolation for curved r-adaptation by distortion minimization)
  Eigen::Matrix2d Leigen, dLdxeigen, dLdyeigen;
  Leigen(0,0) = L(0,0);
  Leigen(0,1) = L(0,1);
  Leigen(1,0) = L(1,0);
  Leigen(1,1) = L(1,1);

  dLdxeigen(0,0) = dLdx(0,0);
  dLdxeigen(0,1) = dLdx(0,1);
  dLdxeigen(1,0) = dLdx(1,0);
  dLdxeigen(1,1) = dLdx(1,1);

  dLdyeigen(0,0) = dLdy(0,0);
  dLdyeigen(0,1) = dLdy(0,1);
  dLdyeigen(1,0) = dLdy(1,0);
  dLdyeigen(1,1) = dLdy(1,1);
  Eigen::EigenSolver<Eigen::Matrix2d> es;
  Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType eigenvalues;
  Eigen::EigenSolver<Eigen::Matrix2d>::EigenvectorsType eigenvectors;

  es.compute(Leigen, true);
  eigenvalues = es.eigenvalues();
  Eigen::Vector2d lambda = eigenvalues.real();

  eigenvectors = es.eigenvectors();
  Eigen::Matrix2d U = eigenvectors.real();
  // std::cout << "eigenvectors = " << eigenvectors <<std::endl;
  // std::cout << "U = " << U <<std::endl;
  Eigen::Vector2d u1 = U.col(0);
  Eigen::Vector2d u2 = U.col(1);

  // std::cout << "lambda = " << lambda <<std::endl;
  // std::cout << "u1 = " << u1 <<std::endl;
  // std::cout << "u2 = " << u2 <<std::endl;

  // Gradient des valeurs propres d_j lambda_l
  double dl1dx = u1.transpose() * dLdxeigen * u1;
  double dl1dy = u1.transpose() * dLdyeigen * u1;

  double dl2dx = u2.transpose() * dLdxeigen * u2;
  double dl2dy = u2.transpose() * dLdyeigen * u2;

  // std::cout << "dl1dx = " << dl1dx <<std::endl;
  // std::cout << "dl2dx = " << dl2dx <<std::endl;
  // std::cout << "dl1dy = " << dl1dy <<std::endl;
  // std::cout << "dl2dy = " << dl2dy <<std::endl;

  // Gradient des vecteurs propres d_j u_l
  Eigen::Matrix2d L1 = Leigen - lambda(0) * Eigen::Matrix2d::Identity();
  Eigen::Matrix2d dL1dx = dLdxeigen - dl1dx * Eigen::Matrix2d::Identity();
  Eigen::Matrix2d dL1dy = dLdyeigen - dl1dy * Eigen::Matrix2d::Identity();

  Eigen::Matrix2d L2 = Leigen - lambda(1) * Eigen::Matrix2d::Identity();
  Eigen::Matrix2d dL2dx = dLdxeigen - dl2dx * Eigen::Matrix2d::Identity();
  Eigen::Matrix2d dL2dy = dLdyeigen - dl2dy * Eigen::Matrix2d::Identity();
  
  Eigen::Vector2d du1dx = L1.completeOrthogonalDecomposition().solve(-dL1dx * u1);
  Eigen::Vector2d du1dy = L1.completeOrthogonalDecomposition().solve(-dL1dy * u1);

  // Correction de https://math.stackexchange.com/questions/2689374/derivative-of-eigenvectors-of-a-symmetric-matrix-valued-function
  Eigen::Vector2d w = du1dx;
  Eigen::Vector2d v = u1;
  du1dx = w - (v.transpose() * w) * v;

  w = du1dy;
  du1dy = w - (v.transpose() * w) * v;

  Eigen::Vector2d du2dx = L2.completeOrthogonalDecomposition().solve(-dL2dx * u2);
  Eigen::Vector2d du2dy = L2.completeOrthogonalDecomposition().solve(-dL2dy * u2);

  w = du2dx;
  v = u2;
  du2dx = w - (v.transpose() * w) * v;

  w = du2dy;
  du2dy = w - (v.transpose() * w) * v;

  // std::cout << "du1dx = " << du1dx <<std::endl;
  // std::cout << "du2dx = " << du2dx <<std::endl;
  // std::cout << "du1dy = " << du1dy <<std::endl;
  // std::cout << "du2dy = " << du2dy <<std::endl;

  Eigen::Matrix2d dUdx = Eigen::Matrix2d::Zero(), dUdy = Eigen::Matrix2d::Zero();
  dUdx.col(0) = du1dx;
  dUdx.col(1) = du2dx;
  dUdy.col(0) = du1dy;
  dUdy.col(1) = du2dy;

  Eigen::Matrix2d expD;
  expD << exp(lambda(0)), 0., 0., exp(lambda(1));
  Eigen::Matrix2d dDdx;
  dDdx << dl1dx, 0., 0., dl2dx;
  Eigen::Matrix2d dDdy;
  dDdy << dl1dy, 0., 0., dl2dy;

  // std::cout << expD <<std::endl;
  // std::cout << dDdx <<std::endl;

  Eigen::Matrix2d dexpDdx = expD * dDdx;
  Eigen::Matrix2d dexpDdy = expD * dDdy;

  dMdx = dUdx * expD * U.transpose() + U * dexpDdx * U.transpose() + U * expD * dUdx.transpose();
  dMdy = dUdy * expD * U.transpose() + U * dexpDdy * U.transpose() + U * expD * dUdy.transpose();
}

// Wrapper to give gmsh without the metric derivatives
void interpolateMetricP2CallbackWithoutDerivatives(void *metricPtr, const double *x, Eigen::Matrix2d &M)
{
  static_cast<feMetric*>(metricPtr)->interpolateMetricP2(x, M);
}

// Wrapper to give gmsh with the metric derivatives
void interpolateMetricP2CallbackWithDerivatives(void *metricPtr, const double *x, Eigen::Matrix2d &M,
                           Eigen::Matrix2d &dMdx, Eigen::Matrix2d &dMdy)
{
  static_cast<feMetric*>(metricPtr)->interpolateMetricP2(x, M, dMdx, dMdy);
}

thread_local double UVW[3];
thread_local int GMSHNODETAGS_P2[6];

//
// Main metric interpolation function
// Interpolate metric on P2 background mesh without computing metric derivatives
//
void feMetric::interpolateMetricP2(const double *x, Eigen::Matrix2d &M)
{
  // Locate point in the feMesh
  int elm;
  bool isFound = static_cast<feMesh2DP1 *>(_recoveredFields[0]->_mesh)->locateVertex(x, elm, UVW);
  if(!isFound) {
    // Point was not found in the mesh
  } else {

    double xsi[2] = {UVW[0], UVW[1]};

    // Get the vertex tags from the feMesh
    for(int i = 0; i < 6; ++i) {
      int vertexTag = _recoveredFields[0]->_cnc->getVertexConnectivity(elm, i);
      // GMSHNODETAGS_P2[i] = _sequentialTag2nodeTag[vertexTag];
      GMSHNODETAGS_P2[i] = _sequentialTag2nodeTagVec[vertexTag];
    }

    this->logEuclidianP2Interpolation(xsi, _logMetricTensorAtNodetags_eigen[GMSHNODETAGS_P2[0]-1],
                                           _logMetricTensorAtNodetags_eigen[GMSHNODETAGS_P2[1]-1],
                                           _logMetricTensorAtNodetags_eigen[GMSHNODETAGS_P2[2]-1],
                                           _logMetricTensorAtNodetags_eigen[GMSHNODETAGS_P2[3]-1],
                                           _logMetricTensorAtNodetags_eigen[GMSHNODETAGS_P2[4]-1],
                                           _logMetricTensorAtNodetags_eigen[GMSHNODETAGS_P2[5]-1], M);
  }
}

thread_local Eigen::Matrix2d L, DLDX, DLDY;
thread_local double L1, L2, DL1DX, DL1DY, DL2DX, DL2DY;
thread_local Eigen::Vector2d U1, U2, DU1DX, DU1DY, DU2DX, DU2DY;

//
// Main metric interpolation function
// Interpolate metric on P2 background mesh and compute metric derivatives
//
void feMetric::interpolateMetricP2(const double *x, Eigen::Matrix2d &M,
                                   Eigen::Matrix2d &dMdx, Eigen::Matrix2d &dMdy)
{
  // Locate point in the feMesh
  int elm;
  bool isFound = static_cast<feMesh2DP1 *>(_recoveredFields[0]->_mesh)->locateVertex(x, elm, UVW);
  if(!isFound) {
    // Point was not found in the mesh
  } else {

    double xsi[2] = {UVW[0], UVW[1]};

    // Get the vertex tags from the feMesh
    for(int i = 0; i < 6; ++i) {
      int vertexTag = _recoveredFields[0]->_cnc->getVertexConnectivity(elm, i);
      // GMSHNODETAGS_P2[i] = _sequentialTag2nodeTag[vertexTag];
      GMSHNODETAGS_P2[i] = _sequentialTag2nodeTagVec[vertexTag];
    }

    this->logEuclidianP2Interpolation(xsi, _logMetricTensorAtNodetags_eigen[GMSHNODETAGS_P2[0]-1],
                                           _logMetricTensorAtNodetags_eigen[GMSHNODETAGS_P2[1]-1],
                                           _logMetricTensorAtNodetags_eigen[GMSHNODETAGS_P2[2]-1],
                                           _logMetricTensorAtNodetags_eigen[GMSHNODETAGS_P2[3]-1],
                                           _logMetricTensorAtNodetags_eigen[GMSHNODETAGS_P2[4]-1],
                                           _logMetricTensorAtNodetags_eigen[GMSHNODETAGS_P2[5]-1], M);

    // Eigen::Matrix2d L, dLdx, dLdy;
    // double l1, l2;
    // double dl1dx;
    // double dl1dy;
    // double dl2dx;
    // double dl2dy;
    // Eigen::Vector2d u1, u2;
    // Eigen::Vector2d du1dx;
    // Eigen::Vector2d du1dy;
    // Eigen::Vector2d du2dx;
    // Eigen::Vector2d du2dy;
    gradLogEuclidianP2Interpolation(xsi, elm, _logMetricTensorAtNodetags_eigen[GMSHNODETAGS_P2[0]-1],
                                              _logMetricTensorAtNodetags_eigen[GMSHNODETAGS_P2[1]-1],
                                              _logMetricTensorAtNodetags_eigen[GMSHNODETAGS_P2[2]-1],
                                              _logMetricTensorAtNodetags_eigen[GMSHNODETAGS_P2[3]-1],
                                              _logMetricTensorAtNodetags_eigen[GMSHNODETAGS_P2[4]-1],
                                              _logMetricTensorAtNodetags_eigen[GMSHNODETAGS_P2[5]-1], dMdx, dMdy, L, DLDX, DLDY,
                                              L1, DL1DX, DL1DY, L2, DL2DX, DL2DY, U1, DU1DX, DU1DY, U2, DU2DX, DU2DY);
  }
}

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
                           Eigen::Vector2d &du2dy)
{
  static_cast<feMetric*>(metricPtr)->interpolateMetricP2Log(x, M, dMdx, dMdy, L, dLdx, dLdy, l1, dl1dx, dl1dy, l2, dl2dx, dl2dy, u1, du1dx, du1dy, u2, du2dx, du2dy);
}
void feMetric::interpolateMetricP2Log(const double *x, Eigen::Matrix2d &M,
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
                                   Eigen::Vector2d &du2dy)
{
  // Locate point in the feMesh
  int elm;
  double UVW[3];
  bool isFound = static_cast<feMesh2DP1 *>(_recoveredFields[0]->_mesh)->locateVertex(x, elm, UVW);
  if(!isFound) {
    feWarning("In interpolateMetricP2 : Point (%f, %f) was not found in the mesh.", x[0], x[1]);
    M = Eigen::Matrix2d::Identity();
    dMdx = Eigen::Matrix2d::Zero();
    dMdy = Eigen::Matrix2d::Zero();
  } else {

    double xsi[2] = {UVW[0], UVW[1]};

    // Get the vertex tags from the feMesh
    int gmshNodeTags[6];
    for(int i = 0; i < 6; ++i) {
      int vertexTag = _recoveredFields[0]->_cnc->getVertexConnectivity(elm, i);
      gmshNodeTags[i] = _sequentialTag2nodeTag[vertexTag];
    }

    Eigen::Matrix2d res;
    this->logEuclidianP2Interpolation(xsi, _logMetricTensorAtNodetags_eigen[gmshNodeTags[0]-1],
                                           _logMetricTensorAtNodetags_eigen[gmshNodeTags[1]-1],
                                           _logMetricTensorAtNodetags_eigen[gmshNodeTags[2]-1],
                                           _logMetricTensorAtNodetags_eigen[gmshNodeTags[3]-1],
                                           _logMetricTensorAtNodetags_eigen[gmshNodeTags[4]-1],
                                           _logMetricTensorAtNodetags_eigen[gmshNodeTags[5]-1], res);
    M(0,0) = res(0,0);
    M(0,1) = res(0,1);
    M(1,0) = res(1,0);
    M(1,1) = res(1,1);

    gradLogEuclidianP2Interpolation(xsi, elm, _logMetricTensorAtNodetags_eigen[gmshNodeTags[0]-1],
                                              _logMetricTensorAtNodetags_eigen[gmshNodeTags[1]-1],
                                              _logMetricTensorAtNodetags_eigen[gmshNodeTags[2]-1],
                                              _logMetricTensorAtNodetags_eigen[gmshNodeTags[3]-1],
                                              _logMetricTensorAtNodetags_eigen[gmshNodeTags[4]-1],
                                              _logMetricTensorAtNodetags_eigen[gmshNodeTags[5]-1], dMdx, dMdy, L, dLdx, dLdy,
                                              l1, dl1dx, dl1dy, l2, dl2dx, dl2dy, u1, du1dx, du1dy, u2, du2dx, du2dy);
  }
}

void getPhiTriP2(double r, double s, double phi[6], double dphidr[6], double dphids[6]){
  phi[0] = (1. - r - s) * (1. - 2.*r - 2.*s);
  phi[1] = r * (2. * r - 1.);
  phi[2] = s * (2. * s - 1.);
  phi[3] = 4. * r * (1. - r - s);
  phi[4] = 4. * r * s;
  phi[5] = 4. * s * (1. - r - s);

  dphidr[0] = 4. * (r+s) - 3.;
  dphidr[1] = 4. * r - 1.;
  dphidr[2] = 0.;
  dphidr[3] = 4.*(1.-2.*r-s);
  dphidr[4] = 4.*s;
  dphidr[5] = -4.*s;

  dphids[0] = 4.*(r+s)-3.;
  dphids[1] = 0.;
  dphids[2] = 4.*s-1.;
  dphids[3] = -4.*r;
  dphids[4] = 4.*r;
  dphids[5] = 4.*(1.-r-2.*s);
}

Eigen::Matrix2d IDENTITY = Eigen::Matrix2d::Identity();

// Gradient of the interpolated log-metric on P1 background mesh
// with respect to front mesh coordinates
void feMetric::gradLogEuclidianP2Interpolation(const double xsi[2],
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
                                               Eigen::Vector2d &du2dyres)
{
  double phi[6], dphidr[6], dphids[6];
  double r = xsi[0];
  double s = xsi[1];

  getPhiTriP2(r, s, phi, dphidr, dphids);

  // Step 1: compute inverse of grad_(xsiHat) of phi_pHat
  // Inverse of jacobian of P1 transformation in the background mesh
  Triangle *t = _recoveredFields[0]->_mesh->_elements[element];
  double x0 = t->getVertex(0)->x();
  double y0 = t->getVertex(0)->y();
  double x1 = t->getVertex(1)->x();
  double y1 = t->getVertex(1)->y();
  double x2 = t->getVertex(2)->x();
  double y2 = t->getVertex(2)->y();

  // FIXME: TRANSFORMATION SHOULD BE QUADRATIC!
  Eigen::Matrix2d Jac;
  Jac(0,0) = x1 - x0;
  Jac(0,1) = x2 - x0;
  Jac(1,0) = y1 - y0;
  Jac(1,1) = y2 - y0;
  Eigen::Matrix2d invJ = Jac.inverse();

  // Step 2: compute gradient of BGM shape functions with respect to 
  // the coordinates in the front mesh. This is dphiHat/dxsiHat * dxsiHat/dx,
  // where hat are the background mesh quantities.

  // Hardcoded gradient of the 6 P2 shape functions with respect to xj, yj in front mesh
  Eigen::RowVector2d dN0dx, dN1dx, dN2dx, dN3dx, dN4dx, dN5dx;
  dN0dx << dphidr[0], dphids[0]; dN0dx = dN0dx * invJ;
  dN1dx << dphidr[1], dphids[1]; dN1dx = dN1dx * invJ;
  dN2dx << dphidr[2], dphids[2]; dN2dx = dN2dx * invJ;
  dN3dx << dphidr[3], dphids[3]; dN3dx = dN3dx * invJ;
  dN4dx << dphidr[4], dphids[4]; dN4dx = dN4dx * invJ;
  dN5dx << dphidr[5], dphids[5]; dN5dx = dN5dx * invJ;

  Eigen::Matrix2d L    = logM0 * phi[0] + logM1 * phi[1] + logM2 * phi[2]
                       + logM3 * phi[3] + logM4 * phi[4] + logM5 * phi[5];
  Eigen::Matrix2d dLdx = logM0 * dN0dx(0) + logM1 * dN1dx(0) + logM2 * dN2dx(0)
                       + logM3 * dN3dx(0) + logM4 * dN4dx(0) + logM5 * dN5dx(0);
  Eigen::Matrix2d dLdy = logM0 * dN0dx(1) + logM1 * dN1dx(1) + logM2 * dN2dx(1)
                       + logM3 * dN3dx(1) + logM4 * dN4dx(1) + logM5 * dN5dx(1);

  // Step 3: compute gradient of eigenvalues and eigenvectors of L　from Aparicio-Estrems et al, High order
  // metric interpolation for curved r-adaptation by distortion minimization)

  Eigen::EigenSolver<Eigen::Matrix2d> es;
  Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType eigenvalues;
  Eigen::EigenSolver<Eigen::Matrix2d>::EigenvectorsType eigenvectors;

  es.compute(L, true);
  eigenvalues = es.eigenvalues();
  Eigen::Vector2d lambda = eigenvalues.real();

  eigenvectors = es.eigenvectors();
  Eigen::Matrix2d U = eigenvectors.real();
  Eigen::Vector2d u1 = U.col(0);
  Eigen::Vector2d u2 = U.col(1);

  // Gradient des valeurs propres d_j lambda_l
  double dl1dx = u1.transpose() * dLdx * u1;
  double dl1dy = u1.transpose() * dLdy * u1;

  double dl2dx = u2.transpose() * dLdx * u2;
  double dl2dy = u2.transpose() * dLdy * u2;

  l1res = lambda(0);
  l2res = lambda(1);
  dl1dxres = dl1dx;
  dl1dyres = dl1dy;
  dl2dxres = dl2dx;
  dl2dyres = dl2dy;

  // // Gradient des vecteurs propres d_j u_l
  // // Solve one at the time to avoid creating a lot of Eigen::Matrix2d
  // Eigen::Matrix2d L1 = L - lambda(0) * IDENTITY;
  // Eigen::Matrix2d dL1dx = dLdx - dl1dx * IDENTITY;
  // Eigen::Matrix2d dL1dy = dLdy - dl1dy * IDENTITY;

  // Eigen::Matrix2d L2 = L - lambda(1) * IDENTITY;
  // Eigen::Matrix2d dL2dx = dLdx - dl2dx * IDENTITY;
  // Eigen::Matrix2d dL2dy = dLdy - dl2dy * IDENTITY;

  // const Eigen::CompleteOrthogonalDecomposition<Eigen::Matrix2d> L1_decomposition = L1.completeOrthogonalDecomposition();
  // Eigen::Vector2d du1dx = L1_decomposition.solve(-dL1dx * u1);
  // Eigen::Vector2d du1dy = L1_decomposition.solve(-dL1dy * u1);

  // // Correction de https://math.stackexchange.com/questions/2689374/derivative-of-eigenvectors-of-a-symmetric-matrix-valued-function
  // Eigen::Vector2d w = du1dx;
  // Eigen::Vector2d v = u1;
  // du1dx = w - (v.transpose() * w) * v;
  // // feInfo("du1dx ref = %+-1.10e - %+-1.10e", du1dx(0), du1dx(1));

  // w = du1dy;
  // du1dy = w - (v.transpose() * w) * v;
  // // feInfo("du1dy ref = %+-1.10e - %+-1.10e", du1dy(0), du1dy(1));

  // const Eigen::CompleteOrthogonalDecomposition<Eigen::Matrix2d> L2_decomposition = L2.completeOrthogonalDecomposition();
  // Eigen::Vector2d du2dx = L2_decomposition.solve(-dL2dx * u2);
  // Eigen::Vector2d du2dy = L2_decomposition.solve(-dL2dy * u2);

  // w = du2dx;
  // v = u2;
  // du2dx = w - (v.transpose() * w) * v;
  // // feInfo("du2dx ref = %+-1.10e - %+-1.10e", du2dx(0), du2dx(1));

  // w = du2dy;
  // du2dy = w - (v.transpose() * w) * v;
  // feInfo("du2dy ref = %+-1.10e - %+-1.10e", du2dy(0), du2dy(1));

  // u1res = u1;
  // u2res = u2;
  // du1dxres = du1dx;
  // du1dyres = du1dy;
  // du2dxres = du2dx;
  // du2dyres = du2dy;

  // Eigen::Matrix2d dUdx = Eigen::Matrix2d::Zero(), dUdy = Eigen::Matrix2d::Zero();
  // dUdx.col(0) = du1dx;
  // dUdx.col(1) = du2dx;
  // dUdy.col(0) = du1dy;
  // dUdy.col(1) = du2dy;

  ////////////////////////////////////////////////////////
  // Gradient des vecteurs propres d_j u_l
  // Solve one at the time to avoid creating a lot of Eigen::Matrix2d
  Eigen::Matrix2d dUdx, dUdy;

  Eigen::Matrix2d L_l;
  L_l(0,0) = L(0,0) - lambda(0);
  L_l(0,1) = L(0,1);
  L_l(1,0) = L(1,0);
  L_l(1,1) = L(1,1) - lambda(0);
  const Eigen::CompleteOrthogonalDecomposition<Eigen::Matrix2d> L_l_decomposition = L_l.completeOrthogonalDecomposition();

  // // Solve for du1dx - gradL is dL1dx - gradU is du1dx
  Eigen::Matrix2d gradL; // = dLdx - dl1dx * IDENTITY;
  gradL(0,0) = dLdx(0,0) - dl1dx;
  gradL(0,1) = dLdx(0,1);
  gradL(1,0) = dLdx(1,0);
  gradL(1,1) = dLdx(1,1) - dl1dx;
  Eigen::Vector2d gradU = L_l_decomposition.solve(-gradL * u1);
  // // Correction de https://math.stackexchange.com/questions/2689374/derivative-of-eigenvectors-of-a-symmetric-matrix-valued-function
  Eigen::Vector2d w = gradU;
  Eigen::Vector2d v = u1;
  gradU = w - (v.transpose() * w) * v;
  dUdx.col(0) = gradU;
  du1dxres = gradU;

  // Solve for du1dy - gradL is dL1dy - gradU is du1dy
  // gradL = dLdy - dl1dy * IDENTITY;
  gradL(0,0) = dLdy(0,0) - dl1dy;
  gradL(0,1) = dLdy(0,1);
  gradL(1,0) = dLdy(1,0);
  gradL(1,1) = dLdy(1,1) - dl1dy;
  gradU = L_l_decomposition.solve(-gradL * u1);
  w = gradU;
  gradU = w - (v.transpose() * w) * v;
  dUdy.col(0) = gradU;
  du1dyres = gradU;

  L_l = L - lambda(1) * IDENTITY;
  const Eigen::CompleteOrthogonalDecomposition<Eigen::Matrix2d> L_l_decomposition2 = L_l.completeOrthogonalDecomposition();

  // Solve for du2dx - gradL is dL2dx - gradU is du2dx
  // gradL = dLdx - dl2dx * IDENTITY;
  gradL(0,0) = dLdx(0,0) - dl2dx;
  gradL(0,1) = dLdx(0,1);
  gradL(1,0) = dLdx(1,0);
  gradL(1,1) = dLdx(1,1) - dl2dx;
  gradU = L_l_decomposition2.solve(-gradL * u2);
  w = gradU;
  v = u2;
  gradU = w - (v.transpose() * w) * v;
  dUdx.col(1) = gradU;
  du2dxres = gradU;

  // Solve for du2dy - gradL is dL2dy - gradU is du2dy
  // gradL = dLdy - dl2dy * IDENTITY;
  gradL(0,0) = dLdy(0,0) - dl2dy;
  gradL(0,1) = dLdy(0,1);
  gradL(1,0) = dLdy(1,0);
  gradL(1,1) = dLdy(1,1) - dl2dy;
  gradU = L_l_decomposition2.solve(-gradL * u2);
  w = gradU;
  gradU = w - (v.transpose() * w) * v;
  dUdy.col(1) = gradU;
  du2dyres = gradU;

  u1res = u1;
  u2res = u2;
  ////////////////////////////////////////////////////////

  // Eigen::Matrix2d expD;
  // expD << exp(lambda(0)), 0., 0., exp(lambda(1));
  // Eigen::Matrix2d dDdx;
  // dDdx << dl1dx, 0., 0., dl2dx;
  // Eigen::Matrix2d dDdy;
  // dDdy << dl1dy, 0., 0., dl2dy;

  // // std::cout << expD <<std::endl;
  // // std::cout << dDdx <<std::endl;

  // Eigen::Matrix2d dexpDdx = expD * dDdx;
  // Eigen::Matrix2d dexpDdy = expD * dDdy;

  // dMdx = dUdx * expD * U.transpose() + U * dexpDdx * U.transpose() + U * expD * dUdx.transpose();
  // dMdy = dUdy * expD * U.transpose() + U * dexpDdy * U.transpose() + U * expD * dUdy.transpose();

  Eigen::Matrix2d expD;
  expD << exp(lambda(0)), 0., 0., exp(lambda(1));

  Eigen::Matrix2d gradExp;
  gradExp << dl1dx, 0., 0., dl2dx;
  gradExp *= expD;
  dMdx = dUdx * expD * U.transpose() + U * gradExp * U.transpose() + U * expD * dUdx.transpose();

  gradExp << dl1dy, 0., 0., dl2dy;
  gradExp *= expD;
  dMdy = dUdy * expD * U.transpose() + U * gradExp * U.transpose() + U * expD * dUdy.transpose();
}

// void feMetric::interpolateMetricP1(const double *x, Eigen::Matrix2d &M, Eigen::Matrix2d &sumduda1M,
//                                    Eigen::Matrix2d &sumduda2M)
// {
//   int elm;
//   bool isFound = static_cast<feMesh2DP1 *>(_recoveredFields[0]->_mesh)->locateVertex(x, elm, UVW);
//   if(!isFound) {
//     feWarning("In interpolateMetricP1 : Point (%f, %f) was not found in the mesh.\n", x[0], x[1]);
//     M = Eigen::Matrix2d::Identity();
//     return;
//   } else {
//     // Interpolate the log-metric and take exponential
//     Triangle *t = _recoveredFields[0]->_mesh->_elements[elm];

//     int t0 = _v2n[t->getVertex(0)];
//     int t1 = _v2n[t->getVertex(1)];
//     int t2 = _v2n[t->getVertex(2)];

//     std::map<int, Eigen::Matrix2d>::iterator it = _metricsOnGmshModel_eigen.find(t0);
//     if(it == _metricsOnGmshModel_eigen.end()) {
//       feInfo("Metric not found at %d", t0);
//     }
//     it = _metricsOnGmshModel_eigen.find(t1);
//     if(it == _metricsOnGmshModel_eigen.end()) {
//       feInfo("Metric not found at %d", t1);
//     }
//     it = _metricsOnGmshModel_eigen.find(t2);
//     if(it == _metricsOnGmshModel_eigen.end()) {
//       feInfo("Metric not found at %d", t2);
//     }

//     Eigen::Matrix2d M0 = _metricsOnGmshModel_eigen[t0];
//     Eigen::Matrix2d M1 = _metricsOnGmshModel_eigen[t1];
//     Eigen::Matrix2d M2 = _metricsOnGmshModel_eigen[t2];

//     Eigen::Matrix2d logM0 = M0.log();
//     Eigen::Matrix2d logM1 = M1.log();
//     Eigen::Matrix2d logM2 = M2.log();

//     M = (1. - UVW[0] - UVW[1]) * logM0 + UVW[0] * logM1 + UVW[1] * logM2;
//     M = M.exp();

//      // The barycentric coordinates (u,v,w) are different from UVW above :

//      //    u = 1 - r - s
//      //    v = r
//      //    w = s

//      //    with UVW[0] = r, UVW[1] = s, UVW[2] = t = 0 in 2D
//      //    Unfortunate notation (-:
    

//     // Gradient of the barycentric coordinates w.r.t. P11^l
//     double xa = t->getVertex(0)->x();
//     double ya = t->getVertex(0)->y();
//     double xb = t->getVertex(1)->x();
//     double yb = t->getVertex(1)->y();
//     double xc = t->getVertex(2)->x();
//     double yc = t->getVertex(2)->y();

//     // Area of the triangle, the 1/2 is be cancelled by the area of the subtriangles
//     // associated to the barycentric coordinates. Should we use fabs ?
//     double A = (xa - xc) * (yb - yc) - (xb - xc) * (ya - yc);

//     GRADUVW[0] = (yb - yc) / A; // grad u = (duda1, duda2)
//     GRADUVW[1] = -(xb - xc) / A;
//     GRADUVW[2] = (yc - ya) / A; // grad v
//     GRADUVW[3] = -(xc - xa) / A;
//     GRADUVW[4] = (ya - yb) / A; // grad w
//     GRADUVW[5] = -(xa - xb) / A;

//     sumduda1M = GRADUVW[0] * logM0 + GRADUVW[2] * logM1 + GRADUVW[4] * logM2;
//     sumduda2M = GRADUVW[1] * logM0 + GRADUVW[3] * logM1 + GRADUVW[5] * logM2;
//   }
// }

// void feMetric::interpolateMetricP1(const double *x, const double *gammaOrth, Eigen::Matrix2d &M,
//                                    Eigen::Matrix2d &sumdudaM)
// {
//   int elm;
//   bool isFound = static_cast<feMesh2DP1 *>(_recovery->_mesh)->locateVertex(x, elm, UVW);
//   if(!isFound) {
//     feWarning("In interpolateMetricP1 : Point (%f, %f) was not found in the mesh.\n", x[0], x[1]);
//     M = Eigen::Matrix2d::Identity();
//     return;
//   } else {
//     // Interpolate the log-metric and take exponential
//     Triangle *t = _recovery->_mesh->_elements[elm];

//     int t0 = _v2n[t->getVertex(0)];
//     int t1 = _v2n[t->getVertex(1)];
//     int t2 = _v2n[t->getVertex(2)];

//     Eigen::Matrix2d M0 = _metricsOnGmshModel_eigen[t0];
//     Eigen::Matrix2d M1 = _metricsOnGmshModel_eigen[t1];
//     Eigen::Matrix2d M2 = _metricsOnGmshModel_eigen[t2];

//     Eigen::Matrix2d logM0 = M0.log();
//     Eigen::Matrix2d logM1 = M1.log();
//     Eigen::Matrix2d logM2 = M2.log();

//     M = (1. - UVW[0] - UVW[1]) * logM0 + UVW[0] * logM1 + UVW[1] * logM2;
//     M = M.exp();

//     // Gradient of the barycentric coordinates w.r.t. P11^l
//     double xa = t->getVertex(0)->x();
//     double ya = t->getVertex(0)->y();
//     double xb = t->getVertex(1)->x();
//     double yb = t->getVertex(1)->y();
//     double xc = t->getVertex(2)->x();
//     double yc = t->getVertex(2)->y();
//     // Area of the triangle, the 1/2 is be cancelled by the area of the subtriangles
//     // associated to the barycentric coordinates. Should we use fabs ?
//     double A = (xa - xc) * (yb - yc) - (xb - xc) * (ya - yc);
//     // Derivative with respect to alpha in 1D, along the bisector
//     double gx = gammaOrth[0];
//     double gy = gammaOrth[1];
//     double duda = (gx * (yb - yc) - gy * (xb - xc)) / A;
//     double dvda = (gx * (yc - ya) - gy * (xc - xa)) / A;
//     double dwda = (gx * (ya - yb) - gy * (xa - xb)) / A;

//     sumdudaM = duda * logM0 + dvda * logM1 + dwda * logM2;
//   }
// }

// void feMetric::interpolateMetricAndDerivativeOnP2Edge(
//   double t, const Eigen::Matrix2d &M11, const Eigen::Matrix2d &M20, const Eigen::Matrix2d &M02,
//   const Eigen::Matrix2d &sumduda1M, const Eigen::Matrix2d &sumduda2M, Eigen::Matrix2d &M,
//   Eigen::Matrix2d &dMda1, Eigen::Matrix2d &dMda2)
// {
//   double phi20 = (1. - t) * (1. - 2. * t);
//   double phi02 = t * (2. * t - 1.);
//   double L = 4. * t * (1. - t); // phi11

//   Eigen::Matrix2d T = phi20 * M20.log() + phi02 * M02.log() + L * M11.log(); // = log M(x(t))

//   // lT = T;

//   // // Scaling and squaring : find m s.t. ||T||/(2^m) <= 1
//   // int m = 0;
//   // if(lT.norm() > 1){
//   //   m = 1;
//   //   while(true){
//   //     lT /= pow(2,m);
//   //     if(lT.norm() <= 1.)
//   //       break;
//   //     m++;
//   //   }
//   //   // feInfo("Found m = %d", m);
//   //   // feInfo("norm = %f", lT.norm());
//   //   // std::cout<<lT<<std::endl;
//   // }

//   Eigen::Matrix2d Tn = T;
//   M = T.exp();

//   Eigen::Matrix2d dTda1 = L * sumduda1M;
//   Eigen::Matrix2d dTda2 = L * sumduda2M;
//   Eigen::Matrix2d dTnda1 = dTda1;
//   Eigen::Matrix2d dTnda2 = dTda2;

//   // Compute derivative of exponential by recursion
//   // Term for k = 0 is the zero matrix (derivative of I)
//   // Term for k = 1 is dTdai
//   int kfac = 1;
//   dMda1 = dTda1;
//   dMda2 = dTda2;

//   for(int k = 2; k < 10; ++k) {
//     dTnda1 = Tn * dTda1 + dTnda1 * T;
//     dTnda2 = Tn * dTda2 + dTnda2 * T;

//     Tn *= T; // T^n

//     kfac *= k;
//     dMda1 += 1. / ((double)kfac) * dTnda1;
//     dMda2 += 1. / ((double)kfac) * dTnda2;
//   }

//   // // Compute exponential and derivative of lT instead
//   // double lambda = pow(2., -m);
//   // Tn = lT;
//   // M = lT.exp();
//   // M = M.pow(1./lambda);

//   // dTda1 = L * sumduda1M;
//   // dTda2 = L * sumduda2M;
//   // dTnda1 = dTda1;
//   // dTnda2 = dTda2;

//   // // Compute derivative of exponential by recursion
//   // // Term for k = 0 is the zero matrix (derivative of I)
//   // // Term for k = 1 is dTdai
//   // int kfac = 1;
//   // dMda1 = dTda1;
//   // dMda2 = dTda2;

//   // for(int k = 2; k < 10; ++k){
//   //   dTnda1 = Tn * dTda1 + dTnda1 * lT;
//   //   dTnda2 = Tn * dTda2 + dTnda2 * lT;

//   //   Tn *= lT; // T^n

//   //   kfac *= k;
//   //   dMda1 += 1./(double) kfac * dTnda1;
//   //   dMda2 += 1./(double) kfac * dTnda2;
//   // }

//   // Amat = M.pow(lambda - 1.);
//   // Amat = Amat.inverse();
//   // Amat *= 1./lambda;

//   // dMda1 = Amat * dMda1;
//   // dMda2 = Amat * dMda2;
// }

// void feMetric::interpolateMetricAndDerivativeOnP2Edge(double t, const Eigen::Matrix2d &M11,
//                                                       const Eigen::Matrix2d &M20,
//                                                       const Eigen::Matrix2d &M02,
//                                                       const Eigen::Matrix2d &sumdudaM,
//                                                       Eigen::Matrix2d &M, Eigen::Matrix2d &dMda)
// {
//   double phi20 = (1. - t) * (1. - 2. * t);
//   double phi02 = t * (2. * t - 1.);
//   double L = 4. * t * (1. - t); // phi11

//   Eigen::Matrix2d T = phi20 * M20.log() + phi02 * M02.log() + L * M11.log(); // = log M(x(t))

//   Eigen::Matrix2d Tn = T;
//   M = T.exp();

//   Eigen::Matrix2d dTda = L * sumdudaM;
//   Eigen::Matrix2d dTnda = dTda;

//   // Compute derivative of exponential by recursion
//   // Term for k = 0 is the zero matrix (derivative of I)
//   // Term for k = 1 is dTdai
//   int kfac = 1;
//   dMda = dTda;

//   for(int k = 2; k < 10; ++k) {
//     dTnda = Tn * dTda + dTnda * T;
//     Tn *= T; // T^n
//     kfac *= k;
//     dMda += 1. / ((double)kfac) * dTnda;
//   }
// }

// void interpolateMetricP1WithDerivativesWrapper(void *metric, const double *x, Eigen::Matrix2d &M,
//                                                Eigen::Matrix2d &dMdx, Eigen::Matrix2d &dMdy)
// {
//   feMetric *m = static_cast<feMetric *>(metric);
//   m->interpolateMetricP1WithDerivatives(x, M, dMdx, dMdy);
// }

// void interpolateMetricP1Wrapper(void *metric, const double *x, Eigen::Matrix2d &M,
//                                 Eigen::Matrix2d &sumduda1M, Eigen::Matrix2d &sumduda2M)
// {
//   feMetric *m = static_cast<feMetric *>(metric);
//   m->interpolateMetricP1(x, M, sumduda1M, sumduda2M);
// }

// void interpolateMetricP1Wrapper1D(void *metric, const double *x, const double *gammaOrth,
//                                   Eigen::Matrix2d &M, Eigen::Matrix2d &sumdudaM)
// {
//   feMetric *m = static_cast<feMetric *>(metric);
//   m->interpolateMetricP1(x, gammaOrth, M, sumdudaM);
// }

// void interpolateMetricAndDerivativeOnP2EdgeWrapper(
//   void *metric, const double t, const Eigen::Matrix2d &M11, const Eigen::Matrix2d &M20,
//   const Eigen::Matrix2d &M02, const Eigen::Matrix2d &sumduda1M, const Eigen::Matrix2d &sumduda2M,
//   Eigen::Matrix2d &M, Eigen::Matrix2d &dMda1, Eigen::Matrix2d &dMda2)
// {
//   feMetric *m = static_cast<feMetric *>(metric);
//   m->interpolateMetricAndDerivativeOnP2Edge(t, M11, M20, M02, sumduda1M, sumduda2M, M, dMda1,
//                                             dMda2);
// }

// void interpolateMetricAndDerivativeOnP2EdgeWrapper1D(void *metric, const double t,
//                                                      const Eigen::Matrix2d &M11,
//                                                      const Eigen::Matrix2d &M20,
//                                                      const Eigen::Matrix2d &M02,
//                                                      const Eigen::Matrix2d &sumdudaM,
//                                                      Eigen::Matrix2d &M, Eigen::Matrix2d &dMda)
// {
//   feMetric *m = static_cast<feMetric *>(metric);
//   m->interpolateMetricAndDerivativeOnP2Edge(t, M11, M20, M02, sumdudaM, M, dMda);
// }