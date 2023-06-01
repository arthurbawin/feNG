#include "feMetric.h"

// #include "../contrib/Eigen/Eigen"
// #include "../contrib/Eigen/Eigenvalues"
// #include "../contrib/Eigen/QR"
#include "../contrib/unsupported/Eigen/MatrixFunctions"

static double UVW[3];
static double GRADUVW[6];

void feMetric::interpolateMetricP1WithDerivatives(const double *x, Eigen::Matrix2d &M,
                                                  Eigen::Matrix2d &dMdx, Eigen::Matrix2d &dMdy)
{
  feErrorMsg(FE_STATUS_ERROR, " SHOULD NOT USE THIS INTERPOLATION. USE LOGEUCLIDIAN INTERPOLATION ON NEW BACKGROUND MESH INSTEAD");
  exit(-1);

  int elm;
  bool isFound = static_cast<feMesh2DP1 *>(_newRecovery->_mesh)->locateVertex(x, elm, UVW);
  if(!isFound) {
    feWarning("In interpolateMetricP1 : Point (%f, %f) was not found in the mesh.\n", x[0], x[1]);
    M = Eigen::Matrix2d::Identity();
    return;
  } else {
    double r = UVW[0];
    double s = UVW[1];

    // Interpolate the log-metric and take exponential
    Triangle *t = _newRecovery->_mesh->_elements[elm];

    int t0 = _v2n[t->getVertex(0)];
    int t1 = _v2n[t->getVertex(1)];
    int t2 = _v2n[t->getVertex(2)];

    Eigen::Matrix2d M0 = _metricsOnGmshModel_eigen[t0];
    Eigen::Matrix2d M1 = _metricsOnGmshModel_eigen[t1];
    Eigen::Matrix2d M2 = _metricsOnGmshModel_eigen[t2];

    // std::cout << "M0 (at " << t0 << ") = " << M0 << std::endl;
    // std::cout << "M1 (at " << t1 << ") = " << M1 << std::endl;
    // std::cout << "M2 (at " << t2 << ") = " << M2 << std::endl;

    // std::cout << "r = " << r <<std::endl;
    // std::cout << "s = " << s <<std::endl;

    Eigen::Matrix2d logM0 = M0.log();
    Eigen::Matrix2d logM1 = M1.log();
    Eigen::Matrix2d logM2 = M2.log();

    Eigen::Matrix2d L = (1. - r - s) * logM0 + r * logM1 + s * logM2;
    M = L.exp();

    // std::cout << "L = " << L <<std::endl;
    // std::cout << "M = " << M <<std::endl;

    // Compute the derivatives of the eigen decomposition (from Aparicio-Estrems et al, High order
    // metric interpolation for curved r-adaptation by distortion minimization, IMR 2021 or 2022)
    Eigen::EigenSolver<Eigen::Matrix2d> es;
    Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType eigenvalues;
    Eigen::EigenSolver<Eigen::Matrix2d>::EigenvectorsType eigenvectors;

    es.compute(L, true);
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

    // Gradient de L ( = sum_i log M_i phi_i)
    double x1 = t->getVertex(0)->x();
    double y1 = t->getVertex(0)->y();
    double x2 = t->getVertex(1)->x();
    double y2 = t->getVertex(1)->y();
    double x3 = t->getVertex(2)->x();
    double y3 = t->getVertex(2)->y();

    // std::cout<<"Coordonnees"<<std::endl;
    // std::cout<<x1<<std::endl;
    // std::cout<<y1<<std::endl;
    // std::cout<<x2<<std::endl;
    // std::cout<<y2<<std::endl;
    // std::cout<<x3<<std::endl;
    // std::cout<<y3<<std::endl;

    double J = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);

    // std::cout << "J = " << J <<std::endl;

    double drdx = (y3 - y1) / J;
    double drdy = -(x3 - x1) / J;
    double dsdx = -(y2 - y1) / J;
    double dsdy = (x2 - x1) / J;

    double dphi0dr = -1.;
    double dphi0ds = -1.;
    double dphi1dr = 1.;
    double dphi1ds = 0.;
    double dphi2dr = 0.;
    double dphi2ds = 1.;

    double dphi0dx = dphi0dr * drdx + dphi0ds * dsdx;
    double dphi0dy = dphi0dr * drdy + dphi0ds * dsdy;
    double dphi1dx = dphi1dr * drdx + dphi1ds * dsdx;
    double dphi1dy = dphi1dr * drdy + dphi1ds * dsdy;
    double dphi2dx = dphi2dr * drdx + dphi2ds * dsdx;
    double dphi2dy = dphi2dr * drdy + dphi2ds * dsdy;

    Eigen::Matrix2d dLdx = logM0 * dphi0dx + logM1 * dphi1dx + logM2 * dphi2dx;
    Eigen::Matrix2d dLdy = logM0 * dphi0dy + logM1 * dphi1dy + logM2 * dphi2dy;

    // std::cout << "dLdx = " << dLdx <<std::endl;
    // std::cout << "dLdy = " << dLdy <<std::endl;

    // Gradient des valeurs propres d_j lambda_l
    double dl1dx = u1.transpose() * dLdx * u1;
    double dl1dy = u1.transpose() * dLdy * u1;

    double dl2dx = u2.transpose() * dLdx * u2;
    double dl2dy = u2.transpose() * dLdy * u2;

    // std::cout << "dl1dx = " << dl1dx <<std::endl;
    // std::cout << "dl2dx = " << dl2dx <<std::endl;
    // std::cout << "dl1dy = " << dl1dy <<std::endl;
    // std::cout << "dl2dy = " << dl2dy <<std::endl;

    // Gradient des vecteurs propres d_j u_l
    Eigen::Matrix2d L1 = L - lambda(0) * Eigen::Matrix2d::Identity();
    Eigen::Matrix2d dL1dx = dLdx - dl1dx * Eigen::Matrix2d::Identity();
    Eigen::Matrix2d dL1dy = dLdy - dl1dy * Eigen::Matrix2d::Identity();

    Eigen::Matrix2d L2 = L - lambda(1) * Eigen::Matrix2d::Identity();
    Eigen::Matrix2d dL2dx = dLdx - dl2dx * Eigen::Matrix2d::Identity();
    Eigen::Matrix2d dL2dy = dLdy - dl2dy * Eigen::Matrix2d::Identity();

    Eigen::Vector2d du1dx = L1.completeOrthogonalDecomposition().solve(-dL1dx * u1);
    Eigen::Vector2d du1dy = L1.completeOrthogonalDecomposition().solve(-dL1dy * u1);

    Eigen::Vector2d du2dx = L2.completeOrthogonalDecomposition().solve(-dL2dx * u2);
    Eigen::Vector2d du2dy = L2.completeOrthogonalDecomposition().solve(-dL2dy * u2);

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

    // std::cout << dMdx <<std::endl;
    // std::cout << dMdy <<std::endl;
    // std::cout<< "ref = " << 2.123 * 2. * x[0] << std::endl;

    // exit(-1);
  }
}

void feMetric::interpolateMetricP1(const double *x, Eigen::Matrix2d &M, Eigen::Matrix2d &sumduda1M,
                                   Eigen::Matrix2d &sumduda2M)
{
  int elm;
  bool isFound = static_cast<feMesh2DP1 *>(_newRecovery->_mesh)->locateVertex(x, elm, UVW);
  if(!isFound) {
    feWarning("In interpolateMetricP1 : Point (%f, %f) was not found in the mesh.\n", x[0], x[1]);
    M = Eigen::Matrix2d::Identity();
    return;
  } else {
    // Interpolate the log-metric and take exponential
    Triangle *t = _newRecovery->_mesh->_elements[elm];

    int t0 = _v2n[t->getVertex(0)];
    int t1 = _v2n[t->getVertex(1)];
    int t2 = _v2n[t->getVertex(2)];

    std::map<int, Eigen::Matrix2d>::iterator it = _metricsOnGmshModel_eigen.find(t0);
    if(it == _metricsOnGmshModel_eigen.end()) {
      feInfo("Metric not found at %d", t0);
    }
    it = _metricsOnGmshModel_eigen.find(t1);
    if(it == _metricsOnGmshModel_eigen.end()) {
      feInfo("Metric not found at %d", t1);
    }
    it = _metricsOnGmshModel_eigen.find(t2);
    if(it == _metricsOnGmshModel_eigen.end()) {
      feInfo("Metric not found at %d", t2);
    }

    Eigen::Matrix2d M0 = _metricsOnGmshModel_eigen[t0];
    Eigen::Matrix2d M1 = _metricsOnGmshModel_eigen[t1];
    Eigen::Matrix2d M2 = _metricsOnGmshModel_eigen[t2];

    Eigen::Matrix2d logM0 = M0.log();
    Eigen::Matrix2d logM1 = M1.log();
    Eigen::Matrix2d logM2 = M2.log();

    M = (1. - UVW[0] - UVW[1]) * logM0 + UVW[0] * logM1 + UVW[1] * logM2;
    M = M.exp();

     // The barycentric coordinates (u,v,w) are different from UVW above :

     //    u = 1 - r - s
     //    v = r
     //    w = s

     //    with UVW[0] = r, UVW[1] = s, UVW[2] = t = 0 in 2D
     //    Unfortunate notation (-:
    

    // Gradient of the barycentric coordinates w.r.t. P11^l
    double xa = t->getVertex(0)->x();
    double ya = t->getVertex(0)->y();
    double xb = t->getVertex(1)->x();
    double yb = t->getVertex(1)->y();
    double xc = t->getVertex(2)->x();
    double yc = t->getVertex(2)->y();

    // Area of the triangle, the 1/2 is be cancelled by the area of the subtriangles
    // associated to the barycentric coordinates. Should we use fabs ?
    double A = (xa - xc) * (yb - yc) - (xb - xc) * (ya - yc);

    GRADUVW[0] = (yb - yc) / A; // grad u = (duda1, duda2)
    GRADUVW[1] = -(xb - xc) / A;
    GRADUVW[2] = (yc - ya) / A; // grad v
    GRADUVW[3] = -(xc - xa) / A;
    GRADUVW[4] = (ya - yb) / A; // grad w
    GRADUVW[5] = -(xa - xb) / A;

    sumduda1M = GRADUVW[0] * logM0 + GRADUVW[2] * logM1 + GRADUVW[4] * logM2;
    sumduda2M = GRADUVW[1] * logM0 + GRADUVW[3] * logM1 + GRADUVW[5] * logM2;
  }
}

void feMetric::interpolateMetricP1(const double *x, const double *gammaOrth, Eigen::Matrix2d &M,
                                   Eigen::Matrix2d &sumdudaM)
{
  int elm;
  bool isFound = static_cast<feMesh2DP1 *>(_recovery->_mesh)->locateVertex(x, elm, UVW);
  if(!isFound) {
    feWarning("In interpolateMetricP1 : Point (%f, %f) was not found in the mesh.\n", x[0], x[1]);
    M = Eigen::Matrix2d::Identity();
    return;
  } else {
    // Interpolate the log-metric and take exponential
    Triangle *t = _recovery->_mesh->_elements[elm];

    int t0 = _v2n[t->getVertex(0)];
    int t1 = _v2n[t->getVertex(1)];
    int t2 = _v2n[t->getVertex(2)];

    Eigen::Matrix2d M0 = _metricsOnGmshModel_eigen[t0];
    Eigen::Matrix2d M1 = _metricsOnGmshModel_eigen[t1];
    Eigen::Matrix2d M2 = _metricsOnGmshModel_eigen[t2];

    Eigen::Matrix2d logM0 = M0.log();
    Eigen::Matrix2d logM1 = M1.log();
    Eigen::Matrix2d logM2 = M2.log();

    M = (1. - UVW[0] - UVW[1]) * logM0 + UVW[0] * logM1 + UVW[1] * logM2;
    M = M.exp();

    // Gradient of the barycentric coordinates w.r.t. P11^l
    double xa = t->getVertex(0)->x();
    double ya = t->getVertex(0)->y();
    double xb = t->getVertex(1)->x();
    double yb = t->getVertex(1)->y();
    double xc = t->getVertex(2)->x();
    double yc = t->getVertex(2)->y();
    // Area of the triangle, the 1/2 is be cancelled by the area of the subtriangles
    // associated to the barycentric coordinates. Should we use fabs ?
    double A = (xa - xc) * (yb - yc) - (xb - xc) * (ya - yc);
    // Derivative with respect to alpha in 1D, along the bisector
    double gx = gammaOrth[0];
    double gy = gammaOrth[1];
    double duda = (gx * (yb - yc) - gy * (xb - xc)) / A;
    double dvda = (gx * (yc - ya) - gy * (xc - xa)) / A;
    double dwda = (gx * (ya - yb) - gy * (xa - xb)) / A;

    sumdudaM = duda * logM0 + dvda * logM1 + dwda * logM2;
  }
}

void feMetric::interpolateMetricAndDerivativeOnP2Edge(
  double t, const Eigen::Matrix2d &M11, const Eigen::Matrix2d &M20, const Eigen::Matrix2d &M02,
  const Eigen::Matrix2d &sumduda1M, const Eigen::Matrix2d &sumduda2M, Eigen::Matrix2d &M,
  Eigen::Matrix2d &dMda1, Eigen::Matrix2d &dMda2)
{
  double phi20 = (1. - t) * (1. - 2. * t);
  double phi02 = t * (2. * t - 1.);
  double L = 4. * t * (1. - t); // phi11

  Eigen::Matrix2d T = phi20 * M20.log() + phi02 * M02.log() + L * M11.log(); // = log M(x(t))

  // lT = T;

  // // Scaling and squaring : find m s.t. ||T||/(2^m) <= 1
  // int m = 0;
  // if(lT.norm() > 1){
  //   m = 1;
  //   while(true){
  //     lT /= pow(2,m);
  //     if(lT.norm() <= 1.)
  //       break;
  //     m++;
  //   }
  //   // feInfo("Found m = %d", m);
  //   // feInfo("norm = %f", lT.norm());
  //   // std::cout<<lT<<std::endl;
  // }

  Eigen::Matrix2d Tn = T;
  M = T.exp();

  Eigen::Matrix2d dTda1 = L * sumduda1M;
  Eigen::Matrix2d dTda2 = L * sumduda2M;
  Eigen::Matrix2d dTnda1 = dTda1;
  Eigen::Matrix2d dTnda2 = dTda2;

  // Compute derivative of exponential by recursion
  // Term for k = 0 is the zero matrix (derivative of I)
  // Term for k = 1 is dTdai
  int kfac = 1;
  dMda1 = dTda1;
  dMda2 = dTda2;

  for(int k = 2; k < 10; ++k) {
    dTnda1 = Tn * dTda1 + dTnda1 * T;
    dTnda2 = Tn * dTda2 + dTnda2 * T;

    Tn *= T; // T^n

    kfac *= k;
    dMda1 += 1. / ((double)kfac) * dTnda1;
    dMda2 += 1. / ((double)kfac) * dTnda2;
  }

  // // Compute exponential and derivative of lT instead
  // double lambda = pow(2., -m);
  // Tn = lT;
  // M = lT.exp();
  // M = M.pow(1./lambda);

  // dTda1 = L * sumduda1M;
  // dTda2 = L * sumduda2M;
  // dTnda1 = dTda1;
  // dTnda2 = dTda2;

  // // Compute derivative of exponential by recursion
  // // Term for k = 0 is the zero matrix (derivative of I)
  // // Term for k = 1 is dTdai
  // int kfac = 1;
  // dMda1 = dTda1;
  // dMda2 = dTda2;

  // for(int k = 2; k < 10; ++k){
  //   dTnda1 = Tn * dTda1 + dTnda1 * lT;
  //   dTnda2 = Tn * dTda2 + dTnda2 * lT;

  //   Tn *= lT; // T^n

  //   kfac *= k;
  //   dMda1 += 1./(double) kfac * dTnda1;
  //   dMda2 += 1./(double) kfac * dTnda2;
  // }

  // Amat = M.pow(lambda - 1.);
  // Amat = Amat.inverse();
  // Amat *= 1./lambda;

  // dMda1 = Amat * dMda1;
  // dMda2 = Amat * dMda2;
}

void feMetric::interpolateMetricAndDerivativeOnP2Edge(double t, const Eigen::Matrix2d &M11,
                                                      const Eigen::Matrix2d &M20,
                                                      const Eigen::Matrix2d &M02,
                                                      const Eigen::Matrix2d &sumdudaM,
                                                      Eigen::Matrix2d &M, Eigen::Matrix2d &dMda)
{
  double phi20 = (1. - t) * (1. - 2. * t);
  double phi02 = t * (2. * t - 1.);
  double L = 4. * t * (1. - t); // phi11

  Eigen::Matrix2d T = phi20 * M20.log() + phi02 * M02.log() + L * M11.log(); // = log M(x(t))

  Eigen::Matrix2d Tn = T;
  M = T.exp();

  Eigen::Matrix2d dTda = L * sumdudaM;
  Eigen::Matrix2d dTnda = dTda;

  // Compute derivative of exponential by recursion
  // Term for k = 0 is the zero matrix (derivative of I)
  // Term for k = 1 is dTdai
  int kfac = 1;
  dMda = dTda;

  for(int k = 2; k < 10; ++k) {
    dTnda = Tn * dTda + dTnda * T;
    Tn *= T; // T^n
    kfac *= k;
    dMda += 1. / ((double)kfac) * dTnda;
  }
}

void interpolateMetricP1WithDerivativesWrapper(void *metric, const double *x, Eigen::Matrix2d &M,
                                               Eigen::Matrix2d &dMdx, Eigen::Matrix2d &dMdy)
{
  feMetric *m = static_cast<feMetric *>(metric);
  m->interpolateMetricP1WithDerivatives(x, M, dMdx, dMdy);
}

void interpolateMetricP1Wrapper(void *metric, const double *x, Eigen::Matrix2d &M,
                                Eigen::Matrix2d &sumduda1M, Eigen::Matrix2d &sumduda2M)
{
  feMetric *m = static_cast<feMetric *>(metric);
  m->interpolateMetricP1(x, M, sumduda1M, sumduda2M);
}

void interpolateMetricP1Wrapper1D(void *metric, const double *x, const double *gammaOrth,
                                  Eigen::Matrix2d &M, Eigen::Matrix2d &sumdudaM)
{
  feMetric *m = static_cast<feMetric *>(metric);
  m->interpolateMetricP1(x, gammaOrth, M, sumdudaM);
}

void interpolateMetricAndDerivativeOnP2EdgeWrapper(
  void *metric, const double t, const Eigen::Matrix2d &M11, const Eigen::Matrix2d &M20,
  const Eigen::Matrix2d &M02, const Eigen::Matrix2d &sumduda1M, const Eigen::Matrix2d &sumduda2M,
  Eigen::Matrix2d &M, Eigen::Matrix2d &dMda1, Eigen::Matrix2d &dMda2)
{
  feMetric *m = static_cast<feMetric *>(metric);
  m->interpolateMetricAndDerivativeOnP2Edge(t, M11, M20, M02, sumduda1M, sumduda2M, M, dMda1,
                                            dMda2);
}

void interpolateMetricAndDerivativeOnP2EdgeWrapper1D(void *metric, const double t,
                                                     const Eigen::Matrix2d &M11,
                                                     const Eigen::Matrix2d &M20,
                                                     const Eigen::Matrix2d &M02,
                                                     const Eigen::Matrix2d &sumdudaM,
                                                     Eigen::Matrix2d &M, Eigen::Matrix2d &dMda)
{
  feMetric *m = static_cast<feMetric *>(metric);
  m->interpolateMetricAndDerivativeOnP2Edge(t, M11, M20, M02, sumdudaM, M, dMda);
}