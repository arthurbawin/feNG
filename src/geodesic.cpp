
#include "geodesic.h"
#include "feMessage.h"

void christoffel(const MetricTensor &M, const MetricTensor &dMdx, const MetricTensor &dMdy, double res[2][2][2])
{
  double g[2][2][2];
  // Compute g = d(g_ab)/dx^c
  for(int i = 0; i < 2; ++i){
    for(int j = 0; j < 2; ++j){
      g[i][j][0] = dMdx(i,j);
      g[i][j][1] = dMdy(i,j);
    }
  }

  // Compute M^-1
  double minv[2][2];
  double det = M(0,0)*M(1,1) - M(0,1)*M(0,1);
  minv[0][0] =  M(1,1)/det;
  minv[0][1] = -M(0,1)/det;
  minv[1][0] = -M(1,0)/det;
  minv[1][1] =  M(0,0)/det;

  // res = {{{0., 0.},{0., 0.}},{{0., 0.},{0., 0.}}};
  for(int i = 0; i < 2; ++i) {
    for(int k = 0; k < 2; ++k) {
      for(int l = 0; l < 2; ++l) { res[i][k][l] = 0.; }
  	}
	}

  for(int i = 0; i < 2; ++i) {
    for(int k = 0; k < 2; ++k) {
      for(int l = 0; l < 2; ++l) {
        for(int m = 0; m < 2; ++m) {
          res[i][k][l] += 0.5 * (minv[i][m]*(g[m][k][l] + g[m][l][k] - g[k][l][m]));
        }
      }
    }
  }
}

void christoffel(const double M[2][2], const double dMdx[2][2], const double dMdy[2][2], double res[2][2][2])
{
  double g[2][2][2];
  // Compute g = d(g_ab)/dx^c
  for(int i = 0; i < 2; ++i){
    for(int j = 0; j < 2; ++j){
      g[i][j][0] = dMdx[i][j];
      g[i][j][1] = dMdy[i][j];
    }
  }

  // Compute M^-1
  double minv[2][2];
  double det = M[0][0]*M[1][1] - M[0][1]*M[0][1];
  minv[0][0] =  M[1][1]/det;
  minv[0][1] = -M[0][1]/det;
  minv[1][0] = -M[1][0]/det;
  minv[1][1] =  M[0][0]/det;

  // res = {{{0., 0.},{0., 0.}},{{0., 0.},{0., 0.}}};
  for(int i = 0; i < 2; ++i) {
    for(int k = 0; k < 2; ++k) {
      for(int l = 0; l < 2; ++l) { res[i][k][l] = 0.; }
  	}
	}

  for(int i = 0; i < 2; ++i) {
    for(int k = 0; k < 2; ++k) {
      for(int l = 0; l < 2; ++l) {
        for(int m = 0; m < 2; ++m) {
          res[i][k][l] += 0.5 * (minv[i][m]*(g[m][k][l] + g[m][l][k] - g[k][l][m]));
        }
      }
    }
  }
}

void getRHS(const double X[4], const MetricTensor &M, const MetricTensor &dMdx, const MetricTensor &dMdy, double F[4])
{
  double res[2][2][2] = {{{0., 0.},{0., 0.}},{{0., 0.},{0., 0.}}};
  christoffel(M, dMdx, dMdy, res);

  F[0] = F[1] = F[2] = F[3] = 0.;
  F[0] = X[2];
  F[1] = X[3];

  for(int b = 0; b < 2; ++b){
    for(int c = 0; c < 2; ++c){
      F[2] -= res[0][b][c] * X[b+2] * X[c+2];
      F[3] -= res[1][b][c] * X[b+2] * X[c+2];
    }
  }
}

void getRHS(const double X[4], const double M[2][2], const double dMdx[2][2], const double dMdy[2][2], double F[4])
{
  double res[2][2][2] = {{{0., 0.},{0., 0.}},{{0., 0.},{0., 0.}}};
  christoffel(M, dMdx, dMdy, res);

  F[0] = F[1] = F[2] = F[3] = 0.;
  F[0] = X[2];
  F[1] = X[3];

  for(int b = 0; b < 2; ++b){
    for(int c = 0; c < 2; ++c){
      F[2] -= res[0][b][c] * X[b+2] * X[c+2];
      F[3] -= res[1][b][c] * X[b+2] * X[c+2];
    }
  }
}

// Integrate geodesic equations until prescribed length
std::vector<SPoint2> integrateRK4(const double x0[2], const double u0[2], 
	const MetricFunctionPtr getMetric, const double targetLength,
	const double ds, double *length = nullptr)
{
  double X[4] = {x0[0], x0[1], u0[0], u0[1]};
  double k1[4], k2[4], k3[4], k4[4], Xcpy[4];

  MetricTensor M0, M, dMdx, dMdy;

  double L = 0.;
  double fPrev = 0., f = 0.;

  std::vector<SPoint2> geodesic;

  // Integrate the geodesic equation until L = targetLength
  while(L <= targetLength){

    getMetric(M, dMdx, dMdy, X);
    M0.assignMatrixFrom(M);
    getRHS(X, M, dMdx, dMdy, k1);

    for(size_t i = 0; i < 4; ++i)
      Xcpy[i] = X[i] + ds/2. * k1[i];

    getMetric(M, dMdx, dMdy, Xcpy);
    getRHS(Xcpy, M, dMdx, dMdy, k2);

    for(size_t i = 0; i < 4; ++i)
      Xcpy[i] = X[i] + ds/2. * k2[i];

    getMetric(M, dMdx, dMdy, Xcpy);
    getRHS(Xcpy, M, dMdx, dMdy, k3);

    for(size_t i = 0; i < 4; ++i)
      Xcpy[i] = X[i] + ds * k3[i];

    getMetric(M, dMdx, dMdy, Xcpy);
    getRHS(Xcpy, M, dMdx, dMdy, k4);

    for(size_t i = 0; i < 4; ++i)
      X[i] += 1./6. * ds * (k1[i] + 2.*k2[i] + 2.*k3[i] + k4[i]);

    geodesic.push_back(SPoint2(X[0],X[1]));

    // Get initial metric to compute edge length
    M.assignMatrixFrom(M0);
    
    // Current velocity vector is (X[2],X[3])
    // Approach the length with a trapeze rule for integration
    f = sqrt(X[2]*X[2]*M(0,0) + 2.*X[2]*X[3]*M(0,1) + X[3]*X[3]*M(1,1));
    L += (f + fPrev)/2. * ds;
    fPrev = f;
  }
  if(length) *length = L;
  return geodesic;
}

void assignArrayFrom(double toAssign[2][2], const double other[2][2])
{
	toAssign[0][0] = other[0][0];
	toAssign[0][1] = other[0][1];
	toAssign[1][0] = other[1][0];
	toAssign[1][1] = other[1][1];
}

// Integrate geodesic equations until prescribed length
std::vector<SPoint2> integrateRK4(const double x0[2], const double u0[2], 
	const MetricFunctionPtrArray getMetric, const double targetLength,
	const double ds, double *length = nullptr)
{
  double X[4] = {x0[0], x0[1], u0[0], u0[1]};
  double k1[4], k2[4], k3[4], k4[4], Xcpy[4];

  double M0[2][2], M[2][2], dMdx[2][2], dMdy[2][2];

  double L = 0.;
  double fPrev = 0., f = 0.;

  std::vector<SPoint2> geodesic;

  // Integrate the geodesic equation until L = targetLength
  while(L <= targetLength){

    getMetric(M, dMdx, dMdy, X);
    assignArrayFrom(M0,M);
    getRHS(X, M, dMdx, dMdy, k1);

    for(size_t i = 0; i < 4; ++i)
      Xcpy[i] = X[i] + ds/2. * k1[i];

    getMetric(M, dMdx, dMdy, Xcpy);
    getRHS(Xcpy, M, dMdx, dMdy, k2);

    for(size_t i = 0; i < 4; ++i)
      Xcpy[i] = X[i] + ds/2. * k2[i];

    getMetric(M, dMdx, dMdy, Xcpy);
    getRHS(Xcpy, M, dMdx, dMdy, k3);

    for(size_t i = 0; i < 4; ++i)
      Xcpy[i] = X[i] + ds * k3[i];

    getMetric(M, dMdx, dMdy, Xcpy);
    getRHS(Xcpy, M, dMdx, dMdy, k4);

    for(size_t i = 0; i < 4; ++i)
      X[i] += 1./6. * ds * (k1[i] + 2.*k2[i] + 2.*k3[i] + k4[i]);

    geodesic.push_back(SPoint2(X[0],X[1]));

    // Get initial metric to compute edge length
    assignArrayFrom(M,M0);
    
    // Current velocity vector is (X[2],X[3])
    // Approach the length with a trapeze rule for integration
    f = sqrt(X[2]*X[2]*M[0][0] + 2.*X[2]*X[3]*M[0][1] + X[3]*X[3]*M[1][1]);
    L += (f + fPrev)/2. * ds;
    fPrev = f;
  }
  if(length) *length = L;
  return geodesic;
}

inline double sign(double x){return x > 0 ? 1. : (x < 0 ? -1. : 0.); };

// Compute the signed distance from a geodesic to a point
double distGeodesicToPoint(const double x[2], const std::vector<SPoint2> &geodesic)
{
	double mindist = 1e10, dist, orient;
	size_t numPoints = geodesic.size();
	for(size_t i = 0; i < numPoints-1; ++i) {
		SPoint2 p = geodesic[i];
		dist = sqrt((x[0] - p[0])*(x[0] - p[0]) + (x[1] - p[1])*(x[1] - p[1]));
		if(dist < mindist) {
			mindist = dist;
			double AB[2] = {geodesic[i+1][0] - p[0], geodesic[i+1][1] - p[1]};
			double Ax[2] = {            x[0] - p[0],             x[1] - p[1]};
			double crossprod = AB[0]*Ax[1] - AB[1]*Ax[0];
			orient = sign(crossprod);
		}
	}
	return orient * mindist;
}

// Integrates geodesic until target length, then return a positive value 
// if the closest of the geodesic is "left" of the target point,
// negative value if "right" and 0 if it passes through the point
double bisectionAngle_callback(const double theta, const double x0[2], const double x1[2], const MetricFunctionPtr getMetric,
	const double targetLength, const double ds)
{
	double u0[2] = {cos(theta), sin(theta)};
	std::vector<SPoint2> geodesic = integrateRK4(x0, u0, getMetric, targetLength, ds);
	return distGeodesicToPoint(x1, geodesic);
}

double bisectionAngle_callback(const double theta, const double x0[2], const double x1[2], const MetricFunctionPtrArray getMetric,
	const double targetLength, const double ds)
{
	double u0[2] = {cos(theta), sin(theta)};
	std::vector<SPoint2> geodesic = integrateRK4(x0, u0, getMetric, targetLength, ds);
	return distGeodesicToPoint(x1, geodesic);
}

// Determine the angle along which to shoot to get from x0 to x1 along a geodesic
double bisectionAngle(const double x0[2], const double x1[2], const MetricFunctionPtr getMetric, 
	const int maxIter, const double tol, const double ds)
{
	double d[2] = {x1[0]-x0[0], x1[1]-x0[1]};
	double theta0 = atan2(d[1], d[0]);
	double targetLength = 1.1;
	double theta_a = theta0 - M_PI/4.; double fa = bisectionAngle_callback(theta_a, x0, x1, getMetric, targetLength, ds);
	double theta_b = theta0 + M_PI/4.; double fb = bisectionAngle_callback(theta_b, x0, x1, getMetric, targetLength, ds);
  if(fa*fb > 0) {
    feErrorMsg(FE_STATUS_ERROR, "Bisection bracket for angle (geodesic) does not work: f(%+-1.3e) = %+-1.3e and f(%+-1.3e) = %+-1.3e",
      theta_a, fa, theta_b, fb);
  }
  int iter = 0;
  double theta;
  while(iter++ <= maxIter) {
    theta = (theta_a + theta_b)/2.;
    double fc = bisectionAngle_callback(theta, x0, x1, getMetric, targetLength, ds);
    feInfo("Bisection for theta: theta = %+-1.3e - f(theta) = %+-1.3e", theta, fc);
    if(fabs(fc) < tol || (theta_b-theta_a)/2. < tol) {
      return theta;
    }
    if(fa*fc < 0) {
      theta_b = theta;
      fb = fc;
    } else {
      theta_a = theta;
      fa = fc;
    }
  }
  feErrorMsg(FE_STATUS_ERROR, "Bisection for theta did not converge :/");
  return theta;
}

double bisectionAngle(const double x0[2], const double x1[2], const MetricFunctionPtrArray getMetric, 
	const int maxIter, const double tol, const double ds)
{
	double d[2] = {x1[0]-x0[0], x1[1]-x0[1]};
	double theta0 = atan2(d[1], d[0]);
	double targetLength = 1.1;
	double theta_a = theta0 - M_PI/4.; double fa = bisectionAngle_callback(theta_a, x0, x1, getMetric, targetLength, ds);
	double theta_b = theta0 + M_PI/4.; double fb = bisectionAngle_callback(theta_b, x0, x1, getMetric, targetLength, ds);
  if(fa*fb > 0) {
    feErrorMsg(FE_STATUS_ERROR, "Bisection bracket for angle (geodesic) does not work: f(%+-1.3e) = %+-1.3e and f(%+-1.3e) = %+-1.3e",
      theta_a, fa, theta_b, fb);
  }
  int iter = 0;
  double theta;
  while(iter++ <= maxIter) {
    theta = (theta_a + theta_b)/2.;
    double fc = bisectionAngle_callback(theta, x0, x1, getMetric, targetLength, ds);
    feInfo("Bisection for theta: theta = %+-1.3e - f(theta) = %+-1.3e", theta, fc);
    if(fabs(fc) < tol || (theta_b-theta_a)/2. < tol) {
      return theta;
    }
    if(fa*fc < 0) {
      theta_b = theta;
      fb = fc;
    } else {
      theta_a = theta;
      fa = fc;
    }
  }
  feErrorMsg(FE_STATUS_ERROR, "Bisection for theta did not converge :/");
  return theta;
}

std::vector<SPoint2> geodesicBetweenTwoPoints(const double x0[2], const double x1[2], const MetricFunctionPtr getMetric,
	const int maxIter, const double tol, const double ds, double &length, double &error)
{
	if(fabs(ds) > fabs(tol)) {
		feErrorMsg(FE_STATUS_ERROR, "Geodesic: step size ds should be smaller than the tolerance for bisections");
		exit(-1);
	}

	double theta = bisectionAngle(x0, x1, getMetric, maxIter, tol, ds);
	double u0[2] = {cos(theta), sin(theta)};
	double targetLength = 1.2; // Un peu plus que 1, il faudrait faire une deuxième bissection
	double L;
	std::vector<SPoint2> geodesic = integrateRK4(x0, u0, getMetric, targetLength, ds, &L);
	length = L;

	// Compute the error: min distance between geodesic and x1
	error = distGeodesicToPoint(x1, geodesic);

	return geodesic;
}

std::vector<SPoint2> geodesicBetweenTwoPoints_arrayPtr(const double x0[2], const double x1[2], const MetricFunctionPtrArray getMetric,
	const int maxIter, const double tol, const double ds, double &length, double &error)
{
	if(fabs(ds) > fabs(tol)) {
		feErrorMsg(FE_STATUS_ERROR, "Geodesic: step size ds should be smaller than the tolerance for bisections");
		exit(-1);
	}

	double theta = bisectionAngle(x0, x1, getMetric, maxIter, tol, ds);
	double u0[2] = {cos(theta), sin(theta)};
	double targetLength = 1.2; // Un peu plus que 1, il faudrait faire une deuxième bissection
	double L;
	std::vector<SPoint2> geodesic = integrateRK4(x0, u0, getMetric, targetLength, ds, &L);
	length = L;

	// Compute the error: min distance between geodesic and x1
	error = distGeodesicToPoint(x1, geodesic);

	return geodesic;
}