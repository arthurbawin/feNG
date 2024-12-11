
#include "geodesic.h"
#include "feMessage.h"

#if defined(HAVE_CERES)
#include "/usr/local/include/ceres/ceres.h"
using ceres::CostFunction;
using ceres::Problem;
using ceres::SizedCostFunction;
using ceres::Solve;
using ceres::Solver;
#endif

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
	const double ds, double *length = nullptr, const double *x1 = nullptr,
  const double *distToPoint = nullptr)
{
  double X[4] = {x0[0], x0[1], u0[0], u0[1]};
  double k1[4], k2[4], k3[4], k4[4], Xcpy[4];

  double M0[2][2], M[2][2], dMdx[2][2], dMdy[2][2];

  double L = 0.;
  double fPrev = 0., f = 0.;

  std::vector<SPoint2> geodesic;

  double mindist = 1e10;
  double toldist = 5*ds;
  if(distToPoint) toldist = 1.2 * (*distToPoint);
  bool reached = false;

  // Integrate the geodesic equation until L = targetLength
  while(L <= targetLength) {

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

    if(isnan(X[0]) || isnan(X[1])) {
      break;
    }

    geodesic.push_back(SPoint2(X[0],X[1]));

    // Get initial metric to compute edge length
    assignArrayFrom(M,M0);
    
    // Current velocity vector is (X[2],X[3])
    // Approach the length with a trapeze rule for integration
    f = sqrt(X[2]*X[2]*M[0][0] + 2.*X[2]*X[3]*M[0][1] + X[3]*X[3]*M[1][1]);
    L += (f + fPrev)/2. * ds;
    fPrev = f;

    // Stop if point is close enough (in Euclidean distance) to given endpoint x1
    if(x1) {
      double dist = sqrt( (x1[0]-X[0])*(x1[0]-X[0]) + (x1[1]-X[1])*(x1[1]-X[1]) );
      mindist = fmin(mindist, dist);
      // tolerance for bisection is 2*ds, so check if dist is a bit greater than that
      if(dist < toldist) {
        reached = true;
        break;
      }
    }
  }
  if(length) *length = L;

  if(x1 && !reached) {
    feInfo("Geodesic did not reach endpoint: min dist to endpoint = %+-1.4e - threshold = %+-1.4e", mindist, toldist);
  }
  return geodesic;
}

inline double sign(double x){return x > 0 ? 1. : (x < 0 ? -1. : 0.); };

// Compute the signed distance from a geodesic to a point
double distGeodesicToPoint(const double x[2], const std::vector<SPoint2> &geodesic)
{
	double mindist = 1e10, dist, dotprod = 1e10, orient = 1.;
	size_t numPoints = geodesic.size();
  if(numPoints > 0) {
  	for(size_t i = 0; i < numPoints; ++i) {
  		SPoint2 p = geodesic[i];
  		dist = sqrt((x[0] - p[0])*(x[0] - p[0]) + (x[1] - p[1])*(x[1] - p[1]));
  		if(dist < mindist) {
  			mindist = dist;
        if(i == numPoints-1) {
          double AB[2] = {p[0] - geodesic[i-1][0], p[1] - geodesic[i-1][1]};
          double Ax[2] = {            x[0] - p[0],             x[1] - p[1]};
          double crossprod = AB[0]*Ax[1] - AB[1]*Ax[0];
          dotprod = AB[0]*Ax[0] + AB[1]*Ax[1];
          orient = sign(crossprod);
        } else {
    			double AB[2] = {geodesic[i+1][0] - p[0], geodesic[i+1][1] - p[1]};
    			double Ax[2] = {            x[0] - p[0],             x[1] - p[1]};
          double crossprod = AB[0]*Ax[1] - AB[1]*Ax[0];
          dotprod = AB[0]*Ax[0] + AB[1]*Ax[1];
          orient = sign(crossprod);
        }
  		}
  	}
    // feInfo("%f - %f - %+-1.10e", orient, mindist, dotprod);
    if(fabs(orient) < 1e-10) {
      // Points are colinear: return dotProduct?
      return sign(dotprod) * mindist;
    } else {
  	 return orient * mindist;
    }
  }
  return 1e10;
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
	return (geodesic.size() > 0) ? distGeodesicToPoint(x1, geodesic) : 1e10;
}

// Determine the angle along which to shoot to get from x0 to x1 along a geodesic
double bisectionAngle(const double x0[2], const double x1[2], const MetricFunctionPtr getMetric, 
	const double maxLength, const int maxIter, const double tol, const double ds)
{
	double d[2] = {x1[0]-x0[0], x1[1]-x0[1]};
	double theta0 = atan2(d[1], d[0]);
	double targetLength = maxLength;
	double theta_a = theta0 - M_PI/4.; double fa = bisectionAngle_callback(theta_a, x0, x1, getMetric, targetLength, ds);
	double theta_b = theta0 + M_PI/4.; double fb = bisectionAngle_callback(theta_b, x0, x1, getMetric, targetLength, ds);
  if(fa*fb > 0) {
    feErrorMsg(FE_STATUS_ERROR, "Bisection bracket for angle (geodesic) does not work: f(%+-1.3e) = %+-1.3e and f(%+-1.3e) = %+-1.3e",
      theta_a, fa, theta_b, fb);
  }
  int iter = 0;
  double theta = 0.;
  while(iter++ <= maxIter) {
    theta = (theta_a + theta_b)/2.;
    double fc = bisectionAngle_callback(theta, x0, x1, getMetric, targetLength, ds);
    // feInfo("Bisection for theta: theta = %+-1.3e - f(theta) = %+-1.3e", theta, fc);
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
	const double maxLength, const int maxIter, const double tol, const double ds, double &distToPoint)
{
	double d[2] = {x1[0]-x0[0], x1[1]-x0[1]};
	double theta0 = atan2(d[1], d[0]);
	double targetLength = maxLength;
	double theta_a = theta0 - M_PI/4.; double fa = bisectionAngle_callback(theta_a, x0, x1, getMetric, targetLength, ds);
	double theta_b = theta0 + M_PI/4.; double fb = bisectionAngle_callback(theta_b, x0, x1, getMetric, targetLength, ds);
  if(fa*fb > 0) {
    // Try again with larger angle
    theta_a = theta0 - M_PI/2.; fa = bisectionAngle_callback(theta_a, x0, x1, getMetric, targetLength, ds);
    theta_b = theta0 + M_PI/2.; fb = bisectionAngle_callback(theta_b, x0, x1, getMetric, targetLength, ds);
    if(fa*fb > 0) {
      feErrorMsg(FE_STATUS_ERROR, "Bisection bracket for angle (geodesic) does not work: f(%+-1.3e) = %+-1.3e and f(%+-1.3e) = %+-1.3e",
        theta_a, fa, theta_b, fb);
    }
  }
  int iter = 0;
  double theta = 0.;
  while(iter++ <= maxIter) {
    theta = (theta_a + theta_b)/2.;
    double fc = bisectionAngle_callback(theta, x0, x1, getMetric, targetLength, ds);
    // feInfo("Bisection for theta: theta = %+-1.3e - f(theta) = %+-1.3e", theta, fc);
    if(fabs(fc) < tol || (theta_b-theta_a)/2. < tol) {
      distToPoint = fabs(fc);
      // feInfo("distToPoint = %f for theta = %f - ds = %f - x0 = %f,%f - x1 = %f,%f", distToPoint, theta, ds, x0[0], x0[1], x1[0], x1[1]);
      // feInfo("Recomputing: %f", bisectionAngle_callback(theta, x0, x1, getMetric, targetLength, ds));
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
	const double maxLength, const int maxIter, const double tol, const double ds, double &length, double &error)
{
	if(fabs(ds) > fabs(tol)) {
		feErrorMsg(FE_STATUS_ERROR, "Geodesic: step size ds should be smaller than the tolerance for bisections");
		exit(-1);
	}

	double theta = bisectionAngle(x0, x1, getMetric, maxLength, maxIter, tol, ds);
	double u0[2] = {cos(theta), sin(theta)};
	double targetLength = maxLength; // Un peu plus que 1, il faudrait faire une deuxième bissection
	double L;
	std::vector<SPoint2> geodesic = integrateRK4(x0, u0, getMetric, targetLength, ds, &L);
	length = L;

	// Compute the error: min distance between geodesic and x1
	error = distGeodesicToPoint(x1, geodesic);

	return geodesic;
}

std::vector<SPoint2> geodesicBetweenTwoPoints_arrayPtr(const double x0[2], const double x1[2], const MetricFunctionPtrArray getMetric,
	const double maxLength, const int maxIter, const double tol, const double ds, double &length, double &error)
{
	if(fabs(ds) > fabs(tol)) {
		feErrorMsg(FE_STATUS_ERROR, "Geodesic: step size ds should be smaller than the tolerance for bisections");
		exit(-1);
	}

  double distToEndpoint = 1e10;
	double theta = bisectionAngle(x0, x1, getMetric, maxLength, maxIter, tol, ds, distToEndpoint);
	double u0[2] = {cos(theta), sin(theta)};
	double targetLength = 2 * maxLength; // Un peu plus que 1, il faudrait faire une deuxième bissection
	double L;
  // feInfo("Computing for theta = %f - ds = %f - x0 = %f,%f - x1 = %f,%f", theta, ds, x0[0], x0[1], x1[0], x1[1]);
  std::vector<SPoint2> geodesic = integrateRK4(x0, u0, getMetric, targetLength, ds, &L, x1, &distToEndpoint);
	length = L;

	// Compute the error: min distance between geodesic and x1
	error = geodesic.size() > 0 ? distGeodesicToPoint(x1, geodesic) : 1e10;

	return geodesic;
}

// /* 10 point rule points */
// static int nQuad = 10;
// static double xQuad[10] = {-9.739065285171716e-01, -8.650633666889845e-01,
//                            -6.794095682990244e-01, -4.333953941292472e-01,
//                            -1.488743389816312e-01, 1.488743389816312e-01,
//                            4.333953941292472e-01,  6.794095682990244e-01,
//                            8.650633666889845e-01,  9.739065285171716e-01};

// /* 10 point rule weights */
// static double wQuad[10] = {6.667134430868774e-02, 1.494513491505805e-01,
//                            2.190863625159822e-01, 2.692667193099962e-01,
//                            2.955242247147529e-01, 2.955242247147529e-01,
//                            2.692667193099962e-01, 2.190863625159822e-01,
//                            1.494513491505805e-01, 6.667134430868774e-02};

/* 20 point rule points */
static int nQuad = 20;
static double xQuad[20] = {
  -9.931285991850949e-01, -9.639719272779138e-01, -9.122344282513259e-01,
  -8.391169718222188e-01, -7.463319064601508e-01, -6.360536807265150e-01,
  -5.108670019508271e-01, -3.737060887154196e-01, -2.277858511416451e-01,
  -7.652652113349733e-02, 7.652652113349733e-02,  2.277858511416451e-01,
  3.737060887154196e-01,  5.108670019508271e-01,  6.360536807265150e-01,
  7.463319064601508e-01,  8.391169718222188e-01,  9.122344282513259e-01,
  9.639719272779138e-01,  9.931285991850949e-01};

/* 20 point rule weights */
static double wQuad[20] = {
  1.761400713915212e-02, 4.060142980038694e-02, 6.267204833410906e-02,
  8.327674157670475e-02, 1.019301198172404e-01, 1.181945319615184e-01,
  1.316886384491766e-01, 1.420961093183821e-01, 1.491729864726037e-01,
  1.527533871307259e-01, 1.527533871307259e-01, 1.491729864726037e-01,
  1.420961093183821e-01, 1.316886384491766e-01, 1.181945319615184e-01,
  1.019301198172404e-01, 8.327674157670475e-02, 6.267204833410906e-02,
  4.060142980038694e-02, 1.761400713915212e-02};

static double METRIC[2][2], DMDX[2][2], DMDY[2][2];

// Assume symmetric matrix
static double dotProduct(const double m[2][2], const double x1[2], const double x2[2])
{
  return x2[0] * (m[0][0] * x1[0] + m[0][1] * x1[1])
       + x2[1] * (m[0][1] * x1[0] + m[1][1] * x1[1]);
}

double computeEdgeLengthInMetric(const double x0[2], const double x1[2], const MetricFunctionPtrArray getMetric)
{
  double L = 0.;
  double x[2], diff[2] = {x1[0]-x0[0], x1[1]-x0[1]};
  for(int i = 0; i < nQuad; ++i)
  {
    // Map from [-1,1] to [0,1]
    double t = (xQuad[i] + 1.)/2.;
    double w = wQuad[i]/2.;
    x[0] = x0[0] + diff[0] * t;
    x[1] = x0[1] + diff[1] * t;
    getMetric(METRIC, DMDX, DMDY, x);
    L += w * sqrt(dotProduct(METRIC,diff,diff));
  }
  return L;
}

double computeEdgeLengthSquaredInMetric(const double x0[2], const double x1[2], const MetricFunctionPtrArray getMetric)
{
  double L = 0.;
  double x[2], diff[2] = {x1[0]-x0[0], x1[1]-x0[1]};
  for(int i = 0; i < nQuad; ++i)
  {
    // Map from [-1,1] to [0,1]
    double t = (xQuad[i] + 1.)/2.;
    double w = wQuad[i]/2.;
    x[0] = x0[0] + diff[0] * t;
    x[1] = x0[1] + diff[1] * t;
    getMetric(METRIC, DMDX, DMDY, x);
    L += w * dotProduct(METRIC,diff,diff);
  }
  return L;
}

double length_costFunction(const std::vector<Vertex*> &vertices, const MetricFunctionPtrArray getMetric)
{
	double cost = 0.;
	for(size_t i = 0; i < vertices.size()-1; ++i) {
		double x0[2] = {vertices[i]->x(), vertices[i]->y()};
		double x1[2] = {vertices[i+1]->x(), vertices[i+1]->y()};
		cost += computeEdgeLengthInMetric(x0, x1, getMetric);
	}
	return cost;
}

double energy_costFunction(const std::vector<Vertex*> &vertices, const MetricFunctionPtrArray getMetric)
{
	double cost = 0.;
	for(size_t i = 0; i < vertices.size()-1; ++i) {
		double x0[2] = {vertices[i]->x(), vertices[i]->y()};
		double x1[2] = {vertices[i+1]->x(), vertices[i+1]->y()};
		cost += 0.5 * computeEdgeLengthSquaredInMetric(x0, x1, getMetric);
	}
	return cost;
}

#if defined(HAVE_CERES)
class costFunctorGeodesic final : public ceres::FirstOrderFunction {
  public:
    size_t _numVertices;
    int _numParameters;
    std::vector<Vertex*> _vertices;
    MetricFunctionPtrArray _getMetric;
    mutable std::vector<double> INITIALPOS;

  public:
    costFunctorGeodesic(std::vector<Vertex*> &vertices, MetricFunctionPtrArray getMetric, int numParameters)
      : _vertices(vertices), _getMetric(getMetric), _numParameters(numParameters)
      {
      	if(vertices.size() <= 2) {
      		feInfo("Cannot create a costFunctorGeodesic for less than 3 points");
      		exit(-1);
      	}
      	_numVertices = vertices.size();
        INITIALPOS.resize(numParameters, 0.);
      };

    // The objective (cost) function and its gradient
    bool Evaluate(const double* parameters,
                  double* cost,
                  double* gradient) const override
    {
      // Save initial positions
      int cnt = 0;
      for(int i = 1; i < _numVertices-1; ++i) {
        Vertex *v = _vertices[i];
        INITIALPOS[cnt] = v->x(); cnt++;
        INITIALPOS[cnt] = v->y(); cnt++;
      }

      // Update the positions
      cnt = 0;
      for(int i = 1; i < _numVertices-1; ++i) {
        Vertex *v = _vertices[i];
        (*v)(0) = parameters[cnt]; cnt++;
        (*v)(1) = parameters[cnt]; cnt++;
      }

      // COMPUTE OBJECTIVE FUNCTION
      cost[0] = length_costFunction(_vertices, _getMetric);
      // cost[0] = energy_costFunction(_vertices, _getMetric);

      // COMPUTE GRADIENT
      if (gradient) {

        // for(int i = 0; i < _numParameters; ++i)
        //   gradient[i] = 0.;

        // cnt = 0;
        // double vertexGradient[2] = {0., 0.};
        // for(int i = 0; i < _numVerticesToModify; ++i) {
        //   computeCostGradient(i, _vertices, _vertNeighbours, _elements, _vert2tri, vertexGradient);
        //   gradient[cnt] = vertexGradient[0]; cnt++;
        //   gradient[cnt] = vertexGradient[1]; cnt++;
        // }
        feInfo("Pas de gradient implémenté");
        exit(-1);
      }

      // Restore initial positions
      cnt = 0;
      for(int i = 1; i < _numVertices-1; ++i) {
        Vertex *v = _vertices[i];
        (*v)(0) = INITIALPOS[cnt]; cnt++;
        (*v)(1) = INITIALPOS[cnt]; cnt++;
      }

      return true;
    }

    // When using numeric differentiation: only compute cost
    bool operator()(const double* parameters, double* cost) const 
    {
      return this->Evaluate(parameters, cost, nullptr);
    }

    int NumParameters() const override { return _numParameters; }
};
#endif

void geodesicBetweenTwoPoints_arrayPtr_minimizeLength(std::vector<Vertex> &vertices, 
	const MetricFunctionPtrArray getMetric, double &length)
{
#if defined(HAVE_CERES)
	std::vector<Vertex*> verticesPtr(vertices.size());
	for(size_t i = 0; i < vertices.size(); ++i) {
		verticesPtr[i] = &vertices[i];
	}	

  int numParameters = 2 * (vertices.size() - 2);
  costFunctorGeodesic *costFun = new costFunctorGeodesic(verticesPtr, getMetric, numParameters);
  
  // Analytic diff
  // ceres::GradientProblem problem(costFun);

  // Numeric diff (only to test)
  auto *foo = new ceres::NumericDiffFirstOrderFunction<costFunctorGeodesic, ceres::CENTRAL>(costFun, numParameters);
  ceres::GradientProblem problem(foo);

  ceres::GradientProblemSolver::Options options;
  // options.function_tolerance = 1e-12;
  // options.gradient_tolerance = 1e-10;
  // options.parameter_tolerance = 1e-12;
  // // options.min_line_search_step_size = 1e-11;
  options.max_num_iterations = 50000;
  // options.max_num_line_search_step_size_iterations = 1000;

  // options.line_search_direction_type = ceres::LBFGS;
  // options.line_search_type = ceres::WOLFE;
  // options.line_search_interpolation_type = ceres::CUBIC;
  // options.max_lbfgs_rank = 5000;

  options.minimizer_progress_to_stdout = true;

  ceres::GradientProblemSolver::Summary summary;

  // Initial positions
  std::vector<double> optimalPos(numParameters, 0.);
  int cnt = 0;
  for(size_t i = 1; i < vertices.size()-1; ++i) {
    Vertex *v = verticesPtr[i];
    optimalPos[cnt] = v->x(); cnt++;
    optimalPos[cnt] = v->y(); cnt++;
  }

  // feInfo("Computing geodesic...");
  ceres::Solve(options, problem, optimalPos.data(), &summary);

  std::cout << summary.FullReport() << "\n";

  if(summary.termination_type == ceres::TerminationType::CONVERGENCE || 
     summary.termination_type == ceres::TerminationType::NO_CONVERGENCE)
  {
    // Converged or ran out of iterations, but did not fail
    // Continue with the last computed candidate
    // return summary.final_cost;
  } else {
    // Failed
    feErrorMsg(FE_STATUS_ERROR, "Minimization did not converge");
  }

  // Apply new position
  cnt = 0;
  for(int i = 1; i < vertices.size()-1; ++i) {
    Vertex *v = verticesPtr[i];
    (*v)(0) = optimalPos[cnt]; cnt++;
    (*v)(1) = optimalPos[cnt]; cnt++;
  }
#else
  UNUSED(vertices, getMetric, length);
#endif
}