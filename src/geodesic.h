#ifndef _FEGEODESIC_
#define _FEGEODESIC_

#include "feNG.h"
#include "feVertex.h"
#include "feMatrixInterface.h"
#include "../contrib/gmsh/SPoint2.h"

typedef void (*MetricFunctionPtr)(MetricTensor &m, MetricTensor &dmdx, MetricTensor &dmdy, double pos[2]);
typedef void (*MetricFunctionPtrArray)(double m[2][2], double dmdx[2][2], double dmdy[2][2], double pos[2]);

// Integrate the geodesics ODE
std::vector<SPoint2> geodesicBetweenTwoPoints(const double x0[2], const double x1[2], const MetricFunctionPtr getMetric,
	const int maxIter, const double tol, const double ds, double &length, double &error);

std::vector<SPoint2> geodesicBetweenTwoPoints_arrayPtr(const double x0[2], const double x1[2], const MetricFunctionPtrArray getMetric,
	const double maxLength, const int maxIter, const double tol, const double ds, double &length, double &error);

// Minimize length functional
void geodesicBetweenTwoPoints_arrayPtr_minimizeLength(std::vector<Vertex> &vertices, 
	const MetricFunctionPtrArray getMetric, double &length);

// Compute the signed distance from a geodesic to a point
double distGeodesicToPoint(const double x[2], const std::vector<SPoint2> &geodesic);

#endif