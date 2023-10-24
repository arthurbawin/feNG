#ifndef _FETRIANGLE_
#define _FETRIANGLE_

#include "feVertex.h"

#include <cmath>

// Adapted from Gmsh's MTriangle

class Triangle
{
protected:
  int _tag; // Global (unique) tag in the mesh
  int _localTag; // Local tag on its Physical Entity
  int _pTag;

public:
  Triangle(int tag = -1, int localTag = -1, int pTag = -1)
    : _tag(tag), _localTag(localTag), _pTag(pTag){}
  ~Triangle() {}

  int getTag() { return _tag; }
  int getLocalTag() { return _localTag; }
  int getPhysicalTag() { return _pTag; }
  virtual int getNumVertices() = 0;
  virtual Vertex *getVertex(int num) = 0;

  // Compute the reference (u,v) coordinates and returns true on success.
  // For P1 triangle simply invert the linear transformation, always a success
  // if the element is valid.
  // For P2 triangle, invert quadratic transformation with a Newton-Raphson method,
  // success if converged to the prescribed tolerance.
  // The tolerance used here should probably be smaller than the tolerance used
  // to check if the (u,v) coordinates are inside the reference triangle (isInside)
  virtual bool xyz2uvw(double xyz[3], double uvw[3], double tol = 1e-5) = 0;

  // Check if physical point xyz is inside triangle.
  // For P1 triangle, we compute the uvw coordinates first to check, although
  // this is not required (could use half-plane check for example).
  // For P2 triangle, we check the orientation of the point w.r.t. the 
  // parabolic edges using their implicit equation.
  // THERE ARE STILL SOME ISSUES, AND ITS POSSIBLY DEPRECATED BECAUSE
  // JUST CHECKING THE NEWTON RESULTS SEEMS ROBUST ENOUGH
  virtual bool isInsidePhysical(double xyz[3], double tol) = 0;

  virtual double sliverness() = 0;

  // Reference element is the same for P1 and P2 triangles
  bool isInsideReference(double u, double v, double w, double tol)
  {
    if(u < (-tol) || v < (-tol) || u > ((1. + tol) - v) || fabs(w) > tol) return false;
    return true;
  }

protected:
  // Compute the reference (u,v) coordinates for P1 triangle.
  // If derived class is P1 triangle, it is the true (u,v) coordinates.
  // If dervied class is P2 traignel, it is the (u,v) coordinates
  // of the P1 triangle formed by the 3 first vertices.
  void getP1ReferenceCoordinates(double xyz[3], double uvw[3], Vertex **v);
};

class TriangleP1 : public Triangle
{
protected:
  Vertex *_v[3];

public:
  TriangleP1(Vertex *v0, Vertex *v1, Vertex *v2, int tag = -1, int localTag = -1, int pTag = -1)
    : Triangle(tag, localTag, pTag)
  {
    _v[0] = v0;
    _v[1] = v1;
    _v[2] = v2;
  }

  int getNumVertices(){ return 3; };
  Vertex *getVertex(int num) { return _v[num]; }
  bool isInsidePhysical(double xyz[3], double tol);
  bool xyz2uvw(double xyz[3], double uvw[3], double tol = 1e-5);
  double sliverness();
};

class TriangleP2 : public Triangle
{
protected:
  Vertex *_v[6];

public:
  TriangleP2(Vertex *v0, Vertex *v1, Vertex *v2, Vertex *v3, Vertex *v4, Vertex *v5,
             int tag = -1, int localTag = -1, int pTag = -1)
    : Triangle(tag, localTag, pTag)
  {
    _v[0] = v0;
    _v[1] = v1;
    _v[2] = v2;
    _v[3] = v3;
    _v[4] = v4;
    _v[5] = v5;
  }

  int getNumVertices(){ return 6; };
  Vertex *getVertex(int num) { return _v[num]; }
  bool isInsidePhysical(double xyz[3], double tol);
  bool xyz2uvw(double xyz[3], double uvw[3], double tol = 1e-5);
  double sliverness() { return -1.; };
};

bool isInConvexRegionOfImplicitParabola(const double implicitParameters[6], double xyz[3]);

#endif