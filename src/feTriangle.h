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
  virtual bool xyz2uvw(double xyz[3], double uvw[3]) = 0;

  // Check if physical point xyz is inside triangle
  // before computing the uvw coordinates.
  virtual bool isInside(double xyz[3], double tol) = 0;

  // Reference element is the same for P1 and P2 triangles
  bool isInside(double u, double v, double w)
  {
    double tol = 1e-5;
    if(u < (-tol) || v < (-tol) || u > ((1. + tol) - v) || fabs(w) > tol) return false;
    return true;
  }
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
  bool isInside(double xyz[3], double tol);
  bool xyz2uvw(double xyz[3], double uvw[3]);
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
  bool isInside(double xyz[3], double tol);
  bool xyz2uvw(double xyz[3], double uvw[3]);
};

#endif