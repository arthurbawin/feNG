#ifndef _FETRIANGLE_
#define _FETRIANGLE_

#include "feVertex.h"

#include <cmath>

// Adapted from Gmsh's MTriangle

class Triangle
{
protected:
  Vertex *_v[3];
  int _tag;
  int _pTag;

public:
  Triangle(Vertex *v0, Vertex *v1, Vertex *v2, int tag = -1, int pTag = -1) : _tag(tag), _pTag(pTag)
  {
    _v[0] = v0;
    _v[1] = v1;
    _v[2] = v2;
  }
  ~Triangle() {}

  int getTag() { return _tag; }
  int getPhysicalTag() { return _pTag; }
  Vertex *getVertex(int num) { return _v[num]; }
  void xyz2uvw(double xyz[3], double uvw[3]);
  bool isInside(double u, double v, double w)
  {
    double tol = 1e-6;
    if(u < (-tol) || v < (-tol) || u > ((1. + tol) - v) || fabs(w) > tol) return false;
    return true;
  }
};

#endif