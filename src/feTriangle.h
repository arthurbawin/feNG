#ifndef _FETRIANGLE_
#define _FETRIANGLE_

#include "feNG.h"
#include "feVertex.h"
#include "feEdge.h"

// To sort 3 integers
static void sort3(int *d) {
#define SWAP(x, y)                                                             \
  if(d[y] < d[x]) {                                                            \
    int tmp = d[x];                                                         \
    d[x] = d[y];                                                               \
    d[y] = tmp;                                                                \
  }
  SWAP(0, 1);
  SWAP(1, 2);
  SWAP(0, 1);
#undef SWAP
}

// Adapted from Gmsh's MTriangle

class Triangle
{
protected:
  int _tag; // Global (unique) tag in the mesh
  int _localTag; // Local tag on its Physical Entity
  int _pTag;
  Vertex *_vlin[3]; // The vertices of the underlying linear triangle to a linear or high-order triangle

  // Boundary is made of 3 edges and their orientation (+1 or -1 for reversed)
  const Edge* _edges[3];
  int _edgesOrientation[3];

public:
  // The former constructor, which may be removed
  Triangle(int tag = -1, int localTag = -1, int pTag = -1)
    : _tag(tag), _localTag(localTag), _pTag(pTag){}
  //
  // Create a generic (topological) triangle by providing its 3 vertices and a unique tag
  // The derived triangle can be linear, in which case they coincide,
  // or high-order, in which case this base class is the underlying linear triangle.
  // Used to create the Edges as boundary entities.
  //
  Triangle(Vertex *v0, Vertex *v1, Vertex *v2, int tag = -1, int localTag = -1, int pTag = -1)
    : _tag(tag), _localTag(localTag), _pTag(pTag)
    {
      _vlin[0] = v0;
      _vlin[1] = v1;
      _vlin[2] = v2;
    }
  ~Triangle() {}

  // Tag of the triangle
  int getTag() const { return _tag; }
  // Tag of the linear vertices
  int getTag(int i) const { return _vlin[i]->getTag(); }
  int getLocalTag() { return _localTag; }
  int getPhysicalTag() { return _pTag; }
  const Edge *getEdge(int num) const { return _edges[num]; }
  int getEdgeOrientation(int num) const { return _edgesOrientation[num]; }
  virtual int getNumVertices() const { return 3; };
  virtual Vertex *getVertex(int num) { return _vlin[num]; };

  // Compute the reference (u,v) coordinates and returns true on success.
  // For P1 triangle simply invert the linear transformation, always a success
  // if the element is valid.
  // For P2 triangle, invert quadratic transformation with a Newton-Raphson method,
  // success if converged to the prescribed tolerance.
  // The tolerance used here should probably be smaller than the tolerance used
  // to check if the (u,v) coordinates are inside the reference triangle (isInside)
  virtual bool xyz2uvw(double xyz[3], double uvw[3], double tol = 1e-5);

  // Reference element is the same for linear and high-order triangles
  bool isInsideReference(double u, double v, double w, double tol)
  {
    if(u < (-tol) || v < (-tol) || u > ((1. + tol) - v) || fabs(w) > tol) return false;
    return true;
  }

  // Check if physical point xyz is inside triangle.
  // For P1 triangle, we compute the uvw coordinates first to check, although
  // this is not required (could use half-plane check for example).
  // For P2 triangle, we check the orientation of the point w.r.t. the 
  // parabolic edges using their implicit equation.
  // THERE ARE STILL SOME ISSUES, AND ITS POSSIBLY DEPRECATED BECAUSE
  // JUST CHECKING THE NEWTON RESULTS SEEMS ROBUST ENOUGH
  // virtual bool isInsidePhysical(double xyz[3], double tol);

  virtual double sliverness();

  // Create the topological edges of the boundary of this triangle.
  // New edges (if any) are added to the set of edges.
  // Returns the number of edges after creation of the boundary,
  // that is, the size of meshEdges.
  void createBoundary(std::set<Edge, EdgeLessThan> &meshEdges, int &numEdges);
};

struct TriangleLessThan {
  bool operator()(const Triangle &t0, const Triangle &t1) const
  {
    int tags0[3] = {t0.getTag(0), t0.getTag(1), t0.getTag(2)};
    int tags1[3] = {t1.getTag(0), t1.getTag(1), t1.getTag(2)};
    sort3(tags0);
    sort3(tags1);
    int diffMin = tags0[0] - tags1[0];
    int diffMid = tags0[1] - tags1[1];
    int diffMax = tags0[2] - tags1[2];
    if(diffMin < 0) return true;
    if(diffMin > 0) return false;
    if(diffMid < 0) return true;
    if(diffMid > 0) return false;
    if(diffMax < 0) return true;
    if(diffMax > 0) return false;
    return false;
  }
};

class TriangleP2 : public Triangle
{
protected:
  // FIXME: Redundant, should only store the high-order vertices
  Vertex *_v[6];

  // Localization falls back to P1 case if the triangle is actually linear
  double _tolP2Straight = 1e-10;
  bool _isP2Straight = false;

public:
  TriangleP2(Vertex *v0, Vertex *v1, Vertex *v2, Vertex *v3, Vertex *v4, Vertex *v5,
             int tag = -1, int localTag = -1, int pTag = -1)
    : Triangle(v0, v1, v2, tag, localTag, pTag)
  {
    _v[0] = v0;
    _v[1] = v1;
    _v[2] = v2;
    _v[3] = v3;
    _v[4] = v4;
    _v[5] = v5;

    // Check if triangle is actually linear
    double dist1 = sqrt( (v3->x() - (v0->x()+v1->x())/2.) * (v3->x() - (v0->x()+v1->x())/2.) + (v3->y() - (v0->y()+v1->y())/2.) * (v3->y() - (v0->y()+v1->y())/2.) );
    double dist2 = sqrt( (v4->x() - (v1->x()+v2->x())/2.) * (v4->x() - (v1->x()+v2->x())/2.) + (v4->y() - (v1->y()+v2->y())/2.) * (v4->y() - (v1->y()+v2->y())/2.) );
    double dist3 = sqrt( (v5->x() - (v2->x()+v0->x())/2.) * (v5->x() - (v2->x()+v0->x())/2.) + (v5->y() - (v2->y()+v0->y())/2.) * (v5->y() - (v2->y()+v0->y())/2.) );
    if(dist1 < _tolP2Straight && dist2 < _tolP2Straight && dist3 < _tolP2Straight) {
      _isP2Straight = true;
    }
  }

  int getNumVertices(){ return 6; };
  Vertex *getVertex(int num) { return _v[num]; }
  bool isInsidePhysical(double xyz[3], double tol);
  bool xyz2uvw(double xyz[3], double uvw[3], double tol = 1e-5);
  double sliverness() { return -1.; };
};

// bool isInConvexRegionOfImplicitParabola(const double implicitParameters[6], double xyz[3]);

#endif