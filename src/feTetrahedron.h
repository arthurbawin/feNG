#ifndef _FETETRAHEDRON_
#define _FETETRAHEDRON_

#include "feNG.h"
#include "feVertex.h"
#include "feTriangle.h"

// To sort 4 integers
static void sort4(int *d) {
#define SWAP(x, y)                                                             \
  if(d[y] < d[x]) {                                                            \
    int tmp = d[x];                                                         \
    d[x] = d[y];                                                               \
    d[y] = tmp;                                                                \
  }
  SWAP(0, 2);
  SWAP(1, 3);
  SWAP(0, 1);
  SWAP(2, 3);
  SWAP(1, 2);
#undef SWAP
}

// Adapted from Gmsh's MTetrahedron

class Tetrahedron
{
protected:
  int _tag; // Global (unique) tag in the mesh
  int _localTag; // Local tag on its Physical Entity
  int _pTag;
  Vertex *_vlin[4]; // The vertices of the underlying linear tetrahedron to a linear or high-order tet

  // Boundary is made of 4 triangles and their orientation (+1 or -1 for reversed)
  const Triangle* _facets[4];
  int _facetsOrientation[4];

public:
  // Create a generic (topological) tetrahedron by providing its 4 vertices and a unique tag
  // The derived tetrahedron can be linear, in which case they coincide,
  // or high-order, in which case this base class is the underlying linear tetrahedron.
  // Used to create the Triangles as boundary entities.
  //
  Tetrahedron(Vertex *v0, Vertex *v1, Vertex *v2, Vertex *v3, int tag = -1, int localTag = -1, int pTag = -1)
    : _tag(tag), _localTag(localTag), _pTag(pTag)
    {
      _vlin[0] = v0;
      _vlin[1] = v1;
      _vlin[2] = v2;
      _vlin[3] = v3;
    }
  ~Tetrahedron() {}

  // Tag of the tetrahedron
  int getTag() const { return _tag; }
  // Tag of the linear vertices
  int getTag(int i) const { return _vlin[i]->getTag(); }
  int getLocalTag() const { return _localTag; }
  int getPhysicalTag() const { return _pTag; }
  const Triangle *getFacet(int num) const { return _facets[num]; }
  int getFacetOrientation(int num) const { return _facetsOrientation[num]; }
  virtual int getNumVertices() const { return 4; };
  virtual Vertex *getVertex(int num) { return _vlin[num]; };

  // Compute the reference (u,v,w) coordinates and returns true on success.
  // For P1 tet simply invert the linear transformation, always a success
  // if the element is valid.
  // The tolerance used here should probably be smaller than the tolerance used
  // to check if the uvw coordinates are inside the reference triangle (isInsideReference)
  virtual bool xyz2uvw(double xyz[3], double uvw[3], double tol = 1e-5);

  // Check if reference point uvw is inside the reference tet
  bool isInsideReference(double u, double v, double w, double tol = 1e-4)
  {
    if(u < (-tol) || v < (-tol) || u > ((1. + tol) - v) || fabs(w) > tol) return false;
    return true;
  }

  virtual double sliverness();

  // From Gmsh's MTetrahedron.h.
  // This numbering is given by the boundary operator del = (-1)^i * [v0, ..., \hat{vi}, ..., vn],
  // where \hat{vi} means vi is omitted on the facet.
  static int faces_tetra(const int face, const int vert)
  {
    static const int f[4][3] = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {3, 1, 2}};
    return f[face][vert];
  }

  // Create the topological triangles of the boundary of this tetrahedron.
  // New triangles (if any) are added to the set of triangles.
  // Returns the number of triangles after creation of the boundary,
  // that is, the size of meshTriangles.
  void createBoundary(std::set<Triangle, TriangleLessThan> &meshTriangles, int &numTriangles,
    std::set<Edge, EdgeLessThan> &meshEdges, int &numEdges);
};

struct TetrahedronLessThan {
  bool operator()(const Tetrahedron &t0, const Tetrahedron &t1) const
  {
    int tags0[4] = {t0.getTag(0), t0.getTag(1), t0.getTag(2), t0.getTag(3)};
    int tags1[4] = {t1.getTag(0), t1.getTag(1), t1.getTag(2), t1.getTag(3)};
    sort4(tags0);
    sort4(tags1);
    int diffMin = tags0[0] - tags1[0];
    int diffMid1 = tags0[1] - tags1[1];
    int diffMid2 = tags0[2] - tags1[2];
    int diffMax = tags0[3] - tags1[3];
    if(diffMin < 0) return true;
    if(diffMin > 0) return false;
    if(diffMid1 < 0) return true;
    if(diffMid1 > 0) return false;
    if(diffMid2 < 0) return true;
    if(diffMid2 > 0) return false;
    if(diffMax < 0) return true;
    if(diffMax > 0) return false;
    return false;
  }
};

#endif