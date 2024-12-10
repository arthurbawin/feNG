#include "feTetrahedron.h"

void Tetrahedron::createBoundary(std::set<Triangle, TriangleLessThan> &meshTriangles, int &numTriangles,
  std::set<Edge, EdgeLessThan> &meshEdges, int &numEdges)
{
  std::pair<std::set<Triangle, TriangleLessThan>::iterator, bool> ret;

  for(int i = 0; i < 4; ++i)
  {
    Vertex *v0 = _vlin[faces_tetra(i,0)];
    Vertex *v1 = _vlin[faces_tetra(i,1)];
    Vertex *v2 = _vlin[faces_tetra(i,2)];

    Triangle tri(v0, v1, v2, numTriangles);
    tri.createBoundary(meshEdges, numEdges);
    ret = meshTriangles.insert(tri);

    // Iterator in ret.first points to the added facet, or to equivalent facet already in set
    _facets[i] = &(*ret.first);

    if(ret.second) { // Facet was added to the set
      _facetsOrientation[i] = 1;
      numTriangles++;
    } else { // Facet is already in the set : orientation is negative
      _facetsOrientation[i] = -1;
    }
  }
}

// From Gmsh's Numeric.cpp
double det3x3(double mat[3][3])
{
  return (mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
          mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
          mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]));
}

int sys3x3(double mat[3][3], double b[3], double res[3], double *det)
{
  double ud;
  int i;

  *det = det3x3(mat);

  if(*det == 0.0) {
    res[0] = res[1] = res[2] = 0.0;
    return (0);
  }

  ud = 1. / (*det);

  res[0] = b[0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
           mat[0][1] * (b[1] * mat[2][2] - mat[1][2] * b[2]) +
           mat[0][2] * (b[1] * mat[2][1] - mat[1][1] * b[2]);

  res[1] = mat[0][0] * (b[1] * mat[2][2] - mat[1][2] * b[2]) -
           b[0] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
           mat[0][2] * (mat[1][0] * b[2] - b[1] * mat[2][0]);

  res[2] = mat[0][0] * (mat[1][1] * b[2] - b[1] * mat[2][1]) -
           mat[0][1] * (mat[1][0] * b[2] - b[1] * mat[2][0]) +
           b[0] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);

  for(i = 0; i < 3; i++) res[i] *= ud;
  return (1);
}

bool Tetrahedron::xyz2uvw(double xyz[3], double uvw[3], double tol)
{
  double mat[3][3], b[3], det;
  mat[0][0] = _vlin[1]->x() - _vlin[0]->x();
  mat[0][1] = _vlin[2]->x() - _vlin[0]->x();
  mat[0][2] = _vlin[3]->x() - _vlin[0]->x();
  mat[1][0] = _vlin[1]->y() - _vlin[0]->y();
  mat[1][1] = _vlin[2]->y() - _vlin[0]->y();
  mat[1][2] = _vlin[3]->y() - _vlin[0]->y();
  mat[2][0] = _vlin[1]->z() - _vlin[0]->z();
  mat[2][1] = _vlin[2]->z() - _vlin[0]->z();
  mat[2][2] = _vlin[3]->z() - _vlin[0]->z();
  b[0] = xyz[0] - _vlin[0]->x();
  b[1] = xyz[1] - _vlin[0]->y();
  b[2] = xyz[2] - _vlin[0]->z();
  sys3x3(mat, b, uvw, &det);
  return true;
}

double Tetrahedron::sliverness()
{
  printf("Implement sliverness for TetrahedronP1!\n");
  return -1.;
}