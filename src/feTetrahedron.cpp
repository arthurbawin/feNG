#include "feMessage.h"
#include "feTetrahedron.h"
#include <numeric>
#include <queue>

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
static double det3x3(double mat[3][3])
{
  return (mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
          mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
          mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]));
}

static int sys3x3(double mat[3][3], double b[3], double res[3], double *det)
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
  UNUSED(tol);
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

//
// Coefficients w.r.t. first barycentric coordinate of the inner reference
// vertices numbering in spiral on a high-order Lagrange triangle
//
void createSpiral(std::vector<int> &res, int n, int minVal, int numMinVal)
{
  int N = n-2;
  if(N <= 0) return;
  for(int n1 = N; n1 > minVal-1; --n1) {
    res.push_back(n1);
  }
  for(int n1 = 0; n1 < numMinVal; ++n1) {
    res.push_back(minVal);
  }
  for(int n1 = minVal; n1 < fmax(0,N); ++n1) {
    res.push_back(n1);
  }
  createSpiral(res, n-2, minVal+1, numMinVal-3);
}

static int sum4Digits(int val) {
  int v1 =  val / 1000;
  int v2 = (val % 1000) / 100;
  int v3 = (val % 100)  / 10;
  int v4 =  val % 10;
  return v1+v2+v3+v4;
}

//
// This is slow for very high polynomial degrees n (>12)
//
void getTetLagrangeBarycentricCoord(int n, std::vector<double> &barycentricCoord)
{
  barycentricCoord.clear();
  barycentricCoord.reserve(4*(n+1)*(n+2)*(n+3)/6);

  // Add the 4 vertices
  for(int i = 0; i < 4; ++i) {
    double pt[4] = {0.,0.,0.,0.};
    pt[i] = (double) n;
    barycentricCoord.insert(barycentricCoord.end(), std::begin(pt), std::end(pt));
  }

  // Edges: add n-1 reference vertices per edge
  if(n >= 2) {
    for(int i = 0; i < 6; ++i) {
      for(int j = 0; j < n-1; ++j) {
        double pt[4] = {0.,0.,0.,0.};
        for(int k = 0; k < 4; ++k) {
          pt[k] = (n-(j+1))/ (double) n * barycentricCoord[4*_edges[i][0]+k] +
                     (j+1) / (double) n * barycentricCoord[4*_edges[i][1]+k];
        }
        barycentricCoord.insert(barycentricCoord.end(), std::begin(pt), std::end(pt));
      }
    }
  }

  // Faces: add max(0, (n-1)(n-2)/2) reference vertices per face
  if(n >= 3) {
    // Apply spiral pattern
    std::vector<int> spiralCoeff;
    createSpiral(spiralCoeff, n, 1, n-4);

    // Determine number of spiral layers
    int n1 = fmax(0, 3*n-9); //3n-9 is the size of a layer of side (n-2) ( = 3(n-2)-3)
    int nLayers = 0;
    int N = n-2;
    std::vector<int> layerLength = {n1};
    while(n1 > 0) {
      nLayers++;
      N -= 3;
      n1 = fmax(0, 3*N-3); // 3N-3 is the size of a layer of side N
      layerLength.push_back(n1);
    }
    // int sum = std::reduce(layerLength.begin(),layerLength.end());
    // feInfo("There are %d layers - total length = %d - total length of layers = %d", nLayers, spiralCoeff.size(), sum);

    for(int i = 0; i < 4; ++i)
    {
      // Shift coefficients to the right by chunks of an outer spiral?
      std::vector<int> c2(spiralCoeff);
      int start = 0, end = 0, shift = fmax(0,n-3);
      for(int j = 0; j < nLayers; ++j) {
        end += layerLength[j];
        std::rotate(c2.begin()+start,
                    c2.begin()+end-shift,
                    c2.begin()+end);
        start = end;
        shift = fmax(0, shift-3);
      }
      std::vector<int> c3(c2);
      start = end = 0;
      shift = fmax(0,n-3);
      for(int j = 0; j < nLayers; ++j) {
        end += layerLength[j];
        std::rotate(c3.begin()+start,
                    c3.begin()+end-shift,
                    c3.begin()+end);
        start = end;
        shift = fmax(0, shift-3);
      }
      
      std::vector<std::vector<int>> all(4, std::vector<int>(spiralCoeff.size(), 0));
      all[_faces[i][0]] = spiralCoeff;
      all[_faces[i][1]] = c2;
      all[_faces[i][2]] = c3;

      for(size_t k = 0; k < spiralCoeff.size(); ++k) {
        // feInfo("%d - %d - %d - %d - sum = %d", all[0][k], all[1][k], all[2][k], all[3][k], all[0][k]+all[1][k]+all[2][k]+all[3][k]);
        double pt[4] = {(double) all[0][k], (double) all[1][k], (double) all[2][k], (double) all[3][k]};
        barycentricCoord.insert(barycentricCoord.end(), std::begin(pt), std::end(pt));
      }
    }
  }

  // Volume reference vertices.
  // The order does not matter for those vertices,
  // as they are not shared between two tets and there are
  // no orientation issues during the DOF numbering.
  // We just need to generate the complete list of barycentric
  // components in any order.
  if(n >= 4) {
    if(n < 13) {
      // If n < 13, we can represent a point in barycentric coordinates
      // using an integer with 4 digits which sum to up to 12,
      // this is faster than using vectors and sets.
      std::vector<int> allPts;
      std::queue<int> q;
      int ini = 1111;
      q.push(ini);
      while(!q.empty()) {
        int ft = q.front();
        int sum = sum4Digits(ft);
        if(sum > n) {
          // Do nothing and remove point
        } else if(sum == n) {
          // Keep point
          allPts.push_back(ft);
        } else {
          // Add all possible increments to the queue
          q.push(ft+1000);
          q.push(ft+100);
          q.push(ft+10);
          q.push(ft+1);
        }
        q.pop();
      }

      std::sort(allPts.begin(), allPts.end());
      allPts.erase(std::unique(allPts.begin(), allPts.end()), allPts.end());

      for(auto &val : allPts) {
        int v1 =  val / 1000;
        int v2 = (val % 1000) / 100;
        int v3 = (val % 100)  / 10;
        int v4 =  val % 10; 
        barycentricCoord.push_back(v1);
        barycentricCoord.push_back(v2);
        barycentricCoord.push_back(v3);
        barycentricCoord.push_back(v4);
      }
    } else {
      feWarning("Constructing reference vertices of degree %d with a set of vectors, which is slow.", n);
      std::set<std::vector<int>> allPts;
      std::queue<std::vector<int>> q;
      std::vector<int> pt = {1,1,1,1};
      q.push(pt);
      while(!q.empty()) {
        std::vector<int> &ft = q.front();
        int sum = std::accumulate(ft.begin(), ft.end(), 0);
        if(sum > n) {
          // Do nothing and remove point
        } else if(sum == n) {
          // Keep point
          allPts.insert(ft);
        } else {
          // Add all possible increments to the queue
          for(int i = 0; i < 4; ++i) {
            std::vector<int> pt2 = ft;
            pt2[i]++;
            q.push(pt2);
          }
        }
        q.pop();
      }

      for(auto &val : allPts) {
        barycentricCoord.insert(barycentricCoord.end(),val.begin(), val.end());
      }
    }
  }
}