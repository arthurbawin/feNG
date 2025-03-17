#include "feMetricTools.h"
#include "feNewRecovery.h"

#if defined(HAVE_GMSH)
#include "gmsh.h"
#endif
#if defined(HAVE_SOPLEX)
#include "soplex.h"
using namespace soplex;
#endif

#include "fullMatrix.h"
#include "ellipseToolbox.h"
#include "../contrib/root_finder/root_finder.hpp"
#include <Eigen/Eigen>
#include <unsupported/Eigen/MatrixFunctions>

extern int FE_VERBOSE;

double evaluateFieldFromRecoveryCallback(int indexDerivative, void *recUserPtr, double *x)
{
  return ((feRecovery *)recUserPtr)->evalDerivative(indexDerivative, x);
}

template <class MetricType>
void drawSingleEllipse(FILE *file, const double *x, MetricType &M, double sizeFactor, int nPoints)
{
  int nt = nPoints;
  std::vector<double> xP(nt, 0.);
  std::vector<double> yP(nt, 0.);

  getEllipsePoints(sizeFactor * M(0, 0), sizeFactor * 2.0 * M(0, 1), sizeFactor * M(1, 1), 
  	x[0], x[1], xP, yP);
  for(int i = 0; i < nt; ++i) {
    if(i != nt - 1) {
      fprintf(file, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[i], yP[i], 0.,
              xP[i + 1], yP[i + 1], 0., 1, 1);
    } else {
      fprintf(file, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[i], yP[i], 0., xP[0],
              yP[0], 0., 1, 1);
    }
  }
}

// Explicit instantiation
template void drawSingleEllipse<MetricTensor>(FILE *file, const double *x, MetricTensor &M, double sizeFactor, int nPoints);
template void drawSingleEllipse<SMetric3>(FILE *file, const double *x, SMetric3 &M, double sizeFactor, int nPoints);

template <class MetricType>
void drawEllipsoids(const std::string &posFile, std::map<int, MetricType> &metrics,
                    const std::vector<std::size_t> &nodeTags, const std::vector<double> &coord,
                    double sizeFactor, int nPoints)
{
  FILE *f = fopen(posFile.data(), "w");
  fprintf(f, "View \"%s\"{\n", posFile.data());

  double factor = sizeFactor;
  int nt = nPoints;
  std::vector<double> xP(nt, 0.);
  std::vector<double> yP(nt, 0.);

  for(size_t i = 0; i < nodeTags.size(); i++) {
    MetricType &m = metrics.at(nodeTags[i]);
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];
    bool OK = getEllipsePoints(factor * m(0, 0), factor * 2.0 * m(0, 1), factor * m(1, 1), x, y, xP, yP);
    if(!OK) {
      feWarning("Metric: m11 = %+-1.6e - m12 = %+-1.6e - m21 = %+-1.6e - m22 = %+-1.6e - det(M) = %+-1.6e vs %+-1.6e", 
        m(0,0), m(0,1), m(1,0), m(1,1), m(0,0) * m(1,1) - m(0,1)*m(1,0), m.determinant());
    }
    for(int j = 0; j < nt; ++j) {
      if(j != nt - 1) {
        fprintf(f, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0.,
                xP[j + 1], yP[j + 1], 0., 1, 1);
      } else {
        fprintf(f, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0., xP[0],
                yP[0], 0., 1, 1);
      }
    }
  }
  fprintf(f, "};");
  fclose(f);

  // Also print the metric tensors to text file
  std::string textfile = posFile + ".txt";
  f = fopen(textfile.data(), "w");
  for(size_t i = 0; i < nodeTags.size(); i++) {
    MetricType &m = metrics[nodeTags[i]];
    fprintf(f, "%+-1.16e \t %+-1.16e \t %+-1.16e \n", m(0,0), m(0,1), m(1,1));
  }
  fclose(f);
}

// Explicit instantiation
template void drawEllipsoids<MetricTensor>(const std::string &posFile,
                                           std::map<int, MetricTensor> &metrics,
                                           const std::vector<std::size_t> &nodeTags,
                                           const std::vector<double> &coord, double sizeFactor,
                                           int nPoints);

SVector3 E_X(1., 0., 0.);
SVector3 E_Y(0., 1., 0.);
SVector3 E_Z(0., 0., 1.);

// Metrics are modified if || M - M_grown || < tolerance_gradation
double tolerance_gradation = 1e-1;
// If fabs(non-diag) is below tolerance_diagonal, matrix is diag
double tolerance_diagonal = 1e-4;

double matNorm2(const SMetric3 &m)
{
  double sqr = 0;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      sqr += m(i, j) * m(i, j);
    }
  }
  return sqrt(sqr);
}

double matNorm2(const SMetric3 &m1, const SMetric3 &m2)
{
  double sqr = 0;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      sqr += (m1(i, j) - m2(i, j)) * (m1(i, j) - m2(i, j));
    }
  }
  return sqrt(sqr);
}

SMetric3 intersectionReductionSimultaneeExplicite(const SMetric3 &m1, const SMetric3 &m2)
{
  double a1 = m1(0, 0), b1 = m1(0, 1), c1 = m1(1, 1);
  double a2 = m2(0, 0), b2 = m2(0, 1), c2 = m2(1, 1);
  double v00, v01, vTInv00, vTInv01, vTInv10, vTInv11;
  // Metrics are identical
  double relDiff = matNorm2(m1, m2) / matNorm2(m1);
  if(matNorm2(m1, m2) / matNorm2(m1) < 1e-2) {
    // std::cout<<"Identical matrix : returning m1"<<std::endl;
    return m2;
  }

  if(fabs(b1) < tolerance_diagonal && fabs(b2) < tolerance_diagonal) {
    // Both matrices are diagonal
    return SMetric3(fmax(a1, a2), fmax(c1, c2), 1., E_X, E_Y, E_Z);

  } else if(fabs(b1) < tolerance_diagonal) {
    // m1 is diagonal
    v00 = (a1 * c2 + a2 * c1 -
           sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 + 4 * a1 * b2 * b2 * c1 +
                a2 * a2 * c1 * c1)) /
            (2 * a1 * b2) -
          c2 / b2;
    v01 = (a1 * c2 + a2 * c1 +
           sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 + 4 * a1 * b2 * b2 * c1 +
                a2 * a2 * c1 * c1)) /
            (2 * a1 * b2) -
          c2 / b2;
    SVector3 v0(v00, 1., 0.); // v0.normalize();
    SVector3 v1(v01, 1., 0.); // v1.normalize();
    double l0 = fmax(dot(v0, m1, v0), dot(v0, m2, v0));
    double l1 = fmax(dot(v1, m1, v1), dot(v1, m2, v1));
    vTInv00 = -(a1 * b2) / sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 + 4 * a1 * b2 * b2 * c1 +
                                a2 * a2 * c1 * c1);
    vTInv01 =
      (sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 + 4 * a1 * b2 * b2 * c1 + a2 * a2 * c1 * c1) -
       a1 * c2 + a2 * c1) /
      (2 *
       sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 + 4 * a1 * b2 * b2 * c1 + a2 * a2 * c1 * c1));
    vTInv10 = (a1 * b2) / sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 + 4 * a1 * b2 * b2 * c1 +
                               a2 * a2 * c1 * c1);
    vTInv11 =
      (sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 + 4 * a1 * b2 * b2 * c1 + a2 * a2 * c1 * c1) +
       a1 * c2 - a2 * c1) /
      (2 *
       sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 + 4 * a1 * b2 * b2 * c1 + a2 * a2 * c1 * c1));
    if(isnan(vTInv00) || isnan(vTInv01) || isnan(vTInv10) || isnan(vTInv11) || isinf(vTInv00) ||
       isinf(vTInv01) || isinf(vTInv10) || isinf(vTInv11)) {
      std::cout << relDiff << std::endl;
      m1.print("1 - m1");
      m2.print("1 - m2");
      std::cout << "NAN alert : " << vTInv00 << " " << vTInv01 << " " << vTInv10 << " " << vTInv11
                << std::endl;
      exit(-1);
    }
    v0 = SVector3(vTInv00, vTInv01, 0.); // v0.normalize();
    v1 = SVector3(vTInv10, vTInv11, 0.); // v1.normalize();
    return SMetric3(l0, l1, 1., v0, v1, E_Z);

  } else if(fabs(b2) < tolerance_diagonal) {
    // m2 is diagonal
    v00 = (a1 * c2) / (a2 * b1) - (a1 * c2 + a2 * c1 +
                                   sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 +
                                        a2 * a2 * c1 * c1 + 4 * a2 * b1 * b1 * c2)) /
                                    (2 * a2 * b1);
    v01 = (a1 * c2) / (a2 * b1) - (a1 * c2 + a2 * c1 -
                                   sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 +
                                        a2 * a2 * c1 * c1 + 4 * a2 * b1 * b1 * c2)) /
                                    (2 * a2 * b1);
    SVector3 v0(v00, 1., 0.); // v0.normalize();
    SVector3 v1(v01, 1., 0.); // v1.normalize();
    double l0 = fmax(dot(v0, m1, v0), dot(v0, m2, v0));
    double l1 = fmax(dot(v1, m1, v1), dot(v1, m2, v1));
    vTInv00 = -(a2 * b1) / sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 + a2 * a2 * c1 * c1 +
                                4 * a2 * b1 * b1 * c2);
    vTInv01 =
      (sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 + a2 * a2 * c1 * c1 + 4 * a2 * b1 * b1 * c2) +
       a1 * c2 - a2 * c1) /
      (2 *
       sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 + a2 * a2 * c1 * c1 + 4 * a2 * b1 * b1 * c2));
    vTInv10 = (a2 * b1) / sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 + a2 * a2 * c1 * c1 +
                               4 * a2 * b1 * b1 * c2);
    vTInv11 =
      (sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 + a2 * a2 * c1 * c1 + 4 * a2 * b1 * b1 * c2) -
       a1 * c2 + a2 * c1) /
      (2 *
       sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 + a2 * a2 * c1 * c1 + 4 * a2 * b1 * b1 * c2));
    if(isnan(vTInv00) || isnan(vTInv01) || isnan(vTInv10) || isnan(vTInv11) || isinf(vTInv00) ||
       isinf(vTInv01) || isinf(vTInv10) || isinf(vTInv11)) {
      std::cout << relDiff << std::endl;
      m1.print("2 - m1");
      m2.print("2 - m2");
      std::cout << "NAN alert : " << vTInv00 << " " << vTInv01 << " " << vTInv10 << " " << vTInv11
                << std::endl;
      exit(-1);
    }
    v0 = SVector3(vTInv00, vTInv01, 0.); // v0.normalize();
    v1 = SVector3(vTInv10, vTInv11, 0.); // v1.normalize();
    return SMetric3(l0, l1, 1., v0, v1, E_Z);

  } else {
    // Check if metrics are multiple of one another
    if(fabs(a2 / a1 - c2 / c1) < 1e-3 && fabs(a2 / a1 - b2 / b1) < 1e-3) {
      return (a2 / a1 <= 1) ? m1 : m2;
    } else {
      v00 = (sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 - 4 * a1 * b1 * b2 * c2 +
                  4 * a1 * b2 * b2 * c1 + a2 * a2 * c1 * c1 + 4 * a2 * b1 * b1 * c2 -
                  4 * a2 * b1 * b2 * c1) +
             a1 * c2 + a2 * c1 - 2 * b1 * b2) /
              (2 * (a1 * b2 - a2 * b1)) -
            (a1 * c2 - b1 * b2) / (a1 * b2 - a2 * b1);
      v01 = -(a1 * c2 - b1 * b2) / (a1 * b2 - a2 * b1) -
            (sqrt(a1 * a1 * c2 * c2 - 2. * a1 * a2 * c1 * c2 - 4. * a1 * b1 * b2 * c2 +
                  4. * a1 * b2 * b2 * c1 + a2 * a2 * c1 * c1 + 4. * a2 * b1 * b1 * c2 -
                  4. * a2 * b1 * b2 * c1) -
             a1 * c2 - a2 * c1 + 2. * b1 * b2) /
              (2. * (a1 * b2 - a2 * b1));
      SVector3 v0(v00, 1., 0.); // v0.normalize();
      SVector3 v1(v01, 1., 0.); // v1.normalize();
      double l0 = fmax(dot(v0, m1, v0), dot(v0, m2, v0));
      double l1 = fmax(dot(v1, m1, v1), dot(v1, m2, v1));
      vTInv00 = (a1 * b2 - a2 * b1) /
                sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 - 4 * a1 * b1 * b2 * c2 +
                     4 * a1 * b2 * b2 * c1 + a2 * a2 * c1 * c1 + 4 * a2 * b1 * b1 * c2 -
                     4 * a2 * b1 * b2 * c1);
      vTInv01 = (sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 - 4 * a1 * b1 * b2 * c2 +
                      4 * a1 * b2 * b2 * c1 + a2 * a2 * c1 * c1 + 4 * a2 * b1 * b1 * c2 -
                      4 * a2 * b1 * b2 * c1) +
                 a1 * c2 - a2 * c1) /
                (2 * sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 - 4 * a1 * b1 * b2 * c2 +
                          4 * a1 * b2 * b2 * c1 + a2 * a2 * c1 * c1 + 4 * a2 * b1 * b1 * c2 -
                          4 * a2 * b1 * b2 * c1));
      vTInv10 = -(a1 * b2 - a2 * b1) /
                sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 - 4 * a1 * b1 * b2 * c2 +
                     4 * a1 * b2 * b2 * c1 + a2 * a2 * c1 * c1 + 4 * a2 * b1 * b1 * c2 -
                     4 * a2 * b1 * b2 * c1);
      vTInv11 = (sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 - 4 * a1 * b1 * b2 * c2 +
                      4 * a1 * b2 * b2 * c1 + a2 * a2 * c1 * c1 + 4 * a2 * b1 * b1 * c2 -
                      4 * a2 * b1 * b2 * c1) -
                 a1 * c2 + a2 * c1) /
                (2 * sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 - 4 * a1 * b1 * b2 * c2 +
                          4 * a1 * b2 * b2 * c1 + a2 * a2 * c1 * c1 + 4 * a2 * b1 * b1 * c2 -
                          4 * a2 * b1 * b2 * c1));
      if(isnan(vTInv00) || isnan(vTInv01) || isnan(vTInv10) || isnan(vTInv11) || isinf(vTInv00) ||
         isinf(vTInv01) || isinf(vTInv10) || isinf(vTInv11)) {
        std::cout << relDiff << std::endl;
        m1.print("3 - m1");
        m2.print("3 - m2");
        std::cout << "NAN alert : " << vTInv00 << " " << vTInv01 << " " << vTInv10 << " " << vTInv11
                  << std::endl;
        exit(-1);
      }
      v0 = SVector3(vTInv00, vTInv01, 0.); // v0.normalize();
      v1 = SVector3(vTInv10, vTInv11, 0.); // v1.normalize();
      return SMetric3(l0, l1, 1., v0, v1, E_Z);
    }
  }
}

SMetric3 gradationAlauzet(double grad, const SMetric3 &m, double x1, double y1, double x2,
                                 double y2)
{
  SVector3 pq(x2 - x1, y2 - y1, 0.);
  double eta = (1. + sqrt(dot(pq, m, pq)) * log(grad));
  eta = 1. / (eta * eta);
  SMetric3 tmp(m);
  tmp *= eta;
  tmp.set_m33(1.0);
  return tmp;
}

void gradationMetriques(double gradation, int nmax, std::vector<double> &/*coord1*/,
                        std::map<int, SMetric3> &metricsOnGmshModel)
{
#if defined(HAVE_GMSH)
  std::vector<int> elementTypes;
  std::vector<std::vector<std::size_t> > elementTags;
  std::vector<std::vector<std::size_t> > elemNodeTags;
  gmsh::model::mesh::getElements(elementTypes, elementTags, elemNodeTags, 2);

  // std::vector<std::pair<size_t,size_t>> edges;
  std::set<std::pair<size_t, size_t>, gmshEdgeLessThan> edges;

  std::vector<std::size_t> nodeTags;
  std::vector<double> coord;
  std::vector<double> nodecoord1;
  std::vector<double> nodecoord2;
  std::vector<double> nodecoord3;
  std::vector<double> parametricCoord;
  std::vector<double> nodeparametricCoord;
  gmsh::model::mesh::getNodesByElementType(elementTypes[0], nodeTags, coord, parametricCoord);

  size_t numNodesPerElem = (elementTypes[0] == 2) ? 3 : 6;

  for(auto val : elementTypes) feInfoCond(FE_VERBOSE > 0, "Gradation: Elements classified on dim 2 are of type %d", val);
  feInfoCond(FE_VERBOSE > 0, "Gradation: There are %d elements of dim 2 and type [0]", elementTags[0].size());
  feInfoCond(FE_VERBOSE > 0, "Gradation: There are %d nodes total in elem [0]", nodeTags.size());

  int nodedim, nodetag;

  FILE *fff = fopen("CHECK.pos", "w");
  fprintf(fff, "View \" CHECK \"{\n");

  std::pair<std::set<std::pair<size_t, size_t>, gmshEdgeLessThan>::iterator, bool> ret;

  for(size_t i = 0; i < elementTags[0].size(); i++) {
    for(size_t j = 0; j < 3; ++j) {

      if(numNodesPerElem == 3) {
        // P1 triangles
        // Add each edge
        size_t n0 = elemNodeTags[0][numNodesPerElem * i + j];
        size_t n1 = elemNodeTags[0][numNodesPerElem * i + (j+1) % numNodesPerElem];
        ret = edges.insert(std::make_pair(n0, n1));
      } else {
        // P2 triangles
        // Add the 6 semi-edges from a P1 vertex to a P2 midnode
        size_t n0 = elemNodeTags[0][numNodesPerElem * i + j];
        size_t n1 = elemNodeTags[0][numNodesPerElem * i + j+3];
        ret = edges.insert(std::make_pair(n0, n1));
      }

      // size_t n0 = elemNodeTags[0][numNodesPerElem * i + 0];
      // size_t n1 = elemNodeTags[0][numNodesPerElem * i + 1];
      // size_t n2 = elemNodeTags[0][numNodesPerElem * i + 2];
      // ret = edges.insert(std::make_pair(n0, n1));
      // ret = edges.insert(std::make_pair(n1, n2));
      // ret = edges.insert(std::make_pair(n2, n0));
      // gmsh::model::mesh::getNode(n0, nodecoord1, nodeparametricCoord, nodedim, nodetag);
      // gmsh::model::mesh::getNode(n1, nodecoord2, nodeparametricCoord, nodedim, nodetag);
      // gmsh::model::mesh::getNode(n2, nodecoord3, nodeparametricCoord, nodedim, nodetag);
      // double x1 = nodecoord1[0];
      // double y1 = nodecoord1[1];
      // double x2 = nodecoord2[0];
      // double y2 = nodecoord2[1];
      // double x3 = nodecoord3[0];
      // double y3 = nodecoord3[1];
      // fprintf(fff, "SL(%.16g,%.16g,0.,%.16g,%.16g,0.){%u, %u};\n", x1, y1, x2, y2, 1, 1);
      // fprintf(fff, "SL(%.16g,%.16g,0.,%.16g,%.16g,0.){%u, %u};\n", x2, y2, x3, y3, 1, 1);
      // fprintf(fff, "SL(%.16g,%.16g,0.,%.16g,%.16g,0.){%u, %u};\n", x3, y3, x1, y1, 1, 1);
    }
  }

  fprintf(fff, "};");
  fclose(fff);

  int nt = 30;
  std::vector<double> xP(nt, 0.);
  std::vector<double> yP(nt, 0.);

  SMetric3 mp, mq, mpq, mqp;
  double xp, yp, xq, yq;
  bool correction = true;
  int iter = 0;
  while(correction && iter < nmax) {
    correction = false;
    iter++;

    for(auto edge : edges) {
      size_t n = edge.first;
      size_t neigh = edge.second;
      mp = metricsOnGmshModel[n];
      mq = metricsOnGmshModel[neigh];

      gmsh::model::mesh::getNode(n, nodecoord1, nodeparametricCoord, nodedim, nodetag);
      gmsh::model::mesh::getNode(neigh, nodecoord2, nodeparametricCoord, nodedim, nodetag);

      xp = nodecoord1[0];
      yp = nodecoord1[1];
      xq = nodecoord2[0];
      yq = nodecoord2[1];

      // fprintf(fff, "SL(%.16g,%.16g,0.,%.16g,%.16g,0.){%u, %u};\n", xp, yp, xq, yq, 1, 1);
      // fprintf(fff, "SP(%.16g,%.16g,0.){%u};\n", xp, yp, n);
      // Span mp to neighbouring vertex
      mpq = gradationAlauzet(gradation, mp, xp, yp, xq, yq);
      if(fabs(mpq(2, 2) - 1.0) > 1e-4) {
        std::cout << "NONZERO z entry in mpq";
        exit(-1);
      }
      mpq = intersectionReductionSimultaneeExplicite(mq, mpq);
      if(fabs(mpq(2, 2) - 1.0) > 1e-4) {
        std::cout << "NONZERO z entry in mpq apres gradation";
        exit(-1);
      }

      // mpq = intersection_reductionSimultanee(mq, mpq);
      // metricsOnGmshModel[neigh] = mpq;

      if(matNorm2(metricsOnGmshModel[neigh], mpq) / matNorm2(mpq) > tolerance_gradation) {
        metricsOnGmshModel[neigh] = mpq;
        correction = true;
      }
      // Span mq to n
      mqp = gradationAlauzet(gradation, mq, xq, yq, xp, yp);
      if(fabs(mpq(2, 2) - 1.0) > 1e-4) {
        std::cout << "NONZERO z entry in mqp";
        exit(-1);
      }
      mqp = intersectionReductionSimultaneeExplicite(mp, mqp);
      if(fabs(mpq(2, 2) - 1.0) > 1e-4) {
        std::cout << "NONZERO z entry in mqp apres gradation";
        exit(-1);
      }

      // mqp = intersection_reductionSimultanee(mp, mqp);
      // metricsOnGmshModel[n] = mqp;

      if(matNorm2(metricsOnGmshModel[n], mqp) / matNorm2(mqp) > tolerance_gradation) {
        metricsOnGmshModel[n] = mqp;
        correction = true;
      }

      // mp = gradationAlauzet(gradation, mref, X0, Y0, xp, yp);

      // double factor = 1000.;
      // getEllipsePoints(factor * mp(0, 0), factor * 2.0 * mp(0, 1), factor * mp(1, 1), xp, yp, xP,
      // yP); for(int j = 0; j < nt; ++j) {
      //   if(j != nt - 1) {
      //     fprintf(fff, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
      //             0., xP[j + 1], yP[j + 1], 0., 1, 1);
      //   } else {
      //     fprintf(fff, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
      //             0., xP[0], yP[0], 0., 1, 1);
      //   }
      // }
    }

    feInfoCond(FE_VERBOSE > 0, "Gradation: Passe %d", iter);
    // correction = false;

    // fprintf(fff, "};");
    // fclose(fff);
  }
  // feInfo("%d coord size", coord.size());
  // feInfo("%d edges size", edges.size());
#else
  UNUSED(gradation, nmax, coord1, metricsOnGmshModel);
#endif
}

void smoothDirections(std::map<size_t, double> &C, std::map<size_t, double> &S, int nIter, double tol, bool plot)
{
#if defined(HAVE_GMSH)
  std::vector<int> elementTypes;
  std::vector<std::vector<std::size_t> > elementTags;
  std::vector<std::vector<std::size_t> > nodeTags;
  gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, 2);
  std::multimap<size_t, size_t> graph;
  std::set<size_t> nodes_to_treat;
  std::set<size_t> nodes;
  size_t numNodesPerElem = (elementTypes[0] == 2) ? 3 : 6;

  // All sub-edges in P2 triangle
  size_t nSubEdges = 15;
  size_t e0[] = {0,0,0,0,0,1,1,1,1,2,2,2,3,3,4};
  size_t e1[] = {1,2,3,4,5,2,3,4,5,3,4,5,4,5,5};

  // for(size_t i = 0; i < elementTags[0].size(); i++) {
  //   for(size_t j = 0; j < numNodesPerElem; ++j) {
  //     size_t n0 = nodeTags[0][numNodesPerElem * i + j];
  //     size_t n1 = nodeTags[0][numNodesPerElem * i + ((j + 1) % numNodesPerElem)];
  //     graph.insert(std::make_pair(n0, n1));
  //     nodes.insert(n0);
  //     nodes.insert(n1);
  //   }
  //   // size_t n0 = nodeTags[0][numNodesPerElem * i + 0];
  //   // size_t n1 = nodeTags[0][numNodesPerElem * i + 1];
  //   // size_t n2 = nodeTags[0][numNodesPerElem * i + 2];
  //   // graph.insert(std::make_pair(n0, n1));
  //   // graph.insert(std::make_pair(n1, n2));
  //   // graph.insert(std::make_pair(n2, n0));
  //   // nodes.insert(n0);
  //   // nodes.insert(n1);
  //   // nodes.insert(n2);
  // }

  for(size_t i = 0; i < elementTags[0].size(); i++) {
    if(numNodesPerElem == 3) {
      // P1 triangles - Add each edge
      for(size_t j = 0; j < 3; ++j) {
        size_t n0 = nodeTags[0][numNodesPerElem * i + j];
        size_t n1 = nodeTags[0][numNodesPerElem * i + (j+1) % numNodesPerElem];
        graph.insert(std::make_pair(n0, n1));
        nodes.insert(n0);
        nodes.insert(n1);
      }
    } else {
      // P2 triangles - Add all possible P2 sub-edges
      for(size_t k = 0; k < nSubEdges; ++k) {
        size_t n0 = nodeTags[0][numNodesPerElem * i + e0[k]];
        size_t n1 = nodeTags[0][numNodesPerElem * i + e1[k]];
        graph.insert(std::make_pair(n0, n1));
        nodes.insert(n0);
        nodes.insert(n1);
      }
    }
  }

  double threshold = tol;

  for(auto n : nodes) {
    double c = C[n];
    double s = S[n];
    for(auto it = graph.lower_bound(n); it != graph.upper_bound(n); ++it) {
      size_t neigh = it->second;
      double cn = C[neigh];
      double sn = S[neigh];
      double pv = fabs(cn * c + sn * s);
      // if(pv > threshold && pv < 1. - threshold) nodes_to_treat.insert(n);
      if(pv < threshold) nodes_to_treat.insert(n);
    }
  }
  int iter = 0;
  while(iter++ < nIter) {
    for(auto n : nodes_to_treat) {
      double c4 = 0;
      double s4 = 0;
      for(auto it = graph.lower_bound(n); it != graph.upper_bound(n); ++it) {
        size_t neigh = it->second;
        double cn = C[neigh];
        double sn = S[neigh];
        double theta = atan2(sn, cn);
        c4 += cos(2. * theta);
        s4 += sin(2. * theta);
      }
      double theta = 0.5 * atan2(s4, c4);
      C[n] = cos(theta);
      S[n] = sin(theta);
    }
  }

  if(plot) {
    FILE *fIso = fopen("isoAfterSmoothing.pos", "w");
    fprintf(fIso, "View\"isoAfterSmoothing\"{\n");
    FILE *fGra = fopen("graAfterSmoothing.pos", "w");
    fprintf(fGra, "View\"graAfterSmoothing\"{\n");

    for(auto n : nodes) {
      double c = C[n];
      double s = S[n];
      if(nodes_to_treat.find(n) != nodes_to_treat.end()) {
        c *= 2;
        s *= 2;
      }
      std::vector<double> coord;
      std::vector<double> par;
      int entityDim, entityTag;
      gmsh::model::mesh::getNode(n, coord, par, entityDim, entityTag);

      fprintf(fGra, "VP(%g,%g,0){%g,%g,0};", coord[0], coord[1], c, s);
      fprintf(fIso, "VP(%g,%g,0){%g,%g,0};", coord[0], coord[1], -s, c);
    }

    fprintf(fIso, "};");
    fclose(fIso);
    fprintf(fGra, "};");
    fclose(fGra);
  }
#else
  UNUSED(C,S,nIter,tol,plot);
#endif
}

void smoothDirections(const std::vector<std::size_t> &nodeTags,
                      const std::vector<double> &coord,
                      std::map<int, MetricTensor> &metrics, int nIter, double tol)
{
  // Diagonalize the metrics and get the principal directions
  // Conserve the same sizes?
  // Smooth the sizes
  // Reassemble the metrics
  // MetricTensor M = metrics.at(1).copy();

  // M.print();

  // double v1[2], v2[2], lambda1, lambda2;
  // M.getEigenvectorsAndEigenvalues(v1, v2, lambda1, lambda2);

  // feInfo("v1 = %+-1.4e - %+-1.4e", v1[0], v1[1]);
  // feInfo("v2 = %+-1.4e - %+-1.4e", v2[0], v2[1]);
  // feInfo("l  = %+-1.4e - %+-1.4e", lambda1, lambda2);

  FILE *directions, *dirSmoothed;
  directions = fopen("directions.pos", "w");
  dirSmoothed = fopen("directionsSmoothed.pos", "w");
  fprintf(directions, "View \" directions \"{\n");
  fprintf(dirSmoothed, "View \" directionsSmoothed \"{\n");

  std::map<size_t, double> COS;
  std::map<size_t, double> SIN;

  size_t N = nodeTags.size();
  double v1[2], v2[2], eigenvalues[2], lambda1, lambda2;

  for(size_t i = 0; i < N; ++i) {

    double x = coord[3 * i + 0];
    double y = coord[3 * i + 1];

    metrics.at(nodeTags[i]).getEigenvectorsAndEigenvalues(v1, v2, lambda1, lambda2);

    double C,S;
    if(lambda1 > lambda2) {
      C = v1[0];
      S = v1[1]; 
    } else {
      C = v2[0];
      S = v2[1]; 
    }
    COS[nodeTags[i]] = C;
    SIN[nodeTags[i]] = S;

    fprintf(directions, "VP(%g,%g,%g){%g,%g,%g};\n", x, y, 0.,  C, S, 0.);
    fprintf(directions, "VP(%g,%g,%g){%g,%g,%g};\n", x, y, 0., -S, C, 0.);
  }

  smoothDirections(COS, SIN, nIter, tol);

  for(size_t i = 0; i < N; ++i) {

    double x = coord[3 * i + 0];
    double y = coord[3 * i + 1];

    metrics.at(nodeTags[i]).getEigenvectorsAndEigenvalues(v1, v2, lambda1, lambda2);

    double C = COS[nodeTags[i]];
    double S = SIN[nodeTags[i]];

    if(lambda1 > lambda2) {
      v1[0] = C;
      v1[1] = S;
      v2[0] = -S;
      v2[1] = C;
    } else {
      v2[0] = C;
      v2[1] = S; 
      v1[0] = -S;
      v1[1] = C;
    }

    eigenvalues[0] = lambda1;
    eigenvalues[1] = lambda2;
    MetricTensor Msmoothed(eigenvalues, v1, v2);
    metrics.at(nodeTags[i]).assignMatrixFrom(Msmoothed);

    fprintf(dirSmoothed, "VP(%g,%g,%g){%g,%g,%g};\n", x, y, 0.,  C, S, 0.);
    fprintf(dirSmoothed, "VP(%g,%g,%g){%g,%g,%g};\n", x, y, 0., -S, C, 0.);
  }

  fprintf(directions, "};"); fclose(directions);
  fprintf(dirSmoothed, "};"); fclose(dirSmoothed);

}

#if defined(HAVE_SOPLEX)
double solveErrorFunction(double k, double *v, double *w, double H[2][2], double C[2][2][2],
  std::set<double> &roots, linearProblem &myLP)
{
  double a1 = k * w[0];
  double a2 = k * w[1];
  // double a1 = 0.0;
  // double a2 = 0.0;
  double v1 = v[0];
  double v2 = v[1];

  double H11 = H[0][0];
  double H12 = H[0][1];
  double H21 = H[1][0];
  double H22 = H[1][1];

  double C111 =  C[0][0][0];
  double C112 = (C[0][0][1] + C[0][1][0] + C[1][0][0]) / 3.;
  double C122 = (C[0][1][1] + C[1][0][1] + C[1][1][0]) / 3.;
  double C222 =  C[1][1][1];

  double cs6 = (C111 * a1 * a1 * a1 + 3. * C112 * a1 * a1 * a2 + 3. * C122 * a1 * a2 * a2 +
                C222 * a2 * a2 * a2) / 120.;
  double cs5 = C111 * a1 * a1 * v1 / 20. + C112 * a1 * a1 * v2 / 20. + C122 * a2 * a2 * v1 / 20. +
               C222 * a2 * a2 * v2 / 20. + C112 * a1 * a2 * v1 / 10. + C122 * a1 * a2 * v2 / 10.;
  double cs4 = C111 * a1 * v1 * v1 / 8. + C112 * a2 * v1 * v1 / 8. + C122 * a1 * v2 * v2 / 8. +
               C222 * a2 * v2 * v2 / 8. + C112 * a1 * v1 * v2 / 4. + C122 * a2 * v1 * v2 / 4.;
  double cs3 = C111 * v1 * v1 * v1 / 6. + C112 * v1 * v1 * v2 / 2. + C122 * v1 * v2 * v2 / 2. +
               C222 * v2 * v2 * v2 / 6.;
  cs4 += (H11 * a1 * a1 + H12 * a1 * a2 + H21 * a1 * a2 + H22 * a2 * a2) / 8.;
  cs3 += (H11 * a1 * v1 + H12 * a2 * v1 + H21 * a1 * v2 + H22 * a2 * v2) / 2.;

  // Ajout du développement limité au premier ordre de Hij = Hij(x0) + s*Cijk*vk
  cs5 += (a1 * a1 / 20.) * (C111 * v1 + C112 * v2);
  cs4 += (v1 * a1 / 8.) * (C111 * v1 + C112 * v2);

  cs5 += (a1 * a2 / 20.) * (C112 * v1 + C122 * v2); // C121*v1 + C122*v2
  cs4 += (v1 * a2 / 8.) * (C112 * v1 + C122 * v2);

  cs5 += (a1 * a2 / 20.) * (C112 * v1 + C122 * v2); // C211*v1 + C212*v2
  cs4 += (v2 * a1 / 8.) * (C112 * v1 + C122 * v2);

  cs5 += (a2 * a2 / 20.) * (C122 * v1 + C222 * v2); // C221*v1 + C222*v2
  cs4 += (v2 * a2 / 8.) * (C122 * v1 + C222 * v2);

  myLP.lvl1PolynomialCoeffs_minus(0) = cs6;
  myLP.lvl1PolynomialCoeffs_minus(1) = cs5;
  myLP.lvl1PolynomialCoeffs_minus(2) = cs4;
  myLP.lvl1PolynomialCoeffs_minus(3) = cs3;
  myLP.lvl1PolynomialCoeffs_minus(4) = 0.;
  myLP.lvl1PolynomialCoeffs_minus(5) = 0.;
  myLP.lvl1PolynomialCoeffs_minus(6) = 1.;

  myLP.lvl1PolynomialCoeffs_plus(0) = cs6;
  myLP.lvl1PolynomialCoeffs_plus(1) = cs5;
  myLP.lvl1PolynomialCoeffs_plus(2) = cs4;
  myLP.lvl1PolynomialCoeffs_plus(3) = cs3;
  myLP.lvl1PolynomialCoeffs_plus(4) = 0.;
  myLP.lvl1PolynomialCoeffs_plus(5) = 0.;
  myLP.lvl1PolynomialCoeffs_plus(6) = -1.;

  // myLP.lvl1PolynomialCoeffs_minus = Eigen::VectorXd::Zero(4);
  // myLP.lvl1PolynomialCoeffs_minus(0) = cs3;
  // myLP.lvl1PolynomialCoeffs_minus(1) = 0.;
  // myLP.lvl1PolynomialCoeffs_minus(2) = 0.;
  // myLP.lvl1PolynomialCoeffs_minus(3) = 1.;

  // myLP.lvl1PolynomialCoeffs_plus = Eigen::VectorXd::Zero(4);
  // myLP.lvl1PolynomialCoeffs_plus(0) = cs3;
  // myLP.lvl1PolynomialCoeffs_plus(1) = 0.;
  // myLP.lvl1PolynomialCoeffs_plus(2) = 0.;
  // myLP.lvl1PolynomialCoeffs_plus(3) = -1.;

  double s = 1e10;
  roots = RootFinder::solvePolynomial(myLP.lvl1PolynomialCoeffs_minus, 0., INFINITY, 1e-8);
  for(auto it = roots.begin(); it != roots.end(); it++) {
    s = fmin(s, *it);
  }

  roots = RootFinder::solvePolynomial(myLP.lvl1PolynomialCoeffs_plus, 0., INFINITY, 1e-8);
  for(auto it = roots.begin(); it != roots.end(); it++) {
    s = fmin(s, *it);
  }

  // feInfo("lowest root is %1.16e", s);
  return s;
}

double evaluateHomogeneousErrorPolynomial(const std::vector<double> &errorCoefficients,
  const int degree, const double xLoc, const double yLoc,
  const bool applyBinomialCoefficients, MetricTensor &Hij)
{
	// Using homogeneous polynomial factorization.
	// See e.g. Coulaud & Loseille, Very High Order Anisotropic Metric-Based Mesh Adaptation in 3D
	// Proposition 2.2. The matrix H has been adapted for 2D polynomial (no factor 2 on non-diagonal entries)
  // double res1 = 0.;
  // for(int i = degree-1, j = 0, cnt = 0; i >= 0; --i, ++j, ++cnt){
  // 	if(i+j == degree-1){
  		
	//   	Hij(0,0) = errorCoefficients[ cnt + 0 + 0 ];
	//   	Hij(0,1) = errorCoefficients[ cnt + 0 + 1 ];
	//   	Hij(1,0) = errorCoefficients[ cnt + 1 + 0 ];
	//   	Hij(1,1) = errorCoefficients[ cnt + 1 + 1 ];
	//   	Hij = Hij.absoluteValueEigen();

	//   	res1 += fabs(pow(xLoc, i) * pow(yLoc, j)) * 
	//   		  (xLoc*xLoc * Hij(0,0) + 2. * xLoc*yLoc * Hij(0,1) + yLoc*yLoc * Hij(1,1));
	//   }
  // }
  // return res1;

  // Classical polynomial evaluation
  double res2 = 0.;
  int expx = errorCoefficients.size() - 1, expy = 0;

  if(applyBinomialCoefficients) {
    // Quick and dirty binomial coefficients
    int coeff[6];
    if(degree == 1) {
      coeff[0] = 1.;
      coeff[1] = 2.;
      coeff[2] = 1.;
    } else if(degree == 2) {
      coeff[0] = 1.;
      coeff[1] = 3.;
      coeff[2] = 3.;
      coeff[3] = 1.;
    } else if(degree == 3) {
      coeff[0] = 1.;
      coeff[1] = 4.;
      coeff[2] = 6.;
      coeff[3] = 4.;
      coeff[4] = 1.;
    } else if(degree == 4) {
      coeff[0] = 1.;
      coeff[1] = 5.;
      coeff[2] = 10.;
      coeff[3] = 10.;
      coeff[4] = 5.;
      coeff[5] = 1.;
    } else {
      feInfo("Implement binomial coefficients for P5+ recoveries");
      exit(-1);
    }

    for(int i = 0; i < errorCoefficients.size(); i++, expx--, expy++) {
      res2 += coeff[i] * errorCoefficients[i] * pow(xLoc, expx) * pow(yLoc, expy);
    }
  } else {
    // Compute the polynomial with the coefficients as-is.
    // This is used when adapting for the H1 seminorm for P3+,
    // where the full coefficients of the gradient of the error
    // polynomial are given.
    for(int i = 0; i < errorCoefficients.size(); i++, expx--, expy++) {
      res2 += errorCoefficients[i] * pow(xLoc, expx) * pow(yLoc, expy);
    }
  }

  return res2;
}

//
// Compute a discretization of the error curve 1 for straight aniso adaptation
//
void computeLvl1(const std::vector<double> &errorCoefficients,
  const int degree, const bool applyBinomialCoefficients, linearProblem &myLP)
{
  double theta, xi, yi, error, root;

  for(size_t i = 0; i < myLP.numConstraints; ++i) {
    theta = 2. * M_PI * (double) i / (double) myLP.numConstraints;
    xi = cos(theta);
    yi = sin(theta);
    error = fmax(1e-14, fabs(evaluateHomogeneousErrorPolynomial(errorCoefficients, degree, xi, yi, applyBinomialCoefficients, myLP.Hij)));
    root = fabs(pow(error, 1. / (degree + 1)));
    myLP.lvl1[2 * i] = xi / root;
    myLP.lvl1[2 * i + 1] = yi / root;
  }

  // Split so that ellipse arc = lengthCurve/N approximately
  if(myLP.uniformErrorCurve) {
    for(int iter = 0; iter < myLP.numLoopsUniformErrorCurve; ++iter) {
      double lengthCurve = 0.;
      for(size_t i = 0; i < myLP.numConstraints; ++i) {
        double x0 = myLP.lvl1[2 * i];
        double y0 = myLP.lvl1[2 * i + 1];
        double x1 = myLP.lvl1[2 * ((i + 1) % myLP.numConstraints)];
        double y1 = myLP.lvl1[2 * ((i + 1) % myLP.numConstraints) + 1];
        double distCurve = sqrt( (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) );
        lengthCurve += distCurve;
      }

      // Scaling factor
      double sumTargetTheta = 0.;
      for(size_t i = 0; i < myLP.numConstraints; ++i) {
        double x0 = myLP.lvl1[2 * i];
        double y0 = myLP.lvl1[2 * i + 1];
        sumTargetTheta += lengthCurve / (myLP.numConstraints * sqrt(x0*x0 + y0*y0));
      }

      theta = 0.;
      for(size_t i = 0; i < myLP.numConstraints; ++i) {
        double x0 = myLP.lvl1[2 * i];
        double y0 = myLP.lvl1[2 * i + 1];
        xi = cos(theta);
        yi = sin(theta);
        double dtheta = lengthCurve / (myLP.numConstraints * sqrt(x0*x0 + y0*y0)) / sumTargetTheta * 2. * M_PI;
        theta += dtheta;
        error = fmax(1e-14, fabs(evaluateHomogeneousErrorPolynomial(errorCoefficients, degree, xi, yi, applyBinomialCoefficients, myLP.Hij)));
        root = fabs(pow(error, 1. / (degree + 1)));
        myLP.lvl1[2 * i] = xi / root;
        myLP.lvl1[2 * i + 1] = yi / root;
      }
    }
  }
}

//
// Compute a discretization of the error curve 1 for curved aniso adaptation
//
void computeLvl1Curved(double K1, double K2, double *g1, double *g2, double Hij[2][2],
                        double Cijk[2][2][2], int nThetaPerQuadrant, linearProblem &myLP)
{
  double v1[2], v2[2], k1, k2, theta, k, v[2], w[2], s, x0[2] = {0., 0.};

  std::set<double> roots;

  for(int iDir = 0; iDir < 4; ++iDir) {
    switch(iDir) {
      case 0:
        v1[0] = g1[0];
        v1[1] = g1[1];
        v2[0] = g2[0];
        v2[1] = g2[1];
        k1 = K1;
        k2 = K2;
        break;
      case 1:
        v1[0] = -g1[0];
        v1[1] = -g1[1];
        v2[0] = -g2[0];
        v2[1] = -g2[1];
        k1 = -K1;
        k2 = -K2;
        break;
      case 2:
        v1[0] = -g2[0];
        v1[1] = -g2[1];
        v2[0] = g1[0];
        v2[1] = g1[1];
        k1 = -K2;
        k2 = K1;
        break;
      case 3:
        v1[0] = g2[0];
        v1[1] = g2[1];
        v2[0] = -g1[0];
        v2[1] = -g1[1];
        k1 = K2;
        k2 = -K1;
        break;
    }

    for(int i = 0; i < nThetaPerQuadrant; ++i) {
      theta = M_PI / 2. * (double)i / (double)nThetaPerQuadrant;
      k = k1 * cos(theta) * cos(theta) + k2 * sin(theta) * sin(theta);
      v[0] = v1[0] * cos(theta) + v2[0] * sin(theta);
      v[1] = v1[1] * cos(theta) + v2[1] * sin(theta);
      w[0] = -v[1];
      w[1] =  v[0];

      s = solveErrorFunction(k, v, w, Hij, Cijk, roots, myLP);

      // Shoot along the parameterization in s:
      myLP.lvl1[iDir * 2 * nThetaPerQuadrant + 2 * i + 0] = x0[0] + s * v[0] + s * s / 2. * k * w[0];
      myLP.lvl1[iDir * 2 * nThetaPerQuadrant + 2 * i + 1] = x0[1] + s * v[1] + s * s / 2. * k * w[1];
    }
  }
}

bool solveLP(linearProblem &myLP, Eigen::Matrix2d &L)
{
  // Removing/adding constraint rows seems to be much faster
  // than changing rows from one solve to the next
  double xi, yi, normSquared, lhs;
  for(int i = 0; i < myLP.numConstraints; ++i) {
    xi = myLP.constraints[2 * i];
    yi = myLP.constraints[2 * i + 1];
    normSquared = xi * xi + yi * yi;
    lhs = -normSquared * log(normSquared);
    myLP.row.add(0, xi * xi);
    myLP.row.add(1, 2. * xi * yi);
    myLP.row.add(2, yi * yi);

    // myLP.lprow.setLhs(lhs);
    // myLP.lprow.setRowVector(myLP.row);
    // myLP.lprowset.add(myLP.lprow);

    // Creating lprow (in setupLinearProblem) may create a segmentation fault?
    // Maybe doing something wrong... Replacing by newLPRow.
    // It does not seem to impact performances too much.
    soplex::LPRow newLPRow(myLP.row, soplex::LPRow::GREATER_EQUAL, lhs);
    myLP.lprowset.add(newLPRow);

    myLP.row.clear();
  }
  myLP.problem.removeRowRangeReal(0, myLP.numConstraints - 1);
  myLP.problem.addRowsReal(myLP.lprowset);
  myLP.lprowset.clear();

  /* Solve LP */
  myLP.stat = myLP.problem.optimize();

  if(myLP.stat == SPxSolver::OPTIMAL || myLP.stat == SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS) {

    if(myLP.stat == SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS) {
      feWarning("LP has beed solved to optimality but unscaled solution contains violations.");
    }
    myLP.problem.getPrimal(myLP.prim);
    L(0, 0) = myLP.prim[0];
    L(0, 1) = myLP.prim[1];
    L(1, 0) = myLP.prim[1];
    L(1, 1) = myLP.prim[2];
    // feInfo("L = %+-1.16e - %+-1.16e - %+-1.16e", L(0,0), L(0,1), L(1,1));
    return true;
  } else {
    switch(myLP.stat) {
      case SPxSolver::ERROR:feInfo("Soplex status: an error occured"); break;
      case SPxSolver::NO_RATIOTESTER: feInfo("No ratiotester loaded"); break;
      case SPxSolver::NO_PRICER: feInfo("No pricer loaded"); break;
      case SPxSolver::NO_SOLVER: feInfo("No linear solver loaded"); break;
      case SPxSolver::NOT_INIT: feInfo("not initialised error"); break;
      case SPxSolver::ABORT_EXDECOMP: feInfo("solve() aborted to exit decomposition simplex"); break;
      case SPxSolver::ABORT_DECOMP: feInfo("solve() aborted due to commence decomposition simplex"); break;
      case SPxSolver::ABORT_CYCLING: feInfo("solve() aborted due to detection of cycling."); break;
      case SPxSolver::ABORT_TIME: feInfo("solve() aborted due to time limit."); break;
      case SPxSolver::ABORT_ITER: feInfo("solve() aborted due to iteration limit."); break;
      case SPxSolver::ABORT_VALUE: feInfo("solve() aborted due to objective limit."); break;
      case SPxSolver::SINGULAR: feInfo("Basis is singular, numerical troubles?"); break;
      case SPxSolver::NO_PROBLEM: feInfo("No Problem has been loaded."); break;
      case SPxSolver::REGULAR: feInfo("LP has a usable Basis (maybe LP is changed)."); break;
      case SPxSolver::RUNNING: feInfo("algorithm is running"); break;
      case SPxSolver::UNKNOWN: feInfo("nothing known on loaded problem."); break;
      case SPxSolver::OPTIMAL: feInfo("LP has been solved to optimality."); break;
      case SPxSolver::UNBOUNDED: feInfo("LP has been proven to be primal unbounded."); break;
      case SPxSolver::INFEASIBLE: feInfo("LP has been proven to be primal infeasible."); break;
      case SPxSolver::INForUNBD: feInfo("LP is primal infeasible or unbounded."); break;
      case SPxSolver::OPTIMAL_UNSCALED_VIOLATIONS: feInfo("LP has beed solved to optimality but unscaled solution contains violations."); break;
    }
    return false;
  }
}

// #define LOG_SIMPLEX_DEBUG

bool computeMetricLogSimplexStraight(const double *x, const std::vector<double> &errorCoefficients,
                                     const int degree, const int nThetaPerQuadrant,
                                     const int maxIter, const double tol,
                                     const bool applyBinomialCoefficients, MetricTensor &Qres,
                                     int &numIter, linearProblem &myLP)
{
  #if defined(LOG_SIMPLEX_DEBUG)
	FILE *ff = fopen("lvl1_straight.pos", "w");
  fprintf(ff, "View\"lvl1\"{\n");
  #endif

  computeLvl1(errorCoefficients, degree, applyBinomialCoefficients, myLP);

  #if defined(LOG_SIMPLEX_DEBUG)
  for(size_t i = 0; i < myLP.numConstraints; ++i) {
    fprintf(ff, "SP(%g,%g,0){%g};\n",
    	x[0] + myLP.lvl1[2 * i],
    	x[1] + myLP.lvl1[2 * i + 1],
    	1.);
  }
  fprintf(ff, "};\n"); fclose(ff);
  #endif

  Eigen::Matrix2d Q, Qprev, Q12, Qm12, L, diff;

  Q = Eigen::Matrix2d::Identity();
  L = Eigen::Matrix2d::Identity();

  double xi, yi;

  #if defined(LOG_SIMPLEX_DEBUG)
  FILE *myfile = fopen("logSimplexDebug_straight.pos", "w");
  fprintf(myfile, "View\"logSimplexDebug\"{\n");
  #endif

  for(int iter = 0; iter < maxIter; ++iter) {
    Qprev = Q;
    Q12 = Q.sqrt();

    // Apply transformation Q12 to the initial constraints
    for(size_t i = 0; i < myLP.numConstraints; ++i) {
      xi = myLP.lvl1[2 * i];
      yi = myLP.lvl1[2 * i + 1];
      myLP.constraints[2 * i] = Q12(0, 0) * xi + Q12(0, 1) * yi;
      myLP.constraints[2 * i + 1] = Q12(1, 0) * xi + Q12(1, 1) * yi;

      // Verification: les nouvelles contraintes doivent être sur la courbe 1 de p o Q^(-1/2)
      // Trivial puisque Q^(-1/2) annule la transformation, et p(xi) = 1.
      // fprintf(myfile, "SP(%g,%g,0){%g,%g,0};",
      // 	x[0] + myLP.constraints[2 * i],
      // 	x[1] + myLP.constraints[2 * i + 1],
      // 	(double) iter, (double) iter);
    }

    // Solve the linear optimization problem for L
    bool success = solveLP(myLP, L);

    if(success) {
    	// Recover Q from L
     	Q = Q12 * L.exp() * Q12;
      diff = Q - Qprev;
      if((diff.norm() / Q.norm()) < tol) {

        // Metric converged to a solution
        numIter = iter;
        Qres(0,0) = Q(0,0);
        Qres(0,1) = Q(0,1);
        Qres(1,0) = Q(1,0);
        Qres(1,1) = Q(1,1);

        #if defined(LOG_SIMPLEX_DEBUG)
        fprintf(myfile, "};\n"); fclose(myfile);
        #endif

        return true;
      }
    } else {

      // LP solver returned an error : return the last valid Q
      feWarning("LP solver returned an error at iteration %d.", iter);
      numIter = iter;
      Qres(0,0) = Q(0,0);
      Qres(0,1) = Q(0,1);
      Qres(1,0) = Q(1,0);
      Qres(1,1) = Q(1,1);

      #if defined(LOG_SIMPLEX_DEBUG)
      drawSingleEllipse(myfile, x, Q, 1, 30);
      fprintf(myfile, "};\n"); fclose(myfile);
      #endif

      return false;
    }
  }

  // Warning: maxIter reached
  feWarning("Max number of iterations (%d) reached when computing log-simplex metric.", maxIter);

  Qres(0,0) = Q(0,0);
  Qres(0,1) = Q(0,1);
  Qres(1,0) = Q(1,0);
  Qres(1,1) = Q(1,1);

  #if defined(LOG_SIMPLEX_DEBUG)
  drawSingleEllipse(myfile, x, Q, 1, 30);
  fprintf(myfile, "};\n"); fclose(myfile);
  #endif

  return false;
}

bool logSimplexOneIteration(const double *x, const std::vector<double> &errorCoefficients,
                            const int degree, const int nThetaPerQuadrant, const bool applyBinomialCoefficients,
                            MetricTensor &Qres, linearProblem &myLP)
{
  #if defined(LOG_SIMPLEX_DEBUG)
  FILE *ff = fopen("lvl1_oneIteration.pos", "w");
  fprintf(ff, "View\"lvl1_oneIteration\"{\n");
  #endif

  computeLvl1(errorCoefficients, degree, applyBinomialCoefficients, myLP);

  #if defined(LOG_SIMPLEX_DEBUG)
  for(size_t i = 0; i < myLP.numConstraints; ++i) {
    fprintf(ff, "SP(%g,%g,0){%g};\n",
      x[0] + myLP.lvl1[2 * i],
      x[1] + myLP.lvl1[2 * i + 1],
      1.);
  }
  fprintf(ff, "};\n"); fclose(ff);
  #endif

  Eigen::Matrix2d L = Eigen::Matrix2d::Identity();

  #if defined(LOG_SIMPLEX_DEBUG)
  FILE *myfile = fopen("logSimplexDebug_oneIteration.pos", "w");
  fprintf(myfile, "View\"logSimplexDebug_oneIteration\"{\n");
  #endif

  // We perform a single iteration, so the constraints are
  // the points on the lvl 1 of the input polynomial
  for(size_t i = 0; i < myLP.numConstraints; ++i) {
    myLP.constraints[2 * i] = myLP.lvl1[2 * i];
    myLP.constraints[2 * i + 1] = myLP.lvl1[2 * i + 1];

    #if defined(LOG_SIMPLEX_DEBUG)
    fprintf(myfile, "SP(%g,%g,0){%g,%g,0};",
     x[0] + myLP.constraints[2 * i],
     x[1] + myLP.constraints[2 * i + 1]);
    #endif
  }

  // Solve the linear optimization problem for L
  bool success = solveLP(myLP, L);

  if(success) {
    // Do not recvoer Q: this is done in feMetric::computeMetricsCurvedReferenceSpace
    // Return the solution in the reference space Q_ref = exp(L)
    Eigen::Matrix2d Q = L.exp();

    Qres(0,0) = Q(0,0);
    Qres(0,1) = Q(0,1);
    Qres(1,0) = Q(1,0);
    Qres(1,1) = Q(1,1);

    #if defined(LOG_SIMPLEX_DEBUG)
    fprintf(myfile, "};\n"); fclose(myfile);
    #endif

    return true;
  } else {

    // LP solver returned an error
    feWarning("Could not solve LP.");

    #if defined(LOG_SIMPLEX_DEBUG)
    drawSingleEllipse(myfile, x, Q, 1, 30);
    fprintf(myfile, "};\n"); fclose(myfile);
    #endif

    return false;
  }
}

// #define LOG_SIMPLEX_DEBUG

bool computeMetricLogSimplexCurved(const int vertex, double *x, double directionGradient[2], feNewRecovery *rec,
                                   MetricTensor &Qres, int maxIter, int nThetaPerQuadrant,
                                   double tol, int &numIter, linearProblem &myLP)
{
  double c = directionGradient[0];
  double s = directionGradient[1];
  double g1[2] = { c, s};
  double g2[2] = {-s, c};

  double fx = rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 0, vertex);
  double fy = rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 1, vertex);

  double fxx =  rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 2, vertex);
  double fxy = (rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 3, vertex) +
                rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 4, vertex)) / 2.;
  double fyy =  rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 5, vertex);

  double fxxx =  rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 6, vertex);
  double fxxy = (rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 7, vertex) + 
                 rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 8, vertex) + 
                 rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 10, vertex)) / 3.;
  double fxyy = (rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 9, vertex) + 
                 rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 11, vertex) + 
                 rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 12, vertex)) / 3.;
  double fyyy =  rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 13, vertex);

  double Hij[2][2] = {{fxx, fxy},
                      {fxy, fyy}};
  double Cijk[2][2][2] = {{{fxxx, fxxy}, {fxxy, fxyy}},
                          {{fxxy, fxyy}, {fxyy, fyyy}}};

  double kappa1 = (- fy * fy * fxx + 2.0 * fx * fy * fxy - fx * fx * fyy) / (pow(fx * fx + fy * fy, 3. / 2.));
  double kappa2 =     (fx * fy * (fyy - fxx) + (fx * fx - fy * fy) * fxy) / (pow(fx * fx + fy * fy, 3. / 2.));

  //
  // Apply the log-simplex algorithm
  //
  #if defined(LOG_SIMPLEX_DEBUG)
  FILE *ff = fopen("lvl1_curved.pos", "w");
  fprintf(ff, "View\"lvl1\"{\n");
  #endif

  // computeLvl1(errorCoefficients, degree, myLP);
  computeLvl1Curved(kappa1, kappa2, g1, g2, Hij, Cijk, nThetaPerQuadrant, myLP);

  #if defined(LOG_SIMPLEX_DEBUG)
  for(size_t i = 0; i < myLP.numConstraints; ++i) {
    fprintf(ff, "SP(%g,%g,0){%g};",
      x[0] + myLP.lvl1[2 * i],
      x[1] + myLP.lvl1[2 * i + 1],
      1.);
  }
  fprintf(ff, "};\n"); fclose(ff);
  #endif

  Eigen::Matrix2d Q, Qprev, Q12, Qm12, L, diff;

  Q = Eigen::Matrix2d::Identity();
  L = Eigen::Matrix2d::Identity();

  double xi, yi;

  #if defined(LOG_SIMPLEX_DEBUG)
  FILE *myfile = fopen("logSimplexDebug_curved.pos", "w");
  fprintf(myfile, "View\"logSimplexDebug\"{\n");
  #endif

  for(int iter = 0; iter < maxIter; ++iter) {
    Qprev = Q;
    Q12 = Q.sqrt();

    // Apply transformation Q12 to the initial constraints
    for(size_t i = 0; i < myLP.numConstraints; ++i) {
      xi = myLP.lvl1[2 * i];
      yi = myLP.lvl1[2 * i + 1];
      myLP.constraints[2 * i] = Q12(0, 0) * xi + Q12(0, 1) * yi;
      myLP.constraints[2 * i + 1] = Q12(1, 0) * xi + Q12(1, 1) * yi;

      // Verification: les nouvelles contraintes doivent être sur la courbe 1 de p o Q^(-1/2)
      // Trivial puisque Q^(-1/2) annule la transformation, et p(xi) = 1.
      // fprintf(myfile, "SP(%g,%g,0){%g,%g,0};",
      //  x[0] + myLP.constraints[2 * i],
      //  x[1] + myLP.constraints[2 * i + 1],
      //  (double) iter, (double) iter);
    }

    // Solve the linear optimization problem for L
    bool success = solveLP(myLP, L);

    if(success) {
      // Recover Q from L
      Q = Q12 * L.exp() * Q12;

      #if defined(LOG_SIMPLEX_DEBUG)
      drawSingleEllipse(myfile, x, Q, 1, 30);
      #endif

      diff = Q - Qprev;
      if((diff.norm() / Q.norm()) < tol) {

        // Metric converged to a solution
        numIter = iter;
        Qres(0,0) = Q(0,0);
        Qres(0,1) = Q(0,1);
        Qres(1,0) = Q(1,0);
        Qres(1,1) = Q(1,1);

        #if defined(LOG_SIMPLEX_DEBUG)
        fprintf(myfile, "};\n"); fclose(myfile);
        #endif

        return true;
      }
    } else {

      // LP solver returned an error : return the last valid Q
      feWarning("LP solver returned an error at iteration %d.", iter);
      numIter = iter;
      Qres(0,0) = Q(0,0);
      Qres(0,1) = Q(0,1);
      Qres(1,0) = Q(1,0);
      Qres(1,1) = Q(1,1);

      #if defined(LOG_SIMPLEX_DEBUG)
      drawSingleEllipse(myfile, x, Q, 1, 30);
      fprintf(myfile, "};\n"); fclose(myfile);
      #endif

      return false;
    }
  }

  // Warning: maxIter reached
  feWarning("Max number of iterations (%d) reached when computing log-simplex metric.", maxIter);

  Qres(0,0) = Q(0,0);
  Qres(0,1) = Q(0,1);
  Qres(1,0) = Q(1,0);
  Qres(1,1) = Q(1,1);

  #if defined(LOG_SIMPLEX_DEBUG)
  drawSingleEllipse(myfile, x, Q, 1, 30);
  fprintf(myfile, "};\n"); fclose(myfile);
  #endif

  return false;
}
#endif

void computeDirectionFieldFromGradient(const double x[2], const int vertex, const double firstDerivatives[2],
  const double tol, double grad[2])
{
  UNUSED(x, vertex);

  // if(_options.useAnalyticDerivatives) {
  //   // Evaluate the analytic gradient
  //   POS[0] = x[0];
  //   POS[1] = x[1];
  //   POS[2] = 0.;
  //   _options.firstDerivatives->eval(0, POS, DU);
  //   gx = DU[0];
  //   gy = DU[1];
  // } else {
  //   gx = recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 0, vertex);
  //   gy = recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 1, vertex);
  // }

  double gx = firstDerivatives[0];
  double gy = firstDerivatives[1];

  double normGrad = sqrt(gx * gx + gy * gy);

  if(normGrad > tol) {
    double theta1 = atan2(gy, gx);
    grad[0] = cos(theta1);
    grad[1] = sin(theta1);
  } else {
    // Gradient norm is too small : directions will be smoothed during gradation
    grad[0] = 1.;
    grad[1] = 0.;
  }
}

inline double maxFabs(const double x, const double threshold)
{
  if(x < 0 && fabs(x) < threshold) return -threshold;
  if(x > 0 && fabs(x) < threshold) return threshold;
  return x;
}

inline double contract(const double G[2], const double u[2])
{
  return G[0]*u[0] + G[1]*u[1];
}

inline double contract(const double H[2][2], const double u[2], const double v[2])
{
  double res = 0.;
  for(int i = 0; i < 2; ++i) {
    for(int j = 0; j < 2; ++j) {
      res += H[i][j] * u[i] * v[j];
    }
  }
  return res;
}

thread_local std::set<double> ROOTS_L_LINEAR;
thread_local Eigen::VectorXd COEFF_L_LINEAR = Eigen::VectorXd::Zero(5);
thread_local Eigen::VectorXd coeffCubic = Eigen::VectorXd::Zero(4);
thread_local Eigen::VectorXd coeffQuadratic = Eigen::VectorXd::Zero(3);

static double smallestRootOfPolynomial_degree4(Eigen::VectorXd &coeff, int &numRoots)
{
  double minRoot = 1e10;
  double tolCoeff = 1e-8;
  double solveTol = 1e-14;

  // feInfo(" ");
  // feInfo("p = [%+-1.6e %+-1.6e %+-1.6e %+-1.6e %+-1.6e]; roots(p)", coeff(0), coeff(1), coeff(2), coeff(3), coeff(4));

  // Remove all small leading coefficients
  // Return a large root for each of those coefficients
  if(fabs(coeff(0)) < tolCoeff) {
    if(fabs(coeff(1)) < tolCoeff) {
      if(fabs(coeff(2)) < tolCoeff) {
        // All coefficients are small.
        // No linear coefficient since polynomial involves L^4,L^3,L^2 only: return large value
        // feInfo("One root: 1e10");
        numRoots = 1;
        return minRoot;
      }

      // Solve for scaled quadratic polynomial
      coeffQuadratic(0) = coeff(2);
      coeffQuadratic(1) = coeff(3);
      coeffQuadratic(2) = coeff(4);
      for(int i = 1; i < coeffQuadratic.size(); ++i) {
        coeffQuadratic(i) /= coeffQuadratic(0);
      }
      coeffQuadratic(0) = 1.;
      // Remove all small coeff multiplying powers of L
      for(int i = 1; i < coeffQuadratic.size() - 1; ++i) {
        if(fabs(coeffQuadratic(i)) < tolCoeff) coeffQuadratic(i) = 0.;
      }
      // Bound the independent term
      coeffQuadratic(2) = maxFabs(coeffQuadratic(2), solveTol);
      ROOTS_L_LINEAR = RootFinder::solvePolynomial(coeffQuadratic, -INFINITY, INFINITY, solveTol);
      numRoots = ROOTS_L_LINEAR.size();
      for(auto it = ROOTS_L_LINEAR.begin(); it != ROOTS_L_LINEAR.end(); it++) {
        minRoot = fmin(minRoot, fabs(*it));
      }
      // feInfo("Two small coeff");
      // feInfo("Solving for quadratic p = [%+-1.6e %+-1.6e %+-1.6e]; roots(p)", coeffQuadratic(0), coeffQuadratic(1), coeffQuadratic(2));
      // for(auto val : ROOTS_L_LINEAR)
      //   feInfo("%+-1.6e", val);
      return minRoot;
    }

    // Solve for scaled cubic polynomial
    coeffCubic(0) = coeff(1);
    coeffCubic(1) = coeff(2);
    coeffCubic(2) = coeff(3);
    coeffCubic(3) = coeff(4);
    for(int i = 1; i < coeffCubic.size(); ++i) {
      coeffCubic(i) /= coeffCubic(0);
    }
    coeffCubic(0) = 1.;
    // Remove all small coeff multiplying powers of L
    for(int i = 1; i < coeffCubic.size() - 1; ++i) {
      if(fabs(coeffCubic(i)) < tolCoeff) coeffCubic(i) = 0.;
    }
    // Bound the independent term
    coeffCubic(3) = maxFabs(coeffCubic(3), solveTol);

    ROOTS_L_LINEAR = RootFinder::solvePolynomial(coeffCubic, -INFINITY, INFINITY, solveTol);
    numRoots = ROOTS_L_LINEAR.size();
    for(auto it = ROOTS_L_LINEAR.begin(); it != ROOTS_L_LINEAR.end(); it++) {
      minRoot = fmin(minRoot, fabs(*it));
    }
    // feInfo("One small coeff");
    // feInfo("Solving for cubic p = [%+-1.6e %+-1.6e %+-1.6e %+-1.6e]; roots(p)", coeffCubic(0), coeffCubic(1), coeffCubic(2), coeffCubic(3));
    // for(auto val : ROOTS_L_LINEAR)
    //   feInfo("%+-1.6e", val);
    return minRoot;
  }

  // Scale polynomial based on leading coefficient
  for(int i = 1; i < coeff.size(); ++i) {
    coeff(i) /= coeff(0);
  }
  coeff(0) = 1.;

  // Remove all small coeff multiplying powers of L
  for(int i = 1; i < coeff.size() - 1; ++i) {
    if(fabs(coeff(i)) < tolCoeff) coeff(i) = 0.;
  }

  // Bound the independent term
  coeff(4) = maxFabs(coeff(4), solveTol);

  ROOTS_L_LINEAR = RootFinder::solvePolynomial(coeff, -INFINITY, INFINITY, solveTol);
  numRoots = ROOTS_L_LINEAR.size();
  for(auto it = ROOTS_L_LINEAR.begin(); it != ROOTS_L_LINEAR.end(); it++) {
    minRoot = fmin(minRoot, fabs(*it));
  }
  // feInfo("No  small coeff");
  // feInfo("Solving for quartic p = [%+-1.6e %+-1.6e %+-1.6e %+-1.6e %+-1.6e]; roots(p)", coeff(0), coeff(1), coeff(2), coeff(3), coeff(4));
  // for(auto val : ROOTS_L_LINEAR)
  //   feInfo("%+-1.6e", val);
  return minRoot;
}

void solveSizePolynomialLinear(const double targetError, const double directionGrad[2], const double du[2], const double d2u[4], double &Liso, double &Lgrad)
{
  // double TOL = 1e-12;

  double fx = du[0]; //maxFabs(du[0], TOL);
  double fy = du[1]; //maxFabs(du[1], TOL);

  double fxx = d2u[0]; //maxFabs(d2u[0], TOL);
  double fxy = d2u[1]; //maxFabs(d2u[1], TOL);
  double fyy = d2u[3]; //maxFabs(d2u[3], TOL);

  double G[2] = {fx,fy};
  double H[2][2] = {{fxx,fxy},{fxy,fyy}};
  // double normGrad = sqrt(fx*fx+fy*fy);

  //
  // Only compute sizes if this is not an isotropic region
  // Acutally: already checked in feMetric.cpp
  //
  // if(normGrad > TOL) {
    // double g[2]     = { fx/normGrad, fy/normGrad};
    // double gOrth[2] = {-fy/normGrad, fx/normGrad};

    // Use the smoothed directions for g and gOrth
    // These directions are unit because computed from cos,sin
    double g[2]     = { directionGrad[0], directionGrad[1]};
    double gOrth[2] = {-directionGrad[1], directionGrad[0]};

    double  kiso = (-fy * fy * fxx + 2.0 * fx * fy * fxy - fx * fx * fyy) / (pow(fx * fx + fy * fy, 3. / 2.));
    double kgrad = (fx * fy * (fyy - fxx) + (fx * fx - fy * fy) * fxy)    / (pow(fx * fx + fy * fy, 3. / 2.));

    // double solveTol = 1e-14;
    double c2, c3, c4;

    ///////////////////////////////////////////////////////
    // Size along isocontour
    ///////////////////////////////////////////////////////
    Liso = 1e10;
    double minRoot;
    int nRootsIso = 0, nTotalRootsIso = 0;

    // Coefficients for integral of error model along quadratic approximation of isoline
    c4 = kiso*kiso/12. * contract(H,g,g);
    c3 = kiso/3. * contract(H,g,gOrth);
    c2 = 1./2. * contract(H,gOrth,gOrth) + kiso/2. * (G[0]*g[0] + G[1]*g[1]);

    // Solve p(L) = error
    COEFF_L_LINEAR(0) = c4;
    COEFF_L_LINEAR(1) = c3;
    COEFF_L_LINEAR(2) = c2;
    COEFF_L_LINEAR(3) = 0.;
    COEFF_L_LINEAR(4) = -targetError;
    minRoot = smallestRootOfPolynomial_degree4(COEFF_L_LINEAR, nRootsIso);
    nTotalRootsIso += nRootsIso;
    Liso = fmin(Liso, minRoot);

    // Solve -p(L) = error
    COEFF_L_LINEAR(0) = -c4;
    COEFF_L_LINEAR(1) = -c3;
    COEFF_L_LINEAR(2) = -c2;
    COEFF_L_LINEAR(3) = 0.;
    COEFF_L_LINEAR(4) = -targetError;
    minRoot = smallestRootOfPolynomial_degree4(COEFF_L_LINEAR, nRootsIso);
    nTotalRootsIso += nRootsIso;
    Liso = fmin(Liso, minRoot);

    if(nTotalRootsIso == 0) {
      feErrorMsg(FE_STATUS_ERROR, "No roots found for iso  polynomials: %+-1.3e L^4 %+-1.3e L^3 %+-1.3e L^2 %+-1.3e = 0 and %+-1.3e L^4 %+-1.3e L^3 %+-1.3e L^2 %+-1.3e = 0",
      c4, c3, c2, -targetError, -c4, -c3, -c2, -targetError);
      feInfo("Derivatives: Grad = [%+-1.4e, %+-1.4e] - Hessian = [%+-1.4e, %+-1.4e, %+-1.4e]", fx, fy, fxx, fxy, fyy);
      feInfo("Curvatures: kiso = %+-1.4e - kgrad = %+-1.4e", kiso, kgrad);
      exit(-1);
    }

    ///////////////////////////////////////////////////////
    // Size along gradient
    ///////////////////////////////////////////////////////
    Lgrad = 1e10;
    int nRootsGrad = 0, nTotalRootsGrad = 0;

    // Coefficients for integral of error model along quadratic approximation of gradient curves
    c4 = kgrad*kgrad/12. * contract(H,gOrth,gOrth);
    c3 = kgrad/3. * contract(H,gOrth,g);
    c2 = 1./2. * contract(H,g,g) + kgrad/2. * (G[0]*gOrth[0] + G[1]*gOrth[1]);

    // Solve p(L) = error
    COEFF_L_LINEAR(0) = c4;
    COEFF_L_LINEAR(1) = c3;
    COEFF_L_LINEAR(2) = c2;
    COEFF_L_LINEAR(3) = 0.;
    COEFF_L_LINEAR(4) = -targetError;
    minRoot = smallestRootOfPolynomial_degree4(COEFF_L_LINEAR, nRootsGrad);
    nTotalRootsGrad += nRootsGrad;
    Lgrad = fmin(Lgrad, minRoot);

    // Solve -p(L) = error
    COEFF_L_LINEAR(0) = -c4;
    COEFF_L_LINEAR(1) = -c3;
    COEFF_L_LINEAR(2) = -c2;
    COEFF_L_LINEAR(3) = 0.;
    COEFF_L_LINEAR(4) = -targetError;
    minRoot = smallestRootOfPolynomial_degree4(COEFF_L_LINEAR, nRootsGrad);
    nTotalRootsGrad += nRootsGrad;
    Lgrad = fmin(Lgrad, minRoot);

    if(nTotalRootsGrad == 0) {
      feErrorMsg(FE_STATUS_ERROR, "No roots found for grad polynomials: %+-1.3e L^4 %+-1.3e L^3 %+-1.3e L^2 %+-1.3e = 0 and %+-1.3e L^4 %+-1.3e L^3 %+-1.3e L^2 %+-1.3e = 0",
      c4, c3, c2, -targetError, -c4, -c3, -c2, -targetError);
      exit(-1);
    }

    // TEST: Use sizes from the original paper (R. Zhang, J.-F. Remacle), without integral Taylor remainder. 
    // Pointwise evaluation of derivatives:
    // double eiso  = fabs(contract(H, gOrth, gOrth) + contract(G,g    ));
    // double egrad = fabs(contract(H, g    , g    ) + contract(G,gOrth));
    // Liso  = sqrt(2. * targetError / eiso );
    // Lgrad = sqrt(2. * targetError / egrad);

    // Compare the size with the one obtained from the error model along the straight directions of the optimal metric.
    double cgrad = 1./2. * contract(H, g, g);
    double ciso  = 1./2. * contract(H, gOrth, gOrth);
    double Lgrad_straight = sqrt( fabs(targetError/cgrad) );
    double Liso_straight  = sqrt( fabs(targetError/ciso ) );

    // feInfo(" ");
    // feInfo("Avec courbure - directions iso : Lgrad = %1.3e - Liso = %1.3e", Lgrad, Liso);
    // feInfo("Sans courbure - directions iso : Lgrad = %1.3e - Liso = %1.3e", Lgrad_straight, Liso_straight);

    Liso = fmax(Liso, Liso_straight);
    Lgrad = fmax(Lgrad, Lgrad_straight);

  // } else {
  //   // Gradient is zero: region is locally isotropic
  //   // Return maximum size


  // }


  if(isnan(Liso) || isnan(Lgrad)) {
    feErrorMsg(FE_STATUS_ERROR, "Size is NaN: hIso = %+-1.3e - hGrad = %+-1.3e", Liso, Lgrad);
    exit(-1);
  }
}

inline double contract(const double C[2][2][2], const double u[2], const double v[2], const double w[2])
{
  double res = 0.;
  for(int i = 0; i < 2; ++i) {
    for(int j = 0; j < 2; ++j) {
      for(int k = 0; k < 2; ++k) {
        res += C[i][j][k] * u[i] * v[j] * w[k];
      }
    }
  }
  return res;
}

// static double analyze2dMetric(const MetricTensor &M, double &C, double &S, double &h1, double &h2)
// {
//   double a = M(0,0);
//   double c = M(0,1);
//   double b = M(1,1);
//   double l1 = 0.5 * (a + b + sqrt((a - b) * (a - b) + 4 * c * c));
//   double l2 = 0.5 * (a + b - sqrt((a - b) * (a - b) + 4 * c * c));
//   h1 = 1. / sqrt(l1);
//   h2 = 1. / sqrt(l2);
//   double theta = 0.5 * atan2(2 * c, a - b);
//   C = cos(theta);
//   S = sin(theta);
//   return theta;
// }

thread_local std::set<double> ROOTS_L_QUADRATIC, ROOTS_TMP;
thread_local Eigen::VectorXd COEFF_L_QUADRATIC = Eigen::VectorXd::Zero(7);
thread_local Eigen::VectorXd coeffQuintic = Eigen::VectorXd::Zero(6);
thread_local Eigen::VectorXd coeffQuartic = Eigen::VectorXd::Zero(5);

static double smallestRootOfPolynomial_degree6(Eigen::VectorXd &coeff, int &numRoots)
{
  double minRoot = 1e10;
  double tolCoeff = 1e-8;
  double solveTol = 1e-14;

  // Remove all small leading coefficients
  // Return a large root for each of those coefficients
  if(fabs(coeff(0)) < tolCoeff) {
    if(fabs(coeff(1)) < tolCoeff) {
      if(fabs(coeff(2)) < tolCoeff) {
        if(fabs(coeff(3)) < tolCoeff) {
          // All coefficients are small.
          // No linear coefficient since polynomial involves L^6, L^5, L^4, L^3 only: return large value
          numRoots = 1;
          return minRoot;
        }

        // Coefficients for L^6, L^5, L^4 are small, but not L^3.
        // Solve for scaled cubic polynomial
        coeffCubic(0) = coeff(3);
        coeffCubic(1) = coeff(4);
        coeffCubic(2) = coeff(5);
        coeffCubic(3) = coeff(6);
        for(int i = 1; i < coeffCubic.size(); ++i) {
          coeffCubic(i) /= coeffCubic(0);
        }
        coeffCubic(0) = 1.;
        // Remove all small coeff multiplying powers of L
        for(int i = 1; i < coeffCubic.size() - 1; ++i) {
          if(fabs(coeffCubic(i)) < tolCoeff) coeffCubic(i) = 0.;
        }
        // Bound the independent term
        coeffCubic(3) = maxFabs(coeffCubic(3), solveTol);
        ROOTS_TMP = RootFinder::solvePolynomial(coeffCubic, -INFINITY, INFINITY, solveTol);
        numRoots = ROOTS_TMP.size();
        for(auto it = ROOTS_TMP.begin(); it != ROOTS_TMP.end(); it++) {
          minRoot = fmin(minRoot, fabs(*it));
        }
        // feInfo("Two small coeff");
        // feInfo("Solving for quadratic p = [%+-1.6e %+-1.6e %+-1.6e]; roots(p)", coeffCubic(0), coeffCubic(1), coeffCubic(2));
        // for(auto val : ROOTS_TMP)
        //   feInfo("%+-1.6e", val);
        return minRoot;
      }

      // Coefficients for L^6, L^5 are small, but not L^4.
      // Solve for scaled quartic polynomial
      coeffQuartic(0) = coeff(2);
      coeffQuartic(1) = coeff(3);
      coeffQuartic(2) = coeff(4);
      coeffQuartic(3) = coeff(5);
      coeffQuartic(4) = coeff(6);
      for(int i = 1; i < coeffQuartic.size(); ++i) {
        coeffQuartic(i) /= coeffQuartic(0);
      }
      coeffQuartic(0) = 1.;
      // Remove all small coeff multiplying powers of L
      for(int i = 1; i < coeffQuartic.size() - 1; ++i) {
        if(fabs(coeffQuartic(i)) < tolCoeff) coeffQuartic(i) = 0.;
      }
      // Bound the independent term
      coeffQuartic(4) = maxFabs(coeffQuartic(4), solveTol);

      ROOTS_TMP = RootFinder::solvePolynomial(coeffQuartic, -INFINITY, INFINITY, solveTol);
      numRoots = ROOTS_TMP.size();
      for(auto it = ROOTS_TMP.begin(); it != ROOTS_TMP.end(); it++) {
        minRoot = fmin(minRoot, fabs(*it));
      }
      return minRoot;
    }

    // Coefficient for L^6 is small, but not L^5.
    // Solve for scaled quintic polynomial
    coeffQuintic(0) = coeff(1);
    coeffQuintic(1) = coeff(2);
    coeffQuintic(2) = coeff(3);
    coeffQuintic(3) = coeff(4);
    coeffQuintic(4) = coeff(5);
    coeffQuintic(5) = coeff(6);
    for(int i = 1; i < coeffQuintic.size(); ++i) {
      coeffQuintic(i) /= coeffQuintic(0);
    }
    coeffQuintic(0) = 1.;
    // Remove all small coeff multiplying powers of L
    for(int i = 1; i < coeffQuintic.size() - 1; ++i) {
      if(fabs(coeffQuintic(i)) < tolCoeff) coeffQuintic(i) = 0.;
    }
    // Bound the independent term
    coeffQuintic(5) = maxFabs(coeffQuintic(5), solveTol);

    ROOTS_TMP = RootFinder::solvePolynomial(coeffQuintic, -INFINITY, INFINITY, solveTol);
    numRoots = ROOTS_TMP.size();
    for(auto it = ROOTS_TMP.begin(); it != ROOTS_TMP.end(); it++) {
      minRoot = fmin(minRoot, fabs(*it));
    }
    return minRoot;
  }

  // All coefficients are not negligible
  // Scale polynomial based on leading coefficient
  for(int i = 1; i < coeff.size(); ++i) {
    coeff(i) /= coeff(0);
  }
  coeff(0) = 1.;

  // Remove all small coeff multiplying powers of L
  for(int i = 1; i < coeff.size() - 1; ++i) {
    if(fabs(coeff(i)) < tolCoeff) coeff(i) = 0.;
  }

  // Bound the independent term
  coeff(6) = maxFabs(coeff(6), solveTol);

  ROOTS_TMP = RootFinder::solvePolynomial(coeff, -INFINITY, INFINITY, solveTol);
  numRoots = ROOTS_TMP.size();
  for(auto it = ROOTS_TMP.begin(); it != ROOTS_TMP.end(); it++) {
    minRoot = fmin(minRoot, fabs(*it));
  }
  // feInfo("No  small coeff");
  // feInfo("Solving for quartic p = [%+-1.6e %+-1.6e %+-1.6e %+-1.6e %+-1.6e]; roots(p)", coeff(0), coeff(1), coeff(2), coeff(3), coeff(4));
  // for(auto val : ROOTS_TMP)
  //   feInfo("%+-1.6e", val);
  return minRoot;
}

void solveSizePolynomialQuadratic(const double targetError, const double directionGrad[2], const double du[2], 
  const double d2u[4], const double d3u[8], double &Liso, double &Lgrad)
{
  // double TOL = 1e-12;

  double fx = du[0]; //maxFabs(du[0], TOL);
  double fy = du[1]; //maxFabs(du[1], TOL);

  double fxx = d2u[0]; // maxFabs(d2u[0], TOL);
  double fxy = d2u[1]; // maxFabs(d2u[1], TOL);
  double fyy = d2u[3]; // maxFabs(d2u[3], TOL);

  double fxxx = d3u[0]; //maxFabs(d3u[0], TOL);
  double fxxy = d3u[1]; //maxFabs(d3u[1], TOL);
  double fxyy = d3u[3]; //maxFabs(d3u[3], TOL);
  double fyyy = d3u[7]; //maxFabs(d3u[7], TOL);

  double H[2][2] = {{fxx,fxy},{fxy,fyy}};
  double C[2][2][2] = {{{fxxx,fxxy},{fxxy,fxyy}},{{fxxy,fxyy},{fxyy,fyyy}}};

  // double normGrad = fmax(TOL, sqrt(fx*fx+fy*fy));

  // double g[2]     = { fx/normGrad, fy/normGrad};
  // double gOrth[2] = {-fy/normGrad, fx/normGrad};

   // Use the smoothed directions for g and gOrth
  // These directions are unit because computed from cos,sin
  double g[2]     = { directionGrad[0], directionGrad[1]};
  double gOrth[2] = {-directionGrad[1], directionGrad[0]};

  double  kiso = (-fy * fy * fxx + 2.0 * fx * fy * fxy - fx * fx * fyy) / (pow(fx * fx + fy * fy, 3. / 2.));
  double kgrad = (fx * fy * (fyy - fxx) + (fx * fx - fy * fy) * fxy)    / (pow(fx * fx + fy * fy, 3. / 2.));

  double c3, c4, c5, c6;

  // Solve for the size along isoline (TODO: go in both directions? forward and backward along curve)
  // Coefficients for integral of error model along quadratic approximation of isoline
  // c6 = kiso*kiso*kiso/120. * contract(C,g,g,g);
  // c5 = kiso*kiso/20. * contract(C,g,g,gOrth);
  // c4 = kiso/8. * contract(C,g,gOrth,gOrth) + kiso*kiso/8. * contract(H,g,g);
  // c3 = 1./6. * contract(C,gOrth,gOrth,gOrth) + kiso/2. * contract(H,g,gOrth);

  // Coefficients when the Hessian is expanded as H = H0 + s * Cijk * gammaDot^k,
  // since the third derivatives are known
  c6 = kiso*kiso*kiso/30. * contract(C,g,g,g);
  c5 = 3*kiso*kiso/20. * contract(C,g,g,gOrth);
  c4 = kiso/4. * contract(C,g,gOrth,gOrth)   + kiso*kiso/8. * contract(H,g,g);
  c3 = 1./6. * contract(C,gOrth,gOrth,gOrth) + kiso/2. * contract(H,g,gOrth);

  Liso = 1e10;
  // int numRoots = 0;
  double minRoot;
  int nRootsIso = 0, nTotalRootsIso = 0;

  // Solve p(L) = error
  COEFF_L_QUADRATIC(0) = c6;
  COEFF_L_QUADRATIC(1) = c5;
  COEFF_L_QUADRATIC(2) = c4;
  COEFF_L_QUADRATIC(3) = c3;
  COEFF_L_QUADRATIC(4) = 0.;
  COEFF_L_QUADRATIC(5) = 0.;
  COEFF_L_QUADRATIC(6) = -targetError;
  // ROOTS_L_QUADRATIC = RootFinder::solvePolynomial(COEFF_L_QUADRATIC, -INFINITY, INFINITY, 1e-14);
  // numRoots += ROOTS_L_QUADRATIC.size();
  // for(auto it = ROOTS_L_QUADRATIC.begin(); it != ROOTS_L_QUADRATIC.end(); it++) {
  //   Liso = fmin(Liso, fabs(*it));
  // }
  // nTotalRootsIso = numRoots;
  minRoot = smallestRootOfPolynomial_degree6(COEFF_L_QUADRATIC, nRootsIso);
  nTotalRootsIso += nRootsIso;
  Liso = fmin(Liso, minRoot);

  // Solve -p(L) = error
  COEFF_L_QUADRATIC(0) = -c6;
  COEFF_L_QUADRATIC(1) = -c5;
  COEFF_L_QUADRATIC(2) = -c4;
  COEFF_L_QUADRATIC(3) = -c3;
  COEFF_L_QUADRATIC(4) = 0.;
  COEFF_L_QUADRATIC(5) = 0.;
  COEFF_L_QUADRATIC(6) = -targetError;
  // ROOTS_L_QUADRATIC = RootFinder::solvePolynomial(COEFF_L_QUADRATIC, -INFINITY, INFINITY, 1e-14);
  // numRoots += ROOTS_L_QUADRATIC.size();
  // for(auto it = ROOTS_L_QUADRATIC.begin(); it != ROOTS_L_QUADRATIC.end(); it++) {
  //   Liso = fmin(Liso, fabs(*it));
  // }
  // nTotalRootsIso = numRoots;
  minRoot = smallestRootOfPolynomial_degree6(COEFF_L_QUADRATIC, nRootsIso);
  nTotalRootsIso += nRootsIso;
  Liso = fmin(Liso, minRoot);

  if(nTotalRootsIso == 0) {
    feErrorMsg(FE_STATUS_ERROR, "No roots found for P2 iso polynomials: %+-1.3e L^6 %+-1.3e L^5 %+-1.3e L^4 %+-1.3e L^3 +/- %1.3e = 0",
    c6, c5, c4, c3, targetError);
    exit(-1);
  }

  // Solve for size along gradient curve
  // Coefficients for integral of error model along quadratic approximation of gradient curves
  // c6 = kgrad*kgrad*kgrad/120. * contract(C,gOrth,gOrth,gOrth);
  // c5 = kgrad*kgrad/20. * contract(C,gOrth,gOrth,g);
  // c4 = kgrad/8. * contract(C,gOrth,g,g) + kgrad*kgrad/8. * contract(H,gOrth,gOrth);
  // c3 = 1./6. * contract(C,g,g,g) + kgrad/2. * contract(H,gOrth,g);

  // Coefficients when the Hessian is expanded as H = H0 + s * Cijk * gammaDot^k,
  // since the third derivatives are known
  c6 = kgrad*kgrad*kgrad/30. * contract(C,gOrth,gOrth,gOrth);
  c5 = 3*kgrad*kgrad/20. * contract(C,gOrth,gOrth,g);
  c4 = kgrad/4. * contract(C,gOrth,g,g) + kgrad*kgrad/8. * contract(H,gOrth,gOrth);
  c3 = 1./6. * contract(C,g,g,g) + kgrad/2. * contract(H,gOrth,g);

  Lgrad = 1e10;
  // numRoots = 0;
  int nRootsGrad = 0, nTotalRootsGrad = 0;

  // Solve p(L) = error
  COEFF_L_QUADRATIC(0) = c6;
  COEFF_L_QUADRATIC(1) = c5;
  COEFF_L_QUADRATIC(2) = c4;
  COEFF_L_QUADRATIC(3) = c3;
  COEFF_L_QUADRATIC(4) = 0.;
  COEFF_L_QUADRATIC(5) = 0.;
  COEFF_L_QUADRATIC(6) = -targetError;
  // ROOTS_L_QUADRATIC = RootFinder::solvePolynomial(COEFF_L_QUADRATIC, -INFINITY, INFINITY, 1e-14);
  // numRoots += ROOTS_L_QUADRATIC.size();
  // for(auto it = ROOTS_L_QUADRATIC.begin(); it != ROOTS_L_QUADRATIC.end(); it++) {
  //   Lgrad = fmin(Lgrad, fabs(*it));
  // }
  // nTotalRootsGrad = numRoots;
  minRoot = smallestRootOfPolynomial_degree6(COEFF_L_QUADRATIC, nRootsGrad);
  nTotalRootsGrad += nRootsGrad;
  Lgrad = fmin(Lgrad, minRoot);

  // Solve -p(L) = error
  COEFF_L_QUADRATIC(0) = -c6;
  COEFF_L_QUADRATIC(1) = -c5;
  COEFF_L_QUADRATIC(2) = -c4;
  COEFF_L_QUADRATIC(3) = -c3;
  COEFF_L_QUADRATIC(4) = 0.;
  COEFF_L_QUADRATIC(5) = 0.;
  COEFF_L_QUADRATIC(6) = -targetError;
  // ROOTS_L_QUADRATIC = RootFinder::solvePolynomial(COEFF_L_QUADRATIC, -INFINITY, INFINITY, 1e-14);
  // numRoots += ROOTS_L_QUADRATIC.size();
  // for(auto it = ROOTS_L_QUADRATIC.begin(); it != ROOTS_L_QUADRATIC.end(); it++) {
  //   Lgrad = fmin(Lgrad, fabs(*it));
  // }
  // nTotalRootsGrad = numRoots;
  minRoot = smallestRootOfPolynomial_degree6(COEFF_L_QUADRATIC, nRootsGrad);
  nTotalRootsGrad += nRootsGrad;
  Lgrad = fmin(Lgrad, minRoot);

  if(nTotalRootsGrad == 0) {
    feErrorMsg(FE_STATUS_ERROR, "No roots found for P2 grad polynomials: %+-1.3e L^6 %+-1.3e L^5 %+-1.3e L^4 %+-1.3e L^3 +/- %1.3e = 0",
    c6, c5, c4, c3, targetError);
    exit(-1);
  }

  // Compare the size with the one obtained from the error model along the straight directions
  double cgrad = 1./6. * contract(C, g, g, g);
  double ciso  = 1./6. * contract(C, gOrth, gOrth, gOrth);
  double Lgrad_straight = pow( fabs(targetError/cgrad), 1./3. );
  double Liso_straight  = pow( fabs(targetError/ciso ), 1./3. );

  // feInfo(" ");
  // feInfo("Avec courbure - directions iso : Lgrad = %1.3e - Liso = %1.3e", Lgrad, Liso);
  // feInfo("Sans courbure - directions iso : Lgrad = %1.3e - Liso = %1.3e", Lgrad_straight, Liso_straight);

  Liso = fmax(Liso, Liso_straight);
  Lgrad = fmax(Lgrad, Lgrad_straight);
}

// double secondDerivativesAlongCurve(const double x[2], const int vertex, const double directionV1[2], const int direction, feNewRecovery *recoveredField)
// {
//   double fx = recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 0, vertex);
//   double fy = recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 1, vertex);

//   double fxx = recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 2, vertex);
//   double fxy = (recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 3, vertex) +
//                 recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 4, vertex)) /
//                2.;
//   double fyy = recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 5, vertex);

//   double TOL = 1e-12;

//   if(_options.useAnalyticDerivatives) {
//     // Evaluate the analytic high-order derivatives instead
//     POS[0] = x[0];
//     POS[1] = x[1];
//     POS[2] = 0.;
//     _options.firstDerivatives->eval(0, POS, DU);
//     fx = maxFabs(DU[0], TOL);
//     fy = maxFabs(DU[1], TOL);
//     _options.secondDerivatives->eval(0, POS, D2U);
//     fxx = maxFabs(D2U[0], TOL);
//     fxy = maxFabs(D2U[1], TOL);
//     fyy = maxFabs(D2U[3], TOL);
//   }

//   double kappa1 = (-fy * fy * fxx + 2.0 * fx * fy * fxy - fx * fx * fyy) / (pow(fx * fx + fy * fy, 3. / 2.));
//   double kappa2 = (fx * fy * (fyy - fxx) + (fx * fx - fy * fy) * fxy) / (pow(fx * fx + fy * fy, 3. / 2.));
//   double kappa[2] = {kappa1, kappa2};

//   double C = directionV1[0];
//   double S = directionV1[1];
//   double g[2], gOrth[2];

//   if(direction == 0) {
//     g[0] = C;
//     g[1] = S;
//     gOrth[0] = -S;
//     gOrth[1] = C;
//   } else {
//     gOrth[0] = C;
//     gOrth[1] = S;
//     g[0] = -S;
//     g[1] = C;
//   }

//   double d2 = g[0] * g[0] * fxx + 2. * g[0] * g[1] * fxy + g[1] * g[1] * fyy;
//   // double d1 = 2. * kappa[direction] * (gOrth[0] * fx + gOrth[1] * fy);
//   double d1 = kappa[direction] * (g[0] * fx + g[1] * fy);
//   feInfo("dir = %d - k = %f - 1/r = %f - %f - %f - %1.4e", direction, kappa[direction], 1./sqrt(x[0]*x[0] + x[1]*x[1]), d2, d1, d2+d1);

//   return d2 + d1;
// }

// double thirdDerivativesAlongCurve(const double x[2], const int vertex, const double directionV1[2], const int direction, feNewRecovery *recoveredField)
// {
//   double fx = recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 0, vertex);
//   double fy = recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 1, vertex);

//   double fxx = recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 2, vertex);
//   double fxy = (recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 3, vertex) +
//                 recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 4, vertex)) /
//                2.;
//   double fyy = recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 5, vertex);

//   double fxxx = recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 6, vertex);
//   double fxxy = (recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 7, vertex) +
//                  recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 8, vertex) +
//                  recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 10, vertex)) /
//                 3.;
//   double fxyy = (recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 9, vertex) +
//                  recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 11, vertex) +
//                  recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 12, vertex)) /
//                 3.;
//   double fyyy = recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 13, vertex);

//   double TOL = 1e-12;

//   if(_options.useAnalyticDerivatives) {
//     // Evaluate the analytic high-order derivatives instead
//     POS[0] = x[0];
//     POS[1] = x[1];
//     POS[2] = 0.;
//     _options.firstDerivatives->eval(0, POS, DU);
//     fx = maxFabs(DU[0], TOL);
//     fy = maxFabs(DU[1], TOL);
//     _options.secondDerivatives->eval(0, POS, D2U);
//     fxx = maxFabs(D2U[0], TOL);
//     fxy = maxFabs(D2U[1], TOL);
//     fyy = maxFabs(D2U[3], TOL);
//     _options.thirdDerivatives->eval(0, POS, D3U);
//     fxxx = maxFabs(D3U[0], TOL);
//     fxxy = maxFabs(D3U[1], TOL);
//     fxyy = maxFabs(D3U[3], TOL);
//     fyyy = maxFabs(D3U[7], TOL);
//   }

//   double kappa1 = (-fy * fy * fxx + 2.0 * fx * fy * fxy - fx * fx * fyy) / (pow(fx * fx + fy * fy, 3. / 2.));
//   double kappa2 = (fx * fy * (fyy - fxx) + (fx * fx - fy * fy) * fxy) / (pow(fx * fx + fy * fy, 3. / 2.));
//   double kappa[2] = {kappa1, kappa2};

//   double C = directionV1[0];
//   double S = directionV1[1];
//   double g[2], gOrth[2];

//   if(direction == 0) {
//     g[0] = C;
//     g[1] = S;
//     gOrth[0] = -S;
//     gOrth[1] = C;
//   } else {
//     gOrth[0] = C;
//     gOrth[1] = S;
//     g[0] = -S;
//     g[1] = C;
//   }

//   double d3 = g[0] * g[0] * g[0] * fxxx + 3. * g[0] * g[0] * g[1] * fxxy +
//               3. * g[0] * g[1] * g[1] * fxyy + g[1] * g[1] * g[1] * fyyy;

//   double d2 = 3. * kappa[direction] * (g[0] * gOrth[0] * fxx + (g[0] * gOrth[1] + g[1] * gOrth[0]) * fxy + g[1] * gOrth[1] * fyy);

//   return d3 + d2;
// }


//
// Reconstruct the underlying 2D manifold of the given Riemannian metric,
// assuming the metric is induced by the surface.
// At each point, reconstruct the tangent space from M, then choose the next point
// on the surface by moving on the tangent space by a distance determined by
// the element size, in the parameter space.
//
void reconstructManifoldFromMetric(const std::vector<std::size_t> &nodeTags,
  const std::vector<double> &coord, const std::map<int, MetricTensor> &metricsAtNodeTags,
  const std::vector<int> sequentialTag2nodeTag,
  feNewRecovery *recoveredField)
{
  FILE *myFile = fopen("reconstructedManifold.pos", "w");
  fprintf(myFile, "View \" reconstructedManifold \"{\n");
  FILE *myFile2 = fopen("treated.pos", "w");
  fprintf(myFile2, "View \" treated \"{\n");

  size_t numVertices = nodeTags.size();
  
  // An arbitrary first node on the surface, since only its derivatives are 
  // known by the induced metric
  double z0 = 0.;
  double nextz0 = 0.;

  double dx = 0.01;
  double dy = 0.01;

  // Find a scaling factor so that c*M is an induced metric,
  // that is, there is no diagonal entry such that sqrt(Mii^2-1) is nan
  double c = 1e10;
  for(size_t i = 0; i < numVertices; ++i) {
    const MetricTensor &M = metricsAtNodeTags.at(nodeTags[i]);
    c = fmin(c, M(0,0));
    c = fmin(c, M(1,1));
  }

  c = c * 0.99;

  std::map<int, std::set<int>> &verticesPatches = recoveredField->getVerticesPatches_singleLayer();


  std::vector<double> allZ0(numVertices);

  std::vector<size_t> treated;
  size_t numTreated = 0;
  size_t i = 0;
  treated.push_back(i);

  while(numTreated < numVertices) {

    // allZ0[i] = z0;

    size_t tag = sequentialTag2nodeTag[i];
    const double x = coord[3 * (tag-1) + 0];
    const double y = coord[3 * (tag-1) + 1];

    fprintf(myFile2, "SP(%f,%f,0.){1};\n", x, y);

    // const MetricTensor &M = metricsAtNodeTags.at(nodeTags[i]);
    MetricTensor M(1.);
    c = 1.;
    double FX = 2.*x;
    M(0,0) = 1. + FX*FX;

    // Determine fx and fyL solve nonlinear least-squares with Newton method
    int iter = 0;
    double normRes, residual[2];
    double A = sqrt(M(0,0)/c - 1.);
    double B =      M(0,1)/c;
    double C = sqrt(M(1,1)/c - 1.);

    if(isnan(A) || isnan(C)) {
      feInfo("Nan: %f - %f - %f", M(0,0)/c, M(0,1)/c, M(1,1)/c);
      exit(-1);
    }
    double fx = sqrt(M(0,0)/c - 1.);
    double fy = sqrt(M(1,1)/c - 1.);

    do {
      double a = fx;
      double b = fy;
      fx += (A-a+A*(a*a)+(B*B)*a+a*(b*b)-a*a*a+(a*a*a)*(b*b)+B*C-C*a*b*2.0-B*(a*a)*b*2.0)/(-B*B+a*a+b*b-(a*a)*(b*b)*3.0+B*a*b*4.0+1.0);
      fy += (C-b+(B*B)*b+C*(b*b)+(a*a)*b-b*b*b+(a*a)*(b*b*b)+A*B-A*a*b*2.0-B*a*(b*b)*2.0)/(-B*B+a*a+b*b-(a*a)*(b*b)*3.0+B*a*b*4.0+1.0);
      residual[0] = fx - A - fy * (B - fx*fy);
      residual[1] = fy - C - fx * (B - fx*fy);
      normRes = sqrt(residual[0]*residual[0] + residual[1]*residual[1]);
      feInfo("iter %2d - fx = %+-1.3e - fy = %+-1.3e - res = (%+-1.3e,%+-1.3e) - norm = %+-1.3e", iter, fx, fy, residual[0], residual[1], normRes);
      iter++;
    } while(iter < 20 && normRes > 1e-2);

    if(iter < 20) {
      feInfo("");
      feInfo("Converged: %+-1.3e - %+-1.3e  vs %+-1.3e - %+-1.3e", fx, fy, sqrt(M(0,0)/c - 1.), sqrt(M(1,1)/c - 1.));
      feInfo("");
    }

    // Diagonal term is unused: 3 equations for 2 variables
    // double fx = sqrt(M(0,0)/c - 1.);
    // double fy = sqrt(M(1,1)/c - 1.);
    // double fx = 1.;
    // double fy = 1.;
    // double fx = 2.*x;
    // double fy = 2.*y;

    fprintf(myFile, "SP(%f,%f,%f){1};\n", x, y, z0);

    int next = -1;
    for(auto v : verticesPatches.at(i)) {

      const double xj = coord[3 * (sequentialTag2nodeTag[v]-1) + 0];
      const double yj = coord[3 * (sequentialTag2nodeTag[v]-1) + 1];
      dx = xj - x;
      dy = yj - y;

      // Plot the tangent space at p0
      double SCALING = 1.;
      fx /= SCALING;
      fy /= SCALING;
      fprintf(myFile, "SP(%f,%f,%f){1};\n", xj, yj, z0);
      fprintf(myFile, "SQ(%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f){1,1,1,1};\n",
        x - dx, y - dy, (z0 - fx*dx - fy*dy),
        x + dx, y - dy, (z0 + fx*dx - fy*dy),
        x + dx, y + dy, (z0 + fx*dx + fy*dy),
        x - dx, y + dy, (z0 - fx*dx + fy*dy));

      // Set next vertex if it was not treated yet
      // Set the next z0
      if(next == -1 && std::find(treated.begin(), treated.end(), v) == treated.end()) {
        treated.push_back(v);
        next = v;
        nextz0 = z0 + fx * dx + fy * dy;
        allZ0[v] = nextz0;
      }
    }

    // All neighbours were already treated
    // Restart from any treated vertex that has untreated neighbours
    if(next == -1) {

      bool OK = false;
      for(size_t j = 0; j < numVertices; ++j) {
        // Check if j was treated
        if(std::find(treated.begin(), treated.end(), j) != treated.end()) {
          for(auto v : verticesPatches.at(j)) {
            if(std::find(treated.begin(), treated.end(), v) == treated.end()) {
              OK = true; // At least one neighbour was not treated
              next = j;
              nextz0 = allZ0[j];
              goto label;
            }
          }
        }
      }

    label:
      if(!OK) {
        feInfo("Erreur: aucun voisin n'a de voisin libre");
        fprintf(myFile, "};"); fclose(myFile);
        fprintf(myFile2, "};"); fclose(myFile2);
        exit(-1);
      }
    }

    i = next;
    z0 = nextz0;
    numTreated = treated.size();
  }


  // size_t sizeTreated = treated.size();
  // size_t i = 0;
  // while(sizeTreated < numVertices) {

  //   const double x = coord[3 * i + 0];
  //   const double y = coord[3 * i + 1];

  //   const MetricTensor &M = metricsAtNodeTags.at(nodeTags[i]);
  //   // Diagonal term is unused: 3 equations for 2 variables
  //   double fx = sqrt(M(0,0)/c - 1.);
  //   double fy = sqrt(M(1,1)/c - 1.);
  //   // double fx = 2.*x;
  //   // double fy = 2.*y;
  //   // feInfo("Check: %f", fx*fy - M(0,1));

  //   // Update z0 and dx,dy by moving towards the next (CLOSEST! otherwise Taylor expansion is wrong) vertex on the mesh (parameter space)
  //   // Determine next untreated closest vertex
  //   double distMin = 1e10;
  //   size_t bestj;
  //   for(size_t j = 0; j < numVertices; ++j) {
  //     const double xj = coord[3 * j + 0];
  //     const double yj = coord[3 * j + 1];
  //     double dist = sqrt((x-xj)*(x-xj) + (y-yj)*(y-yj));
  //     if(dist < distMin && j != i && std::find(treated.begin(), treated.end(), j) == treated.end()) {
  //       distMin = dist;
  //       bestj = j;
  //     }
  //   }
  //   treated.push_back(bestj);
  //   const double xNext = coord[3 * bestj + 0];
  //   const double yNext = coord[3 * bestj + 1];
  //   dx = xNext - x;
  //   dy = yNext - y;
  //   feInfo("z0 = %f", z0);

  //   // Plot the tangent space at p0
  //   double SCALING = 100.;
  //   fx /= SCALING;
  //   fy /= SCALING;
  //   fprintf(myFile, "SP(%f,%f,%f){1};\n", x, y, z0);
  //   fprintf(myFile, "SQ(%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f){1,1,1,1};\n",
  //     x - dx, y - dy, (z0 - fx*dx - fy*dy),
  //     x + dx, y - dy, (z0 + fx*dx - fy*dy),
  //     x + dx, y + dy, (z0 + fx*dx + fy*dy),
  //     x - dx, y + dy, (z0 - fx*dx + fy*dy));

  //   z0 += fx * dx + fy * dy;
  //   i = bestj;
  //   sizeTreated++;
  // } 

  fprintf(myFile, "};"); fclose(myFile);
  fprintf(myFile2, "};"); fclose(myFile2);
}

inline double dotprod(const double v1[3], const double v2[3])
{
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

void computeGraphSurfaceMetric(const double *x, const double *df, const double *d2f, const double *d3f, MetricTensor &M, FILE *file)
{
  UNUSED(d3f);
  
  double fu = df[0];
  double fv = df[1];

  double fuu = d2f[0];
  double fuv = d2f[1];
  double fvv = d2f[3];

  // double fuuu = d3f[0];
  // double fuuv = d3f[1];
  // double fuvv = d3f[3];
  // double fvvv = d3f[7];

  double xu[3] = {1.,0.,fu};
  double xv[3] = {0.,1.,fv};

  double xuu[3] = {0.,0.,fuu};
  double xuv[3] = {0.,0.,fuv};
  double xvv[3] = {0.,0.,fvv};

  // double xuuu[3] = {0.,0.,fuuu};
  // double xuuv[3] = {0.,0.,fuuv};
  // double xuvv[3] = {0.,0.,fuvv};
  // double xvvv[3] = {0.,0.,fvvv};

  // Unit normal
  double nN = sqrt(1. + fu*fu + fv*fv);
  double N[3] = {fu/nN, fv/nN, -1./nN};

  // Compute the principal curvatures and directions at the given point
  // double I[2][2] = {{1. + fu*fu, fu*fv},
  //                   {fu*fv, 1. + fv*fv}};
  // double II[2][2] = {{dotprod(xuu,N), dotprod(xuv,N)},
  //                    {dotprod(xuv,N), dotprod(xvv,N)}};
  // double invI_times_II[2][2] =  {{-(fuu*fv*fv - fu*fuv*fv + fuu)/(fu*fu + fv*fv + 1), -(fuv*fv*fv - fu*fvv*fv + fuv)/(fu*fu + fv*fv + 1)}
  //                                {-(fuv*fu*fu - fuu*fv*fu + fuv)/(fu*fu + fv*fv + 1), -(fvv*fu*fu - fuv*fv*fu + fvv)/(fu*fu + fv*fv + 1)}};

  // double G = (fuu*fvv - fuv*fuv)/((1 + fu*fu+fv*fv)*(1 + fu*fu+fv*fv));
  // double H = ((1+fv*fv)*fuu - 2*fu*fv*fuv + (1+fu*fu)*fvv)/(2*pow(1+fu*fu+fv*fv, 3./2.));
  // double kmax = H + sqrt(H*H-G);
  // double kmin = H - sqrt(H*H-G);

  // Curvature vectors in the basis (xu,xv)
  // double v1[2] = {(fuu-fvv+sqrt(fuu*fvv*-2.0+fuu*fuu+(fuv*fuv)*4.0+fvv*fvv+(fu*fu)*(fuv*fuv)*4.0+(fuu*fuu)*(fv*fv)*2.0+(fuu*fuu)*(fv*fv*fv*fv)+(fu*fu)*(fvv*fvv)*2.0+(fuv*fuv)*(fv*fv)*4.0+(fu*fu*fu*fu)*(fvv*fvv)+(fu*fu)*(fuv*fuv)*(fv*fv)*4.0-(fu*fu)*fuu*fvv*2.0-fuu*(fv*fv)*fvv*2.0-fu*fuu*fuv*(fv*fv*fv)*4.0-(fu*fu*fu)*fuv*fv*fvv*4.0+(fu*fu)*fuu*(fv*fv)*fvv*2.0-fu*fuu*fuv*fv*4.0-fu*fuv*fv*fvv*4.0)+fuu*(fv*fv)-(fu*fu)*fvv)/(fuv*2.0+(fu*fu)*fuv*2.0-fu*fuu*fv*2.0), 1.};
  // double n1 = sqrt(v1[0]*v1[0]+v1[1]*v1[1]);
  // v1[0] /= n1;
  // v1[1] /= n1;
  // double v2[2] = {((-fuu+fvv+sqrt(fuu*fvv*-2.0+fuu*fuu+(fuv*fuv)*4.0+fvv*fvv+(fu*fu)*(fuv*fuv)*4.0+(fuu*fuu)*(fv*fv)*2.0+(fuu*fuu)*(fv*fv*fv*fv)+(fu*fu)*(fvv*fvv)*2.0+(fuv*fuv)*(fv*fv)*4.0+(fu*fu*fu*fu)*(fvv*fvv)+(fu*fu)*(fuv*fuv)*(fv*fv)*4.0-(fu*fu)*fuu*fvv*2.0-fuu*(fv*fv)*fvv*2.0-fu*fuu*fuv*(fv*fv*fv)*4.0-(fu*fu*fu)*fuv*fv*fvv*4.0+(fu*fu)*fuu*(fv*fv)*fvv*2.0-fu*fuu*fuv*fv*4.0-fu*fuv*fv*fvv*4.0)-fuu*(fv*fv)+(fu*fu)*fvv)*(-1.0/2.0))/(fuv+(fu*fu)*fuv-fu*fuu*fv), 1.};
  // double n2 = sqrt(v2[0]*v2[0]+v2[1]*v2[1]);
  // v2[0] /= n2;
  // v2[1] /= n2;

  SquareMatrix A(2);
  A(0,0) = -(fuu*fv*fv - fu*fuv*fv + fuu)/pow(fu*fu + fv*fv + 1, 3./2.);
  A(0,1) = -(fuv*fv*fv - fu*fvv*fv + fuv)/pow(fu*fu + fv*fv + 1, 3./2.);
  A(1,0) = -(fuv*fu*fu - fuu*fv*fu + fuv)/pow(fu*fu + fv*fv + 1, 3./2.);
  A(1,1) = -(fvv*fu*fu - fuv*fv*fu + fvv)/pow(fu*fu + fv*fv + 1, 3./2.);
  double l1, l2, v1[2], v2[2];
  A.getEigenvectorsAndEigenvalues(v1, v2, l1, l2);

  double T1[3] = {v1[0] * xu[0] + v1[1] * xv[0], v1[0] * xu[1] + v1[1] * xv[1], v1[0] * xu[2] + v1[1] * xv[2]};
  double nT1 = sqrt(T1[0]*T1[0]+T1[1]*T1[1]+T1[2]*T1[2]);
  T1[0] /= nT1;
  T1[1] /= nT1;
  T1[2] /= nT1;
  double T2[3] = {v2[0] * xu[0] + v2[1] * xv[0], v2[0] * xu[1] + v2[1] * xv[1], v2[0] * xu[2] + v2[1] * xv[2]};
  double nT2 = sqrt(T2[0]*T2[0]+T2[1]*T2[1]+T2[2]*T2[2]);
  T2[0] /= nT2;
  T2[1] /= nT2;
  T2[2] /= nT2;

  // Check orthogonality of (T1,T2,N)
  if(dotprod(T1,N) > 1e-10) { feInfo("T1 et N pas orthogonaux: %+-1.10e", dotprod(T1,N)); exit(-1); }
  if(dotprod(T2,N) > 1e-10) { feInfo("T2 et N pas orthogonaux: %+-1.10e", dotprod(T2,N)); exit(-1); }
  if(dotprod(T1,T2) > 1e-10) { feInfo("T1 et T2 pas orthogonaux: %+-1.10e", dotprod(T1,T2)); exit(-1); }

  double c = 0.1;
  double z0 = x[0]*x[0]+x[1]*x[1];
  fprintf(file, "SP(%f,%f,%f){1};\n", x[0], x[1], z0);
  fprintf(file, "SL(%f,%f,%f,%f,%f,%f){2,2};\n", x[0], x[1], z0, x[0]+c*T1[0], x[1]+c*T1[1], z0+c*T1[2]);
  fprintf(file, "SL(%f,%f,%f,%f,%f,%f){2,2};\n", x[0], x[1], z0, x[0]+c*T2[0], x[1]+c*T2[1], z0+c*T2[2]);
  fprintf(file, "SL(%f,%f,%f,%f,%f,%f){2,2};\n", x[0], x[1], z0, x[0]+c*N[0],  x[1]+c*N[1],  z0+c*N[2]);

  double eps = 1e-2;
  // double k1 = fmax(l1,l2);
  // double k2 = fmin(l1,l2);
  double L1 = pow( 2./ l2 * sqrt(eps*(2.-eps)), -2);
  double L2 = pow( 2./ l1 * sqrt(eps*(l1/l2)*(2.-eps*(l1/l2))), -2);
  double L3 = 1./(eps*eps);

  // The metric for P1 approximation of the graph surface
  M(0,0) = fu*(fu*(L3*N[2]*N[2] + L1*T1[2]*T1[2] + L2*T2[2]*T2[2]) + L3*N[0]*N[2] + L1*T1[0]*T1[2] + L2*T2[0]*T2[2]) + L3*N[0]*N[0] + L1*T1[0]*T1[0] + L2*T2[0]*T2[0] + fu*(L3*N[0]*N[2] + L1*T1[0]*T1[2] + L2*T2[0]*T2[2]);
  M(0,1) = fv*(fu*(L3*N[2]*N[2] + L1*T1[2]*T1[2] + L2*T2[2]*T2[2]) + L3*N[0]*N[2] + L1*T1[0]*T1[2] + L2*T2[0]*T2[2]) + fu*(L3*N[1]*N[2] + L1*T1[1]*T1[2] + L2*T2[1]*T2[2]) + L3*N[0]*N[1] + L1*T1[0]*T1[1] + L2*T2[0]*T2[1];
  M(1,0) = fu*(fv*(L3*N[2]*N[2] + L1*T1[2]*T1[2] + L2*T2[2]*T2[2]) + L3*N[1]*N[2] + L1*T1[1]*T1[2] + L2*T2[1]*T2[2]) + fv*(L3*N[0]*N[2] + L1*T1[0]*T1[2] + L2*T2[0]*T2[2]) + L3*N[0]*N[1] + L1*T1[0]*T1[1] + L2*T2[0]*T2[1];
  M(1,1) = fv*(fv*(L3*N[2]*N[2] + L1*T1[2]*T1[2] + L2*T2[2]*T2[2]) + L3*N[1]*N[2] + L1*T1[1]*T1[2] + L2*T2[1]*T2[2]) + L3*N[1]*N[1] + L1*T1[1]*T1[1] + L2*T2[1]*T2[1] + fv*(L3*N[1]*N[2] + L1*T1[1]*T1[2] + L2*T2[1]*T2[2]);

  double Xu = (v1[0] * (1 + fu*fu) + v1[1] * (    fu*fv))/nT1;
  double Xv = (v1[0] * (    fu*fv) + v1[1] * (1 + fv*fv))/nT1;
  double Yu = (v2[0] * (1 + fu*fu) + v2[1] * (    fu*fv))/nT2;
  double Yv = (v2[0] * (    fu*fv) + v2[1] * (1 + fv*fv))/nT2;

  double dXdu[2][2] = {{Xu, Xv},{Yu, Yv}};
  // dpsi/dX = dXdu^-1
  double det = dXdu[0][0]*dXdu[1][1] - dXdu[0][1]*dXdu[1][0];
  double dpsidX[2][2] = {{ dXdu[1][1]/det, -dXdu[0][1]/det},
                         {-dXdu[1][0]/det,  dXdu[0][0]/det}};

  double psiux = dpsidX[0][0];
  double psiuy = dpsidX[0][1];
  double psivx = dpsidX[1][0];
  double psivy = dpsidX[1][1];

  double Xuu = (v1[0] * fuu*fu + v1[1] * fuu*fv)/nT1;
  double Xuv = (v1[0] * fuv*fu + v1[1] * fuv*fv)/nT1;
  double Xvv = (v1[0] * fvv*fu + v1[1] * fvv*fv)/nT1;
  double Yuu = (v2[0] * fuu*fu + v2[1] * fuu*fv)/nT2;
  double Yuv = (v2[0] * fuv*fu + v2[1] * fuv*fv)/nT2;
  double Yvv = (v2[0] * fvv*fu + v2[1] * fvv*fv)/nT2;

  double d2Xdu2[2][2][2] = {{{Xuu, Xuv},{Xuv, Xvv}}, {{Yuu, Yuv},{Yuv, Yvv}}};

  double dpsi_udx[2][2] = {{0.,0.},{0.,0.}}, dpsi_vdx[2][2] = {{0.,0.},{0.,0.}};
  for(int a = 0; a < 2; ++a) {
    for(int b = 0; b < 2; ++b) {
      for(int j = 0; j < 2; ++j) {
        for(int k = 0; k < 2; ++k) {
          for(int m = 0; m < 2; ++m) {
            dpsi_udx[a][b] -= dpsidX[0][m] * d2Xdu2[m][j][k] * dXdu[k][b] * dXdu[j][a];
            dpsi_vdx[a][b] -= dpsidX[1][m] * d2Xdu2[m][j][k] * dXdu[k][b] * dXdu[j][a];
          }
        }
      }
    }
  }

  // Coefficients du développement de Z(X,Y) au second ordre
  // À comparer avec diag(k1, k2) ?

  // double Zu = 0.;
  // double Zv = 0.;

  double Zuu = dotprod(xuu,N);
  double Zuv = dotprod(xuv,N);
  double Zvv = dotprod(xvv,N);

  // double Zuuu = dotprod(xuuu,N);
  // double Zuuv = dotprod(xuuv,N);
  // double Zuvv = dotprod(xuvv,N);
  // double Zvvv = dotprod(xvvv,N);

  // double C[2][2] = {{0.,0.},{0.,0.}};
  // Linear terms of (u-u0)^2, (u-u0)*(v-v0) and (v-v0)^2
  // ----> PAS BESOIN SI Zu = Zv = 0?
  // C[0][0] =    psiux*psiux + psiux*psivx + psivx*psivx;
  // C[0][1] = (2.*psiux*psiuy + psiux*psivy + psivx*psiuy + 2.*psivx*psivy)/2.;
  // C[1][0] = (2.*psiux*psiuy + psiux*psivy + psivx*psiuy + 2.*psivx*psivy)/2.;
  // C[1][1] =    psiuy*psiuy + psiuy*psivy + psivy*psivy;

  // Quadratic terms in (u-u0) and (v-v0)
  // C[0][0] = Zuu * 0.5 * (dpsi_udx[0][0] + dpsi_vdx[0][0]);
  // C[0][1] = Zuv * 0.5 * (dpsi_udx[0][1] + dpsi_vdx[0][1]);
  // C[1][0] = Zuv * 0.5 * (dpsi_udx[0][1] + dpsi_vdx[0][1]);
  // C[1][1] = Zvv * 0.5 * (dpsi_udx[1][1] + dpsi_vdx[1][1]);

  // C[0][0] = Zuu * (psiux*psiux + psivx*psivx);
  // // C[0][1] = Zuv * 0.5 * (dpsi_udx[0][1] + dpsi_vdx[0][1]);
  // // C[1][0] = Zuv * 0.5 * (dpsi_udx[0][1] + dpsi_vdx[0][1]);
  // C[1][1] = Zvv * (psiuy*psiuy + psivy*psivy);

  SquareMatrix C(2);
  C(0,0) = Zuu * (psiux*psiux + psivx*psivx);
  C(0,1) = Zuv * (psiux*psiuy + psivx*psivy);
  C(1,0) = Zuv * (psiux*psiuy + psivx*psivy);
  C(1,1) = Zvv * (psiuy*psiuy + psivy*psivy);

  feInfo(" ");
  feInfo("k = %+-1.3e - %+-1.3e", l1, l2);
  feInfo("C = %+-1.3e - %+-1.3e - %+-1.3e", C(0,0), C(0,1), C(1,1));

  double lC1, lC2, vC1[2], vC2[2];
  C.getEigenvectorsAndEigenvalues(vC1, vC2, lC1, lC2);

  // Comparison with P1 metric
  L1 = pow( 2./ lC2 * sqrt(eps*(2.-eps)), -2);
  L2 = pow( 2./ lC1 * sqrt(eps*(2.-eps)), -2);
  M(0,0) = fu*(fu*(L3*N[2]*N[2] + L1*T1[2]*T1[2] + L2*T2[2]*T2[2]) + L3*N[0]*N[2] + L1*T1[0]*T1[2] + L2*T2[0]*T2[2]) + L3*N[0]*N[0] + L1*T1[0]*T1[0] + L2*T2[0]*T2[0] + fu*(L3*N[0]*N[2] + L1*T1[0]*T1[2] + L2*T2[0]*T2[2]);
  M(0,1) = fv*(fu*(L3*N[2]*N[2] + L1*T1[2]*T1[2] + L2*T2[2]*T2[2]) + L3*N[0]*N[2] + L1*T1[0]*T1[2] + L2*T2[0]*T2[2]) + fu*(L3*N[1]*N[2] + L1*T1[1]*T1[2] + L2*T2[1]*T2[2]) + L3*N[0]*N[1] + L1*T1[0]*T1[1] + L2*T2[0]*T2[1];
  M(1,0) = fu*(fv*(L3*N[2]*N[2] + L1*T1[2]*T1[2] + L2*T2[2]*T2[2]) + L3*N[1]*N[2] + L1*T1[1]*T1[2] + L2*T2[1]*T2[2]) + fv*(L3*N[0]*N[2] + L1*T1[0]*T1[2] + L2*T2[0]*T2[2]) + L3*N[0]*N[1] + L1*T1[0]*T1[1] + L2*T2[0]*T2[1];
  M(1,1) = fv*(fv*(L3*N[2]*N[2] + L1*T1[2]*T1[2] + L2*T2[2]*T2[2]) + L3*N[1]*N[2] + L1*T1[1]*T1[2] + L2*T2[1]*T2[2]) + L3*N[1]*N[1] + L1*T1[1]*T1[1] + L2*T2[1]*T2[1] + fv*(L3*N[1]*N[2] + L1*T1[1]*T1[2] + L2*T2[1]*T2[2]);

  // // VERIFIER SI IL FAUT AJOUTER LES VECTEURS PROPRES DE C: PROBABLEMENT
  // // CONTINUER ET TESTER LA MÉTRIQUE AVEC TERMES P3: CALCULER C AVEC LOG-SIMPLEXE
  // MetricTensor Q(1.);
  // // computeAnalyticMetricP2ForLpNorm(coeffs, Q, 100.);

  // SquareMatrix C2(2);
  // C2(0,0) = 1e-2;
  // C2(0,1) = 0.;
  // C2(1,0) = 0.;
  // C2(1,1) = 1e-2;

  // double lC21, lC22, vC21[2], vC22[2];
  // C2.getEigenvectorsAndEigenvalues(vC21, vC22, lC21, lC22);

  // // Comparison with P1 metric
  // eps = 1e-7;
  // L1 = pow( 2./ lC22 * sqrt(eps*(2.-eps)), -2);
  // L2 = pow( 2./ lC21 * sqrt(eps*(2.-eps)), -2);
  // L3 = 1./(eps*eps);
  // M(0,0) = fu*(fu*(L3*N[2]*N[2] + L1*T1[2]*T1[2] + L2*T2[2]*T2[2]) + L3*N[0]*N[2] + L1*T1[0]*T1[2] + L2*T2[0]*T2[2]) + L3*N[0]*N[0] + L1*T1[0]*T1[0] + L2*T2[0]*T2[0] + fu*(L3*N[0]*N[2] + L1*T1[0]*T1[2] + L2*T2[0]*T2[2]);
  // M(0,1) = fv*(fu*(L3*N[2]*N[2] + L1*T1[2]*T1[2] + L2*T2[2]*T2[2]) + L3*N[0]*N[2] + L1*T1[0]*T1[2] + L2*T2[0]*T2[2]) + fu*(L3*N[1]*N[2] + L1*T1[1]*T1[2] + L2*T2[1]*T2[2]) + L3*N[0]*N[1] + L1*T1[0]*T1[1] + L2*T2[0]*T2[1];
  // M(1,0) = fu*(fv*(L3*N[2]*N[2] + L1*T1[2]*T1[2] + L2*T2[2]*T2[2]) + L3*N[1]*N[2] + L1*T1[1]*T1[2] + L2*T2[1]*T2[2]) + fv*(L3*N[0]*N[2] + L1*T1[0]*T1[2] + L2*T2[0]*T2[2]) + L3*N[0]*N[1] + L1*T1[0]*T1[1] + L2*T2[0]*T2[1];
  // M(1,1) = fv*(fv*(L3*N[2]*N[2] + L1*T1[2]*T1[2] + L2*T2[2]*T2[2]) + L3*N[1]*N[2] + L1*T1[1]*T1[2] + L2*T2[1]*T2[2]) + L3*N[1]*N[1] + L1*T1[1]*T1[1] + L2*T2[1]*T2[1] + fv*(L3*N[1]*N[2] + L1*T1[1]*T1[2] + L2*T2[1]*T2[2]);
}