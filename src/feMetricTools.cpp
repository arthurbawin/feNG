#include "feMetricTools.h"

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
#include "../contrib/Eigen/Eigen"
#include "../contrib/unsupported/Eigen/MatrixFunctions"

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
    MetricType &m = metrics[nodeTags[i]];
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];
    getEllipsePoints(factor * m(0, 0), factor * 2.0 * m(0, 1), factor * m(1, 1), x, y, xP, yP);
    for(int i = 0; i < nt; ++i) {
      if(i != nt - 1) {
        fprintf(f, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[i], yP[i], 0.,
                xP[i + 1], yP[i + 1], 0., 1, 1);
      } else {
        fprintf(f, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[i], yP[i], 0., xP[0],
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

void gradationMetriques(double gradation, int nmax, std::vector<double> &coord1,
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

  for(auto val : elementTypes) feInfo("Gradation: Elements classified on dim 2 are of type %d", val);
  feInfo("Gradation: There are %d elements of dim 2 and type [0]", elementTags[0].size());
  feInfo("Gradation: There are %d nodes total in elem [0]", nodeTags.size());

  int nodedim, nodetag;

  FILE *fff = fopen("CHECK.pos", "w");
  fprintf(fff, "View \" CHECK \"{\n");

  std::pair<std::set<std::pair<size_t, size_t>, gmshEdgeLessThan>::iterator, bool> ret;

  for(size_t i = 0; i < elementTags[0].size(); i++) {
    size_t n0 = elemNodeTags[0][3 * i + 0];
    size_t n1 = elemNodeTags[0][3 * i + 1];
    size_t n2 = elemNodeTags[0][3 * i + 2];
    ret = edges.insert(std::make_pair(n0, n1));
    ret = edges.insert(std::make_pair(n1, n2));
    ret = edges.insert(std::make_pair(n2, n0));
    gmsh::model::mesh::getNode(n0, nodecoord1, nodeparametricCoord, nodedim, nodetag);
    gmsh::model::mesh::getNode(n1, nodecoord2, nodeparametricCoord, nodedim, nodetag);
    gmsh::model::mesh::getNode(n2, nodecoord3, nodeparametricCoord, nodedim, nodetag);
    double x1 = nodecoord1[0];
    double y1 = nodecoord1[1];
    double x2 = nodecoord2[0];
    double y2 = nodecoord2[1];
    double x3 = nodecoord3[0];
    double y3 = nodecoord3[1];
    fprintf(fff, "SL(%.16g,%.16g,0.,%.16g,%.16g,0.){%u, %u};\n", x1, y1, x2, y2, 1, 1);
    fprintf(fff, "SL(%.16g,%.16g,0.,%.16g,%.16g,0.){%u, %u};\n", x2, y2, x3, y3, 1, 1);
    fprintf(fff, "SL(%.16g,%.16g,0.,%.16g,%.16g,0.){%u, %u};\n", x3, y3, x1, y1, 1, 1);
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

    // FILE* fff = fopen("CHECK.pos","w");
    // fprintf(fff, "View \" CHECK \"{\n");

    // size_t n0 = elemNodeTags[0][400];
    // gmsh::model::mesh::getNode(n0, nodecoord3, nodeparametricCoord, nodedim, nodetag);
    // double X0 = nodecoord3[0];
    // double Y0 = nodecoord3[1];
    // SMetric3 mref = metricsOnGmshModel[n0];
    // fprintf(fff, "SP(%.16g,%.16g,0.){%u};\n", X0, Y0, 1);

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

    feInfo("Gradation: Passe %d", iter);
    // correction = false;

    // fprintf(fff, "};");
    // fclose(fff);
  }
  // feInfo("%d coord size", coord.size());
  // feInfo("%d edges size", edges.size());
#endif
}

Eigen::VectorXd CMINUS(7);
Eigen::VectorXd CPLUS(7);

double solveErrorFunction(double k, double *v, double *w, double H[2][2], double C[2][2][2],
                                 std::set<double> &roots)
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

  double C111 = C[0][0][0];
  double C112 = (C[0][0][1] + C[0][1][0] + C[1][0][0]) / 3.;
  double C122 = (C[0][1][1] + C[1][0][1] + C[1][1][0]) / 3.;
  double C222 = C[1][1][1];

  double cs6 = (C111 * a1 * a1 * a1 + 3. * C112 * a1 * a1 * a2 + 3. * C122 * a1 * a2 * a2 +
                C222 * a2 * a2 * a2) /
               120.;
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

  CMINUS(0) = cs6;
  CMINUS(1) = cs5;
  CMINUS(2) = cs4;
  CMINUS(3) = cs3;
  CMINUS(4) = 0.;
  CMINUS(5) = 0.;
  CMINUS(6) = 1.;

  CPLUS(0) = cs6;
  CPLUS(1) = cs5;
  CPLUS(2) = cs4;
  CPLUS(3) = cs3;
  CPLUS(4) = 0.;
  CPLUS(5) = 0.;
  CPLUS(6) = -1.;

  double s = 1e10;
  roots = RootFinder::solvePolynomial(CMINUS, 0., INFINITY, 1e-8);
  for(auto it = roots.begin(); it != roots.end(); it++) {
    s = fmin(s, *it);
  }

  roots = RootFinder::solvePolynomial(CPLUS, 0., INFINITY, 1e-8);
  for(auto it = roots.begin(); it != roots.end(); it++) {
    s = fmin(s, *it);
  }

  // feInfo("lowest root is %1.16e", s);
  return s;
}

#if defined(HAVE_SOPLEX)

// MetricTensor Hij = MetricTensor(1.0);

double evaluateHomogeneousErrorPolynomial(const std::vector<double> &errorCoefficients,
                                                 const int degree, const double xLoc, const double yLoc, MetricTensor &Hij)
{
	// Using homogeneous polynomial factorization.
	// See e.g. Coulaud & Loseille, Very High Order Anisotropic Metric-Based Mesh Adaptation in 3D
	// Proposition 2.2. The matrix H has been adapted for 2D polynomial (no factor 2 on non-diagonal entries)
  double res1 = 0.;
  for(int i = degree-1, j = 0, cnt = 0; i >= 0; --i, ++j, ++cnt){
  	if(i+j == degree-1){
  		
	  	Hij(0,0) = errorCoefficients[ cnt + 0 + 0 ];
	  	Hij(0,1) = errorCoefficients[ cnt + 0 + 1 ];
	  	Hij(1,0) = errorCoefficients[ cnt + 1 + 0 ];
	  	Hij(1,1) = errorCoefficients[ cnt + 1 + 1 ];
	  	Hij = Hij.absoluteValueEigen();

	  	res1 += fabs(pow(xLoc, i) * pow(yLoc, j)) * 
	  		  (xLoc*xLoc * Hij(0,0) + 2. * xLoc*yLoc * Hij(0,1) + yLoc*yLoc * Hij(1,1));
	  }
  }
  return res1;

  // Classical polynomial evaluation
  double res2 = 0.;
  int expx = errorCoefficients.size() - 1, expy = 0;

  // Quick and dirty binomial coefficients
  int coeff[4];
  if(degree == 1) {
    coeff[0] = 1.;
    coeff[1] = 2.;
    coeff[2] = 1.;
  } else if(degree == 2) {
    coeff[0] = 1.;
    coeff[1] = 3.;
    coeff[2] = 3.;
    coeff[3] = 1.;
  } else {
    feInfo("Implement binomial coefficients for P3+ recoveries");
    exit(-1);
  }

  for(int i = 0; i < errorCoefficients.size(); i++, expx--, expy++) {
    res2 += coeff[i] * errorCoefficients[i] * pow(xLoc, expx) * pow(yLoc, expy);
  }
  return res2;
}

//
// Compute a discretization of the error curve 1 for straight aniso adaptation
//
void computeLvl1(const std::vector<double> &errorCoefficients, const int degree,
                 linearProblem &myLP)
{
  double theta, xi, yi, error, root;

  for(size_t i = 0; i < myLP.numConstraints; ++i) {
    theta = 2. * M_PI * (double) i / (double) myLP.numConstraints;
    xi = cos(theta);
    yi = sin(theta);
    error = fmax(1e-14, fabs(evaluateHomogeneousErrorPolynomial(errorCoefficients, degree, xi, yi, myLP.Hij)));
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
        error = fmax(1e-14, fabs(evaluateHomogeneousErrorPolynomial(errorCoefficients, degree, xi, yi, myLP.Hij)));
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
void computeLvl1(double K1, double K2, double *g1, double *g2, double Hij[2][2],
                        double Cijk[2][2][2], int nPhi, linearProblem &myLP)
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

    for(int i = 0; i < nPhi; ++i) {
      theta = M_PI / 2. * (double)i / (double)nPhi;
      k = k1 * cos(theta) * cos(theta) + k2 * sin(theta) * sin(theta);
      v[0] = v1[0] * cos(theta) + v2[0] * sin(theta);
      v[1] = v1[1] * cos(theta) + v2[1] * sin(theta);
      w[0] = -v[1];
      w[1] = v[0];

      s = solveErrorFunction(k, v, w, Hij, Cijk, roots);

      // Shoot :
      myLP.lvl1[iDir * 2 * nPhi + 2 * i + 0] = x0[0] + s * v[0] + s * s / 2. * k * w[0];
      myLP.lvl1[iDir * 2 * nPhi + 2 * i + 1] = x0[1] + s * v[1] + s * s / 2. * k * w[1];
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
    myLP.lprow.setLhs(lhs);
    myLP.lprow.setRowVector(myLP.row);
    myLP.lprowset.add(myLP.lprow);
    myLP.row.clear();
  }
  myLP.problem.removeRowRangeReal(0, myLP.numConstraints - 1);
  myLP.problem.addRowsReal(myLP.lprowset);
  myLP.lprowset.clear();

  /* Solve LP */
  myLP.stat = myLP.problem.optimize();

  if(myLP.stat == SPxSolver::OPTIMAL) {
    myLP.problem.getPrimal(myLP.prim);
    L(0, 0) = myLP.prim[0];
    L(0, 1) = myLP.prim[1];
    L(1, 0) = myLP.prim[1];
    L(1, 1) = myLP.prim[2];
    return true;
  } else {
    switch(myLP.stat) {
      case SPxSolver::ERROR: feInfo("Soplex status: an error occured"); break;
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
                                     const int maxIter, const double tol, MetricTensor &Qres,
                                     int &numIter, linearProblem &myLP)
{
  #if defined(LOG_SIMPLEX_DEBUG)
	FILE *ff = fopen("lvl1.pos", "w");
  fprintf(ff, "View\"lvl1\"{\n");
  #endif

  computeLvl1(errorCoefficients, degree, myLP);

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
  FILE *myfile = fopen("logSimplexDebug.pos", "w");
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

bool logSimplexCurved(double *x, double k1, double k2, double g1[2], double g2[2],
                             double Hij[2][2], double Cijk[2][2][2], int maxIter,
                             int nThetaPerQuadrant, double tol, Eigen::Matrix2d &Q, int &numIter,
                             linearProblem &myLP)
{
  int nPhi = nThetaPerQuadrant;
  // int size = 2 * myLP.numConstraints; // 2 coordonnees * 4 quadrants * nPhi angles/quadrant
  // std::vector<double> lvl1(size, 0.), constraints(size, 0.);
  computeLvl1(k1, k2, g1, g2, Hij, Cijk, nPhi, myLP);

  Eigen::Matrix2d Qprev, Q12, L, diff;

  Q = Eigen::Matrix2d::Identity();

  double xi, yi;
  bool res;

  for(int iter = 0; iter < maxIter; ++iter) {
    Qprev = Q;
    Q12 = Q.sqrt();

    // Apply transformation Q12 to the initial constraints
    for(int i = 0; i < myLP.numConstraints; ++i) {
      xi = myLP.lvl1[2 * i];
      yi = myLP.lvl1[2 * i + 1];
      myLP.constraints[2 * i] = Q12(0, 0) * xi + Q12(0, 1) * yi;
      myLP.constraints[2 * i + 1] = Q12(1, 0) * xi + Q12(1, 1) * yi;
    }

    // Solve the linear optimization problem for L
    res = solveLP(myLP, L);

    if(res) {
      // Recover Q from L
      Q = Q12 * L.exp() * Q12;

      diff = Q - Qprev;
      if(diff.norm() < tol) {
        numIter = iter;
        break;
      }
    } else {
      // LP solver returned an error : return the last valid Q
      if(iter == 0) {
        feWarning("In logSimplexCurved : LP solver returned an error at iteration %d. Returning "
                  "identity matrix.",
                  iter);
        FILE *f = fopen("lvl1_check.pos", "w");
        fprintf(f, "View\"lvl1\"{\n");
        for(int i = 0; i < myLP.numConstraints; ++i) {
          // feInfo("%f - %f", myLP.lvl1[2*i], myLP.lvl1[2*i+1]);
          fprintf(f, "SP(%g,%g,0){%g,%g,0};", x[0] + myLP.lvl1[2 * i], x[1] + myLP.lvl1[2 * i + 1],
                  1., 1.);
        }
        int nt = 30;
        std::vector<double> xP(nt, 0.);
        std::vector<double> yP(nt, 0.);
        double factor = 1.;
        getEllipsePoints(factor * Q(0, 0), factor * 2.0 * Q(0, 1), factor * Q(1, 1), x[0], x[1], xP,
                         yP);
        for(int j = 0; j < nt; ++j) {
          if(j != nt - 1) {
            fprintf(f, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0.,
                    xP[j + 1], yP[j + 1], 0., 1, 1);
          } else {
            fprintf(f, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0.,
                    xP[0], yP[0], 0., 1, 1);
          }
        }
        fprintf(f, "};\n");
        fclose(f);
        return false;
      }

      return true;
    }
  }
  return true;
}
#endif