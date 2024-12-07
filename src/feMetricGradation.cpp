
#include "feMetric.h"
#include "feMatrixInterface.h"
#include <algorithm>
#include <random>

#if defined(HAVE_GMSH)
#include "gmsh.h"
#endif

extern int FE_VERBOSE;

#define TOL_ZERO 1e-15

struct gmshEdgeLessThan {
  bool operator()(const std::pair<int, int> &e1, const std::pair<int, int> &e2) const
  {
    int tag10 = e1.first, tag11 = e1.second;
    int tag20 = e2.first, tag21 = e2.second;
    int diffMin = fmin(tag10, tag11) - fmin(tag20, tag21);
    int diffMax = fmax(tag10, tag11) - fmax(tag20, tag21);
    if(diffMin < 0) return true;
    if(diffMin > 0) return false;
    if(diffMax < 0) return true;
    if(diffMax > 0) return false;
    return false;
  }
};

inline MetricTensor spanMetric(const GradationSpace space, const double gradation, 
	const MetricTensor &M, const double pq[2])
{
	switch(space) {
    case GradationSpace::Metric:
		  // Homogeneous growth in the metric space (more anisotropic)
	    return M.spanMetricInMetricSpace(gradation, pq);
    case GradationSpace::Physical:
  	  // Homogeneous growth in the physical space (becomes isotropic far from p)
		  return M.spanMetricInPhysicalSpace(gradation, pq);
    case GradationSpace::Mixed:
      // Mixed space with t = 1/8
      return M.spanMetricInMixedSpace(gradation, pq, 0.125);
    case GradationSpace::ExpMetric:
      return M.spanMetricInMetricSpace(gradation, pq);
    default:
      return M.spanMetricInMetricSpace(gradation, pq);
  }
}

double matNorm2(const MetricTensor &m)
{
  double sqr = 0;
  for(int i = 0; i < 2; i++) {
    for(int j = 0; j < 2; j++) {
      sqr += m(i, j) * m(i, j);
    }
  }
  return sqrt(sqr);
}

double matNorm2Diff(const MetricTensor &m1, const MetricTensor &m2)
{
  double sqr = 0;
  for(int i = 0; i < 2; i++) {
    for(int j = 0; j < 2; j++) {
      sqr += (m1(i, j) - m2(i, j)) * (m1(i, j) - m2(i, j));
    }
  }
  return sqrt(sqr);
}

// Metrics are modified if || M - M_grown ||/norm(M) >= REL_TOL_GRADATION
#define REL_TOL_GRADATION 1e-1
// If fabs(non-diag)/norm(M) is below REL_TOL_DIAGONAL, matrix is diag
#define REL_TOL_DIAGONAL 1e-4
#define REL_TOL_EQUAL 1e-2

double EX[2] = {1., 0.};
double EY[2] = {0., 1.};
thread_local double EV[2], V0[2], V1[2]; // eigenvalues

#define PRINT(name) printScalar2(#name, (name))
void printScalar2(const char *name, double val) {
  feInfo("%s = %+-1.16e", name, val);
} 

//
// Compute the intersection of metrics m1 and m2 using their simultaneous reduction inv(m1)*m2
//
MetricTensor intersectionReductionSimultanee(const MetricTensor &m1, const MetricTensor &m2)
{
  double a1 = m1(0, 0), b1 = m1(0, 1), c1 = m1(1, 1);
  double a2 = m2(0, 0), b2 = m2(0, 1), c2 = m2(1, 1);
  double norm1 = matNorm2(m1);
  double norm2 = matNorm2(m2);
  double v00, v01, vTInv00, vTInv01, vTInv10, vTInv11;

  double relDiff = matNorm2Diff(m1, m2) / fmin(norm1, norm2);

  bool matricesAreEqual = relDiff < REL_TOL_EQUAL;

  // If metrics are identical, return one of them
  if(matricesAreEqual) {
    // feInfo("Matrices are equal");
    return m1.copy();
  }

  bool m1IsDiagonal = fabs(b1)/norm1 < REL_TOL_DIAGONAL;
  bool m2IsDiagonal = fabs(b2)/norm2 < REL_TOL_DIAGONAL;

  if(m1IsDiagonal && m2IsDiagonal) {
    // Both matrices are diagonal
    // feInfo("Both diagonal");
    EV[0] = fmax(a1, a2);
    EV[1] = fmax(c1, c2);
    return MetricTensor(EV, EX, EY);

  } else if(m1IsDiagonal) {
    // feInfo("m1 is diagonal");
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

    V0[0] = v00;
    V0[1] = 1.;
    V1[0] = v01;
    V1[1] = 1.;
    double l0 = fmax(m1.dotProduct(V0, V0), m2.dotProduct(V0, V0));
    double l1 = fmax(m1.dotProduct(V1, V1), m2.dotProduct(V1, V1));

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
    	feInfo("Metric intersection - Error in case 1:");
      m1.print();
      m2.print();
      std::cout << "NAN alert : " << vTInv00 << " " << vTInv01 << " " << vTInv10 << " " << vTInv11
                << std::endl;
      exit(-1);
    }

    V0[0] = vTInv00;
    V0[1] = vTInv01;
    V1[0] = vTInv10;
    V1[1] = vTInv11;
    EV[0] = l0;
    EV[1] = l1;
    return MetricTensor(EV, V0, V1);

  } else if(m2IsDiagonal) {
    // feInfo("m2 is diagonal");
    v00 = (a1 * c2) / (a2 * b1) - (a1 * c2 + a2 * c1 +
                                   sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 +
                                        a2 * a2 * c1 * c1 + 4 * a2 * b1 * b1 * c2)) /
                                    (2 * a2 * b1);
    v01 = (a1 * c2) / (a2 * b1) - (a1 * c2 + a2 * c1 -
                                   sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 +
                                        a2 * a2 * c1 * c1 + 4 * a2 * b1 * b1 * c2)) /
                                    (2 * a2 * b1);

    // PRINT(v00);
    // PRINT(v01);

    V0[0] = v00;
    V0[1] = 1.;
    V1[0] = v01;
    V1[1] = 1.;

    // PRINT(m1.dotProduct(V0, V0));
    // PRINT(m2.dotProduct(V0, V0));
    // PRINT(m1.dotProduct(V1, V1));
    // PRINT(m2.dotProduct(V1, V1));

    // m1.print();
    // m2.print();

    double l0 = fmax(m1.dotProduct(V0, V0), m2.dotProduct(V0, V0));
    double l1 = fmax(m1.dotProduct(V1, V1), m2.dotProduct(V1, V1));

    // PRINT(l0);
    // PRINT(l1);

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

    // PRINT(vTInv00);
    // PRINT(vTInv01);
    // PRINT(vTInv10);
    // PRINT(vTInv11);

    if(isnan(vTInv00) || isnan(vTInv01) || isnan(vTInv10) || isnan(vTInv11) || isinf(vTInv00) ||
       isinf(vTInv01) || isinf(vTInv10) || isinf(vTInv11)) {
      std::cout << relDiff << std::endl;
      feInfo("Metric intersection - Error in case 2:");
      m1.print();
      m2.print();
      std::cout << "NAN alert : " << vTInv00 << " " << vTInv01 << " " << vTInv10 << " " << vTInv11
                << std::endl;
      exit(-1);
    }

    V0[0] = vTInv00;
    V0[1] = vTInv01;
    V1[0] = vTInv10;
    V1[1] = vTInv11;
    EV[0] = l0;
    EV[1] = l1;
    return MetricTensor(EV, V0, V1);

  } else {
    // feInfo("None is diagonal");
    // Check if metrics are multiple of one another
    if(fabs(a2 / a1 - c2 / c1) < 1e-3 && fabs(a2 / a1 - b2 / b1) < 1e-3) {
      // feInfo("Multiple of one another");
      return (a2 / a1 <= 1) ? m1.copy() : m2.copy();
    } else {
      // feInfo("Not multiple of one another");
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

      V0[0] = v00;
	    V0[1] = 1.;
	    V1[0] = v01;
	    V1[1] = 1.;
	    double l0 = fmax(m1.dotProduct(V0, V0), m2.dotProduct(V0, V0));
	    double l1 = fmax(m1.dotProduct(V1, V1), m2.dotProduct(V1, V1));

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
        feInfo("Metric intersection - Error in case 3:");
	      m1.print();
	      m2.print();
        std::cout << "NAN alert : " << vTInv00 << " " << vTInv01 << " " << vTInv10 << " " << vTInv11
                  << std::endl;
        exit(-1);
      }

      V0[0] = vTInv00;
	    V0[1] = vTInv01;
	    V1[0] = vTInv10;
	    V1[1] = vTInv11;
	    EV[0] = l0;
	    EV[1] = l1;
	    return MetricTensor(EV, V0, V1);
    }
  }
}

bool gradationOnEdge(const GradationSpace space, const double gradation, const std::pair<size_t, size_t> &edge, 
	const std::vector<double> &coord, std::map<int, MetricTensor> &metricsAtNodeTags)
{
	bool metricChanged = false;
  size_t np = edge.first; 
  size_t nq = edge.second;
  double xp[2] = {coord[3*(np-1)], coord[3*(np-1)+1]};
  double xq[2] = {coord[3*(nq-1)], coord[3*(nq-1)+1]};
  double pq[2] = {xq[0] - xp[0], xq[1] - xp[1]};
  double qp[2] = {xp[0] - xq[0], xp[1] - xq[1]};
  MetricTensor &Mp = metricsAtNodeTags.at(np);
  MetricTensor &Mq = metricsAtNodeTags.at(nq);

  // FILE *pfile = fopen("mp.pos", "w"); fprintf(pfile, "View \" mp \"{\n");
  // FILE *qfile = fopen("mq.pos", "w"); fprintf(qfile, "View \" mq \"{\n");
  // FILE *psfile = fopen("mpatq.pos", "w"); fprintf(psfile, "View \" mpatq \"{\n");
  // FILE *pifile = fopen("p_inter.pos", "w"); fprintf(pifile, "View \" p_inter \"{\n");
  // FILE *qsfile = fopen("mqatp.pos", "w"); fprintf(qsfile, "View \" mqatp \"{\n");
  // FILE *qifile = fopen("q_inter.pos", "w"); fprintf(qifile, "View \" q_inter \"{\n");

  // drawSingleEllipse(pfile, xp, Mp, 10., 100);
  // drawSingleEllipse(qfile, xq, Mq, 10., 100);

  // MetricTensor MpAtq = spanMetric(space, gradation, Mp, pq);
  // if(MpAtq.determinant() < 1e-10) {
  //   feInfo("1 - Determinant nul ou negatif: %+-1.4e", MpAtq.determinant());
  //   MpAtq.print();
  //   exit(-1);
  // }
  // drawSingleEllipse(psfile, xq, MpAtq, 10., 100);
  // MetricTensor MqAtp = spanMetric(space, gradation, Mq, qp);
  // if(MqAtp.determinant() < 1e-10) {
  //   feInfo("2 - Determinant nul ou negatif: %+-1.4e", MqAtp.determinant());
  //   MqAtp.print();
  //   exit(-1);
  // }
  // drawSingleEllipse(qsfile, xp, MqAtp, 10., 100);

  // Span Mp to q, intersect and check if Mq needs to be reduced
  MetricTensor MpAtq = spanMetric(space, gradation, Mp, pq);
  MpAtq = intersectionReductionSimultanee(Mq, MpAtq);
  ////////////////////
  if(MpAtq.determinant() < TOL_ZERO) {

    // Skip this edge
    metricChanged = false;
    return metricChanged;

    // feInfo("3 - Determinant nul ou negatif: %+-1.4e", MpAtq.determinant());
    // Mq.print();
    // MpAtq.print();

    // feInfo("Compare with previous function:");
    // SMetric3 mq(1.), mpatq(1.);
    // mq(0,0) = Mq(0,0);
    // mq(0,1) = Mq(0,1);
    // mq(1,0) = Mq(1,0);
    // mq(1,1) = Mq(1,1);
    // feInfo("Check mq:");
    // mq.print("mq");
    // Mq.print();
    // mpatq(0,0) = MpAtq(0,0);
    // mpatq(0,1) = MpAtq(0,1);
    // mpatq(1,0) = MpAtq(1,0);
    // mpatq(1,1) = MpAtq(1,1);
    // feInfo("Check mpatq:");
    // mpatq.print("mpatq");
    // MpAtq.print();
    // SMetric3 res = intersectionReductionSimultaneeExplicite(mq, mpatq);
    // res.print("res");

    // fprintf(pfile, "};"); fclose(pfile);
    // fprintf(qfile, "};"); fclose(qfile);
    // fprintf(psfile, "};"); fclose(psfile);
    // fprintf(pifile, "};"); fclose(pifile);
    // fprintf(qsfile, "};"); fclose(qsfile);
    // fprintf(qifile, "};"); fclose(qifile);

    exit(-1);
  }
  // drawSingleEllipse(qifile, xq, MpAtq, 10., 100);
  ////////////////////
  if(matNorm2Diff(Mq, MpAtq) / matNorm2(Mq) > REL_TOL_GRADATION) {
    Mq.assignMatrixFrom(MpAtq);
    metricChanged = true;
  }

  // // Span Mq to p
  MetricTensor MqAtp = spanMetric(space, gradation, Mq, qp);
  MqAtp = intersectionReductionSimultanee(Mp, MqAtp);
  ////////////////////
  if(MqAtp.determinant() < TOL_ZERO) {

     // Skip this edge
    metricChanged = false;
    return metricChanged;

    // feInfo("4 - Determinant nul ou negatif: %+-1.4e", MqAtp.determinant());
    // MqAtp.print();
    // exit(-1);
  }
  // drawSingleEllipse(pifile, xp, MqAtp, 10., 100);
  ////////////////////
  if(matNorm2Diff(Mp, MqAtp) / matNorm2(Mp) > REL_TOL_GRADATION) {
    Mp.assignMatrixFrom(MqAtp);
    metricChanged = true;
  }

  // fprintf(pfile, "};"); fclose(pfile);
  // fprintf(qfile, "};"); fclose(qfile);
  // fprintf(psfile, "};"); fclose(psfile);
  // fprintf(pifile, "};"); fclose(pifile);
  // fprintf(qsfile, "};"); fclose(qsfile);
  // fprintf(qifile, "};"); fclose(qifile);

  return metricChanged;
}

bool gradationOnEdgeExpMetric(const GradationSpace space, const double gradation, const std::pair<size_t, size_t> &edge, 
  const std::vector<double> &coord, const std::map<int, MetricTensor> &metricsAtNodeTags, std::map<int, MetricTensor> &mNew)
{
  bool metricChanged = false;
  size_t np = edge.first; 
  size_t nq = edge.second;
  double xp[2] = {coord[3*(np-1)], coord[3*(np-1)+1]};
  double xq[2] = {coord[3*(nq-1)], coord[3*(nq-1)+1]};
  double pq[2] = {xq[0] - xp[0], xq[1] - xp[1]};
  double qp[2] = {xp[0] - xq[0], xp[1] - xq[1]};
  const MetricTensor &Mp = metricsAtNodeTags.at(np);
  const MetricTensor &Mq = metricsAtNodeTags.at(nq);
  MetricTensor &Mnewp = mNew.at(np);
  MetricTensor &Mnewq = mNew.at(nq);

  // Span Mp to q, intersect and check if Mq needs to be reduced
  MetricTensor MpAtq = spanMetric(space, gradation, Mp, pq);
  MpAtq = intersectionReductionSimultanee(Mnewq, MpAtq);
  if(MpAtq.determinant() < TOL_ZERO) {
     // Skip this edge
    metricChanged = false;
    return metricChanged;
    // feInfo("3 - Determinant nul ou negatif: %+-1.4e", MpAtq.determinant());
    // exit(-1);
  }
  if(matNorm2Diff(Mnewq, MpAtq) / matNorm2(Mnewq) > REL_TOL_GRADATION) {
    Mnewq.assignMatrixFrom(MpAtq);
    metricChanged = true;
  }

  // Span Mq to p
  MetricTensor MqAtp = spanMetric(space, gradation, Mq, qp);
  MqAtp = intersectionReductionSimultanee(Mnewp, MqAtp);
  if(MqAtp.determinant() < TOL_ZERO) {
     // Skip this edge
    metricChanged = false;
    return metricChanged;
    // feInfo("4 - Determinant nul ou negatif: %+-1.4e", MqAtp.determinant());
    // exit(-1);
  }
  if(matNorm2Diff(Mnewp, MqAtp) / matNorm2(Mnewp) > REL_TOL_GRADATION) {
    Mnewp.assignMatrixFrom(MqAtp);
    metricChanged = true;
  }

  return metricChanged;
}

// Apply gradation to the metric tensor field.
// Gradation is performed on the Gmsh model, edge by edge according to Alauzet's paper.
feStatus feMetric::newGradation(const std::vector<std::size_t> &nodeTags, const std::vector<double> &coord,
  const double gradation, std::map<int, MetricTensor> &metricsAtNodeTags)
{
#if defined(HAVE_GMSH)

	int dim = 2;
  std::vector<int> elementTypes;
  std::vector<std::vector<std::size_t> > elementTags;
  std::vector<std::vector<std::size_t> > elemNodeTags;
  gmsh::model::mesh::getElements(elementTypes, elementTags, elemNodeTags, dim);

  if(elementTypes[0] != 2 && elementTypes[0] != 9) {
  	return feErrorMsg(FE_STATUS_ERROR, "Can only apply gradation for triangle meshes!");
  }

  size_t numNodesPerElem = (elementTypes[0] == 2) ? 3 : 6;

  // feInfoCond(FE_VERBOSE > 0, "METRIC GRADATION: Elements classified on dim 2 are of type %d", elementTypes[0]);

  FILE *myFile = fopen("edgesForGradation.pos", "w");
  fprintf(myFile, "View \" edgesForGradation \"{\n");

  std::set<std::pair<size_t, size_t>, gmshEdgeLessThan> edges;

  // All sub-edges in P2 triangle
  size_t nSubEdges = 15;
  size_t e0[] = {0,0,0,0,0,1,1,1,1,2,2,2,3,3,4};
  size_t e1[] = {1,2,3,4,5,2,3,4,5,3,4,5,4,5,5};

  for(size_t i = 0; i < elementTags[0].size(); i++) {
    if(numNodesPerElem == 3) {
      // P1 triangles - Add each edge
      for(size_t j = 0; j < 3; ++j) {
        size_t n0 = elemNodeTags[0][numNodesPerElem * i + j];
        size_t n1 = elemNodeTags[0][numNodesPerElem * i + (j+1) % numNodesPerElem];
        edges.insert(std::make_pair(n0, n1));
        fprintf(myFile, "SL(%f,%f,0,%f,%f,0){1,1};", coord[3*(n0-1)], coord[3*(n0-1)+1], coord[3*(n1-1)], coord[3*(n1-1)+1]);
      }
    } else {
      // P2 triangles - Add all possible P2 sub-edges
      for(size_t k = 0; k < nSubEdges; ++k) {
        size_t n0 = elemNodeTags[0][numNodesPerElem * i + e0[k]];
        size_t n1 = elemNodeTags[0][numNodesPerElem * i + e1[k]];
        edges.insert(std::make_pair(n0, n1));
      }
    }
  }

  ////////////////////////////////////////////////////
  // The circle test case from Alauzet
  // for(size_t i = 0; i < nodeTags.size(); ++i) {
  //   double x[2] = {coord[3*i], coord[3*i+1]};
  //   double theta = atan2(x[1], x[0]);
  //   double v1[2] = {cos(theta),sin(theta)};
  //   double v2[2] = {-sin(theta),cos(theta)};
  //   double h1 = 0.0005 + 1.5 * fabs(1. - sqrt(x[0]*x[0] + x[1]*x[1]));
  //   double h2 = 0.1 * sqrt(x[0]*x[0] + x[1]*x[1]) + 1.5 * fabs(1. - sqrt(x[0]*x[0] + x[1]*x[1]));
  //   double evalues[2] = {1./(h1*h1), 1./(h2*h2)};
  //   MetricTensor metric(evalues, v1, v2);
  //   metricsAtNodeTags.at(nodeTags[i]).assignMatrixFrom(metric);
  // }
  ////////////////////////////////////////////////////

  fprintf(myFile, "};"); fclose(myFile);

  int iter = 0, numCorrected;
  bool correction = true;
  
  if(_options.gradationSpace != GradationSpace::ExpMetric) {
    while(correction && iter < _options.gradationMaxIter)
    {
    	numCorrected = 0; correction = false; iter++;
      for(auto edge : edges) {
      	if(gradationOnEdge(_options.gradationSpace, gradation, edge, coord, metricsAtNodeTags)) {
          numCorrected++;
      	  correction = true;
        };
      }
      feInfoCond(FE_VERBOSE > 0, "Gradation: Passe %d - Corrected %d edges", iter, numCorrected);
    }
  } else {

    double alpha = 1.05;
    double beta = gradation;

    std::map<int, MetricTensor> mNew;

    while(correction && iter < _options.gradationMaxIter)
    {
      numCorrected = 0; correction = false; iter++;
      
      for(size_t i = 0; i < nodeTags.size(); ++i) {
        mNew[nodeTags[i]].assignMatrixFrom(metricsAtNodeTags.at(nodeTags[i]));
      }

      for(auto edge : edges) {
        if(gradationOnEdgeExpMetric(_options.gradationSpace, beta, edge, coord, metricsAtNodeTags, mNew)){
          numCorrected++;
          correction = true;
        }
      }

      if(correction) beta *= alpha;

      for(size_t i = 0; i < nodeTags.size(); ++i) {
        metricsAtNodeTags.at(nodeTags[i]).assignMatrixFrom(mNew.at(nodeTags[i]));
      }

      feInfoCond(FE_VERBOSE > 0, "Gradation: Passe %d - Corrected %d edges", iter, numCorrected);
    }

    mNew.clear();
  }

  return FE_STATUS_OK;
#else
  return feErrorMsg(FE_STATUS_ERROR, "Gmsh is required to apply gradation to the metric tensor field!");
#endif
}