
#include "feMetric.h"
#include "feMatrixInterface.h"

#if defined(HAVE_GMSH)
#include "gmsh.h"
#endif

extern int FE_VERBOSE;

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

MetricTensor spanMetric(const GradationSpace space, const double gradation, 
	const MetricTensor &M, const double pq[2])
{
	if(space == GradationSpace::Metric) {
		// Homogeneous growth in the metric space (more anisotropic)
	  double eta = (1. + sqrt(M.dotProduct(pq, pq)) * log(gradation));
	  eta = 1. / (eta * eta);
	  return M * eta;
	} else {
  	// Homogeneous growth in the physical space (becomes isotropic far from p)
		return M.spanMetricInPhysicalSpace(gradation, pq);
  }
}

double matNorm2(const MetricTensor &m)
{
  double sqr = 0;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      sqr += m(i, j) * m(i, j);
    }
  }
  return sqrt(sqr);
}

double matNorm2Diff(const MetricTensor &m1, const MetricTensor &m2)
{
  double sqr = 0;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
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
thread_local double EV[2], V0[2], V1[0]; // eigenvalues

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

  double relDiff = matNorm2Diff(m1, m2) / fmin(matNorm2(m1), matNorm2(m2));

  bool matricesAreEqual = relDiff < REL_TOL_EQUAL;

  // If metrics are identical, return one of them
  if(matricesAreEqual) {
    return m1.copy();
  }

  bool m1IsDiagonal = fabs(b1)/matNorm2(m1) < REL_TOL_DIAGONAL;
  bool m2IsDiagonal = fabs(b2)/matNorm2(m2) < REL_TOL_DIAGONAL;

  if(m1IsDiagonal && m2IsDiagonal) {
    // Both matrices are diagonal
    EV[0] = fmax(a1, a2);
    EV[1] = fmax(c1, c2);
    return MetricTensor(EV, EX, EY);

  } else if(m1IsDiagonal) {
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
    v00 = (a1 * c2) / (a2 * b1) - (a1 * c2 + a2 * c1 +
                                   sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 +
                                        a2 * a2 * c1 * c1 + 4 * a2 * b1 * b1 * c2)) /
                                    (2 * a2 * b1);
    v01 = (a1 * c2) / (a2 * b1) - (a1 * c2 + a2 * c1 -
                                   sqrt(a1 * a1 * c2 * c2 - 2 * a1 * a2 * c1 * c2 +
                                        a2 * a2 * c1 * c1 + 4 * a2 * b1 * b1 * c2)) /
                                    (2 * a2 * b1);

    V0[0] = v00;
    V0[1] = 1.;
    V1[0] = v01;
    V1[1] = 1.;
    double l0 = fmax(m1.dotProduct(V0, V0), m2.dotProduct(V0, V0));
    double l1 = fmax(m1.dotProduct(V1, V1), m2.dotProduct(V1, V1));

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
    // Check if metrics are multiple of one another
    if(fabs(a2 / a1 - c2 / c1) < 1e-3 && fabs(a2 / a1 - b2 / b1) < 1e-3) {
      return (a2 / a1 <= 1) ? m1.copy() : m2.copy();
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

//
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

  TRACER ICI LES METRIQUES INITIALES ET LES SPANNED

  // Span Mp to q, intersect and check if Mq needs to be reduced
  MetricTensor MpAtq = spanMetric(space, gradation, Mp, pq);
  MpAtq = intersectionReductionSimultanee(Mq, MpAtq);
  if(matNorm2Diff(Mq, MpAtq) / matNorm2(Mq) > REL_TOL_GRADATION) {
    Mq.assignMatrixFrom(MpAtq);
    metricChanged = true;
  }

  // Span Mq to p
  MetricTensor MqAtp = spanMetric(space, gradation, Mq, qp);
  MqAtp = intersectionReductionSimultanee(Mp, MqAtp);

  if(matNorm2Diff(Mp, MqAtp) / matNorm2(Mp) > REL_TOL_GRADATION) {
    Mp.assignMatrixFrom(MqAtp);
    metricChanged = true;
  }

  return metricChanged;
}

// Apply gradation to the metric tensor field.
// Gradation is performed on the Gmsh model, edge by edge according to Alauzet's paper.
feStatus feMetric::newGradation(std::vector<std::size_t> &nodeTagsFoo, std::vector<double> &coord,
	std::map<int, MetricTensor> &metricsAtNodeTags)
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

  for(size_t i = 0; i < elementTags[0].size(); i++) {
    for(size_t j = 0; j < 3; ++j) {

      if(numNodesPerElem == 3) {
        // P1 triangles
        // Add each edge
        size_t n0 = elemNodeTags[0][numNodesPerElem * i + j];
        size_t n1 = elemNodeTags[0][numNodesPerElem * i + (j+1) % numNodesPerElem];
        edges.insert(std::make_pair(n0, n1));
      } else {
        // P2 triangles
        // Add the 6 semi-edges, from a P1 vertex to a P2 midnode
        size_t n0 = elemNodeTags[0][numNodesPerElem * i + j];
        size_t n1 = elemNodeTags[0][numNodesPerElem * i + j+3];
        size_t n2 = elemNodeTags[0][numNodesPerElem * i + (j+1) % 3];
        edges.insert(std::make_pair(n0, n1));
        edges.insert(std::make_pair(n1, n2));

        fprintf(myFile, "SL(%.16g,%.16g,0.,%.16g,%.16g,0.){%f, %f};\n",
        	coord[3 * (n0-1)],
        	coord[3 * (n0-1) + 1],
        	coord[3 * (n1-1)],
        	coord[3 * (n1-1) + 1], 1., 1.);
       	fprintf(myFile, "SL(%.16g,%.16g,0.,%.16g,%.16g,0.){%f, %f};\n",
        	coord[3 * (n1-1)],
        	coord[3 * (n1-1) + 1],
        	coord[3 * (n2-1)],
        	coord[3 * (n2-1) + 1], 1., 1.);
      }
    }
  }

  fprintf(myFile, "};");
  fclose(myFile);

  int iter = 0, numCorrected;
  bool correction = true, edgeWasCorrected;
  while(correction && iter < _options.gradationMaxIter)
  {
  	numCorrected = 0;
    correction = false;
    iter++;

    for(auto edge : edges) {
    	edgeWasCorrected = gradationOnEdge(_options.gradationSpace, _options.gradation, edge, coord, metricsAtNodeTags);
    	if(edgeWasCorrected) numCorrected++;
    	// correction |= edgeWasCorrected;
    }

    feInfoCond(FE_VERBOSE > 0, "Gradation: Passe %d - Corrected %d edges", iter, numCorrected);
  }

  return FE_STATUS_OK;
#else
  return feErrorMsg(FE_STATUS_ERROR, "Gmsh is required to apply gradation to the metric tensor field!");
#endif
}