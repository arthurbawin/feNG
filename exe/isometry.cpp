
#include "feAPI.h"
#include "geodesic.h"

#if defined(HAVE_CERES)
#include "/usr/local/include/ceres/ceres.h"
#endif

using ceres::CostFunction;
using ceres::Problem;
using ceres::SizedCostFunction;
using ceres::Solve;
using ceres::Solver;

/* 10 point rule points */
static int nQuad = 10;
static double xQuad[10] = {-9.739065285171716e-01, -8.650633666889845e-01,
                           -6.794095682990244e-01, -4.333953941292472e-01,
                           -1.488743389816312e-01, 1.488743389816312e-01,
                           4.333953941292472e-01,  6.794095682990244e-01,
                           8.650633666889845e-01,  9.739065285171716e-01};

/* 10 point rule weights */
static double wQuad[10] = {6.667134430868774e-02, 1.494513491505805e-01,
                           2.190863625159822e-01, 2.692667193099962e-01,
                           2.955242247147529e-01, 2.955242247147529e-01,
                           2.692667193099962e-01, 2.190863625159822e-01,
                           1.494513491505805e-01, 6.667134430868774e-02};


double TARGET_LENGTH_IN_METRIC_SUBTRIANGLE;
double TARGET_AREA_IN_METRIC_SUBTRIANGLE;

double A0 = sqrt(3.)/4.;

bool SET_RATIO = true;
double EDGE2TRI_RATIO = 1.;

bool PRINT_INFO = false;

MetricTensor METRIC, DMDX, DMDY, DMDX_DUMMY, DMDY_DUMMY;
double METRIC_ARRAY[2][2], DMDX_ARRAY[2][2], DMDY_ARRAY[2][2], DMDX_DUMMY_ARRAY[2][2], DMDY_DUMMY_ARRAY[2][2];

double determinant(const double m[2][2])
{
  return m[0][0]*m[1][1] - m[0][1]*m[1][0];
}

double dotProduct(const double m[2][2], const double x1[2], const double x2[2])
{
  return x2[0] * (m[0][0] * x1[0] + m[0][1] * x1[1])
       + x2[1] * (m[0][1] * x1[0] + m[1][1] * x1[1]);
}

void getMetric(MetricTensor &m, MetricTensor &dmdx, MetricTensor &dmdy, double pos[2])
{
  // m(0,0) = 2.0;
  // m(0,1) = 1.0;
  // m(1,0) = 1.0;
  // m(1,1) = 3.0;

  // dmdx(0,0) = 0.;
  // dmdx(0,1) = 0.;
  // dmdx(1,0) = 0.;
  // dmdx(1,1) = 0.;

  // dmdy(0,0) = 0.;
  // dmdy(0,1) = 0.;
  // dmdy(1,0) = 0.;
  // dmdy(1,1) = 0.;

  double x = pos[0];
  double y = pos[1];

  m(0,0) = 3./( (1.+x)*(1.+x) );
  m(0,1) = 0.; //1./7.;
  m(1,0) = 0.; //1./7.;
  m(1,1) = 3./( (1.+y)*(1.+y) );

  // if(m.determinant() <= 0) {
  //   feErrorMsg(FE_STATUS_ERROR, "Determinant négatif");
  //   exit(-1);
  // }

  dmdx(0,0) = 3. * (-2./((1+x)*(1+x)*(1+x)));
  dmdx(0,1) = 0.0;
  dmdx(1,0) = 0.0;
  dmdx(1,1) = 0.0;

  dmdy(0,0) = 0.0;
  dmdy(0,1) = 0.0;
  dmdy(1,0) = 0.0;
  dmdy(1,1) = 3. * (-2./((1+y)*(1+y)*(1+y)));

  // Iso et grad
  // double l1 = 1/(0.5*0.5); 
  // double l2 = 1/(0.1*0.1);

  // m(0,0) = 1./(x*x+y*y) * (l2*x*x + l1*y*y);
  // m(0,1) = 1./(x*x+y*y) * (l2-l1) * x*y;
  // m(1,0) = 1./(x*x+y*y) * (l2-l1) * x*y;
  // m(1,1) = 1./(x*x+y*y) * (l1*x*x + l2*y*y);

  // dmdx(0,0) = (l1-l2)/((x*x+y*y)*(x*x+y*y)) * -2.*x*y*y;
  // dmdx(0,1) = (l1-l2)/((x*x+y*y)*(x*x+y*y)) * y*(x*x-y*y);
  // dmdx(1,0) = (l1-l2)/((x*x+y*y)*(x*x+y*y)) * y*(x*x-y*y);
  // dmdx(1,1) = (l1-l2)/((x*x+y*y)*(x*x+y*y)) *  2.*x*y*y;

  // dmdy(0,0) = (l1-l2)/((x*x+y*y)*(x*x+y*y)) *  2.*x*x*y;
  // dmdy(0,1) = (l1-l2)/((x*x+y*y)*(x*x+y*y)) * -x*(x*x-y*y);
  // dmdy(1,0) = (l1-l2)/((x*x+y*y)*(x*x+y*y)) * -x*(x*x-y*y);
  // dmdy(1,1) = (l1-l2)/((x*x+y*y)*(x*x+y*y)) * -2.*x*x*y;
}

void getMetric(double m[2][2], double dmdx[2][2], double dmdy[2][2], double pos[2])
{
  // m(0,0) = 2.0;
  // m(0,1) = 1.0;
  // m(1,0) = 1.0;
  // m(1,1) = 3.0;

  // dmdx(0,0) = 0.;
  // dmdx(0,1) = 0.;
  // dmdx(1,0) = 0.;
  // dmdx(1,1) = 0.;

  // dmdy(0,0) = 0.;
  // dmdy(0,1) = 0.;
  // dmdy(1,0) = 0.;
  // dmdy(1,1) = 0.;

  double x = pos[0];
  double y = pos[1];

  m[0][0] = 3./( (1.+x)*(1.+x) );
  m[0][1] = 0.; //1./7.;
  m[1][0] = 0.; //1./7.;
  m[1][1] = 3./( (1.+y)*(1.+y) );

  // if(m.determinant() <= 0) {
  //   feErrorMsg(FE_STATUS_ERROR, "Determinant négatif");
  //   exit(-1);
  // }

  dmdx[0][0] = 3. * (-2./((1+x)*(1+x)*(1+x)));
  dmdx[0][1] = 0.0;
  dmdx[1][0] = 0.0;
  dmdx[1][1] = 0.0;

  dmdy[0][0] = 0.0;
  dmdy[0][1] = 0.0;
  dmdy[1][0] = 0.0;
  dmdy[1][1] = 3. * (-2./((1+y)*(1+y)*(1+y)));

  // Iso et grad
  // double l1 = 1/(0.5*0.5); 
  // double l2 = 1/(0.1*0.1);

  // m(0,0) = 1./(x*x+y*y) * (l2*x*x + l1*y*y);
  // m(0,1) = 1./(x*x+y*y) * (l2-l1) * x*y;
  // m(1,0) = 1./(x*x+y*y) * (l2-l1) * x*y;
  // m(1,1) = 1./(x*x+y*y) * (l1*x*x + l2*y*y);

  // dmdx(0,0) = (l1-l2)/((x*x+y*y)*(x*x+y*y)) * -2.*x*y*y;
  // dmdx(0,1) = (l1-l2)/((x*x+y*y)*(x*x+y*y)) * y*(x*x-y*y);
  // dmdx(1,0) = (l1-l2)/((x*x+y*y)*(x*x+y*y)) * y*(x*x-y*y);
  // dmdx(1,1) = (l1-l2)/((x*x+y*y)*(x*x+y*y)) *  2.*x*y*y;

  // dmdy(0,0) = (l1-l2)/((x*x+y*y)*(x*x+y*y)) *  2.*x*x*y;
  // dmdy(0,1) = (l1-l2)/((x*x+y*y)*(x*x+y*y)) * -x*(x*x-y*y);
  // dmdy(1,0) = (l1-l2)/((x*x+y*y)*(x*x+y*y)) * -x*(x*x-y*y);
  // dmdy(1,1) = (l1-l2)/((x*x+y*y)*(x*x+y*y)) * -2.*x*x*y;
}

// If the derivatives are not needed
void getMetric(MetricTensor &m, double pos[2])
{
  getMetric(m, DMDX_DUMMY, DMDY_DUMMY, pos);
}

// If the derivatives are not needed
void getMetric(double m[2][2], double pos[2])
{
  getMetric(m, DMDX_DUMMY_ARRAY, DMDY_DUMMY_ARRAY, pos);
}

double getEdgeLengthInMetricSquared(MetricTensor &metricAtMidEdge, double x0[2], double x1[2])
{
  double d[2] = {x1[0] - x0[0], x1[1] - x0[1]};
  return metricAtMidEdge.dotProduct(d,d);
}

double getEdgeLengthInMetricSquared(double metricAtMidEdge[2][2], double x0[2], double x1[2])
{
  double d[2] = {x1[0] - x0[0], x1[1] - x0[1]};
  return dotProduct(metricAtMidEdge,d,d);
}

void getEdgeLengthInMetricSquared_gradient_x0(const MetricTensor &metricAtMidEdge,
  const MetricTensor &dMdx, const MetricTensor &dMdy, const double x0[2], const double x1[2],
  double grad[2])
{
  double d[2] = {x1[0] - x0[0], x1[1] - x0[1]};
  double dd_dx0[2] = {-1., 0.};
  double dd_dy0[2] = {0., -1.};

  double dlSquared_dx0 = metricAtMidEdge.dotProduct(dd_dx0,d) + 0.5 * dMdx.dotProduct(d,d) + metricAtMidEdge.dotProduct(d,dd_dx0);
  double dlSquared_dy0 = metricAtMidEdge.dotProduct(dd_dy0,d) + 0.5 * dMdy.dotProduct(d,d) + metricAtMidEdge.dotProduct(d,dd_dy0);

  grad[0] = dlSquared_dx0;
  grad[1] = dlSquared_dy0;
}

void getEdgeLengthInMetricSquared_gradient_x0(const double metricAtMidEdge[2][2],
  const double dMdx[2][2], const double dMdy[2][2], const double x0[2], const double x1[2], double grad[2])
{
  double d[2] = {x1[0] - x0[0], x1[1] - x0[1]};
  double dd_dx0[2] = {-1., 0.};
  double dd_dy0[2] = {0., -1.};

  double dlSquared_dx0 = dotProduct(metricAtMidEdge,dd_dx0,d) + 0.5 * dotProduct(dMdx,d,d) + dotProduct(metricAtMidEdge,d,dd_dx0);
  double dlSquared_dy0 = dotProduct(metricAtMidEdge,dd_dy0,d) + 0.5 * dotProduct(dMdy,d,d) + dotProduct(metricAtMidEdge,d,dd_dy0);

  grad[0] = dlSquared_dx0;
  grad[1] = dlSquared_dy0;
}

double computeTriJacobianDet(const std::vector<Vertex*> &tri)
{
  return (tri[1]->x() - tri[0]->x()) * (tri[2]->y() - tri[0]->y()) - (tri[2]->x() - tri[0]->x()) * (tri[1]->y() - tri[0]->y());
}

double computeAreaTri(const std::vector<Vertex*> &tri)
{
  return 0.5 * computeTriJacobianDet(tri);
}

double computeCostFunction(const std::vector<Vertex*> &vertices,
  const std::map<Vertex*, std::set<Edge, EdgeLessThan>> &vertNeighbours,
  const std::vector<std::vector<Vertex*>> &elements)
{
  double costEdges = 0.;
  for(size_t i = 0; i < vertices.size(); ++i) {
    Vertex *v = vertices[i];
    for(auto e : vertNeighbours.at(v)) {
      double xi = e.getVertex(0)->x();
      double yi = e.getVertex(0)->y();
      double xj = e.getVertex(1)->x();
      double yj = e.getVertex(1)->y();
      double x0[2] = {xi, yi};
      double x1[2] = {xj, yj};

      // This sets the Euclidean length to the target. When the metric is not diagonal constant,
      // it is not always possible to write the length in the metric as a function of the Euclidean length,
      // so we do the opposite and impose the length in the metric
      // costEdges += 0.5 * pow(lengthSquared - TARGET_LENGTH_SUBTRIANGLE*TARGET_LENGTH_SUBTRIANGLE, 2);

      // TARGET_LENGTH_IN_METRIC_SUBTRIANGLE is always 1/N_SUBTRIANGLES
      double xMid[2] = {(xi+xj)/2., (yi+yj)/2.};
      getMetric(METRIC_ARRAY, DMDX_ARRAY, DMDY_ARRAY, xMid);
      double lengthInMetricSquared = getEdgeLengthInMetricSquared(METRIC_ARRAY, x0, x1);
      costEdges += 0.5 * pow(lengthInMetricSquared - TARGET_LENGTH_IN_METRIC_SUBTRIANGLE*TARGET_LENGTH_IN_METRIC_SUBTRIANGLE, 2);
    }
  }

  double costTri = 0.;
  for(size_t i = 0; i < elements.size(); ++i) {
    double euclideanArea = computeAreaTri(elements[i]);
    Vertex *v0 = elements[i][0];
    Vertex *v1 = elements[i][1];
    Vertex *v2 = elements[i][2];
    double xBary[2] = {(v0->x() + v1->x() + v2->x())/3., (v0->y() + v1->y() + v2->y())/3.};
    getMetric(METRIC_ARRAY, DMDX_ARRAY, DMDY_ARRAY, xBary);
    double areaInMetric = euclideanArea * sqrt(determinant(METRIC_ARRAY));

    double X = areaInMetric/TARGET_AREA_IN_METRIC_SUBTRIANGLE;
    costTri += log(X) * log(X);
  }

  if(SET_RATIO) {
    SET_RATIO = false;
    EDGE2TRI_RATIO = costEdges/costTri;
  }

  return costEdges + EDGE2TRI_RATIO * costTri;
}

// Compute [dE/dxi, dE/dyi]
void computeCostGradient(int iVertex, const std::vector<Vertex*> &verticesToModify,
  const std::map<Vertex*, std::set<Edge, EdgeLessThan>> &vertNeighbours,
  const std::vector<std::vector<Vertex*>> &elements,
  const std::map<Vertex*, std::set<size_t>> &vert2tri,
  double grad[2])
{
  grad[0] = 0.;
  grad[1] = 0.;

  Vertex *v = verticesToModify[iVertex];

  for(auto e : vertNeighbours.at(v)) {
    double xi = e.getVertex(0)->x(); // xi est toujours le vertex courant (voir construction des aretes)
    double yi = e.getVertex(0)->y(); // et xj est l'autre extremite
    double xj = e.getVertex(1)->x();
    double yj = e.getVertex(1)->y();
    double x0[2] = {xi, yi};
    double x1[2] = {xj, yj};

    double xMid[2] = {(xi+xj)/2., (yi+yj)/2.};
    getMetric(METRIC_ARRAY, DMDX_ARRAY, DMDY_ARRAY, xMid);
    double lengthInMetricSquared = getEdgeLengthInMetricSquared(METRIC_ARRAY, x0, x1);
    double gradLSquared[2];
    getEdgeLengthInMetricSquared_gradient_x0(METRIC_ARRAY, DMDX_ARRAY, DMDY_ARRAY, x0, x1, gradLSquared);

    // Double the contribution because each edge is counted twice in the cost,
    // this cancels the 1/2 of the cost function
    grad[0] += 2. * (lengthInMetricSquared + (-1.) * TARGET_LENGTH_IN_METRIC_SUBTRIANGLE*TARGET_LENGTH_IN_METRIC_SUBTRIANGLE) * gradLSquared[0];
    grad[1] += 2. * (lengthInMetricSquared + (-1.) * TARGET_LENGTH_IN_METRIC_SUBTRIANGLE*TARGET_LENGTH_IN_METRIC_SUBTRIANGLE) * gradLSquared[1];
  }

  for(auto iTri : vert2tri.at(v)) {
    double ddetJdx, ddetJdy;

    double x1 = elements[iTri][0]->x();
    double y1 = elements[iTri][0]->y();
    double x2 = elements[iTri][1]->x();
    double y2 = elements[iTri][1]->y();
    double x3 = elements[iTri][2]->x();
    double y3 = elements[iTri][2]->y();

    if(v == elements[iTri][0]) {
      ddetJdx = -(y3-y1) + (y2-y1);
      ddetJdy = -(x2-x1) + (x3-x1);
    } else if(v == elements[iTri][1]) {
      ddetJdx =  (y3-y1);
      ddetJdy = -(x3-x1);
    } else if(v == elements[iTri][2]) {
      ddetJdx = -(y2-y1);
      ddetJdy =  (x2-x1);
    }

    double d_areaEuclidian_x0 = 0.5 * ddetJdx;
    double d_areaEuclidian_y0 = 0.5 * ddetJdy;

    double euclideanArea = computeAreaTri(elements[iTri]);

    Vertex *v0 = elements[iTri][0];
    Vertex *v1 = elements[iTri][1];
    Vertex *v2 = elements[iTri][2];
    double xBary[2] = {(v0->x() + v1->x() + v2->x())/3., (v0->y() + v1->y() + v2->y())/3.};
    getMetric(METRIC_ARRAY, DMDX_ARRAY, DMDY_ARRAY, xBary);
    double sq_det = sqrt(determinant(METRIC_ARRAY));
    double areaInMetric = euclideanArea * sq_det;
    double X = areaInMetric/TARGET_AREA_IN_METRIC_SUBTRIANGLE;

    // Use Jacobi's identity to compute derivative of determinant
    // double dDetdx0 = 1./3. * (METRIC(1,1)*DMDX(0,0) - 2.*METRIC(0,1)*DMDX(0,1) + METRIC(0,0)*DMDX(1,1));
    // double dDetdy0 = 1./3. * (METRIC(1,1)*DMDY(0,0) - 2.*METRIC(0,1)*DMDY(0,1) + METRIC(0,0)*DMDY(1,1));
    double dDetdx0 = 1./3. * (METRIC_ARRAY[1][1]*DMDX_ARRAY[0][0] - 2.*METRIC_ARRAY[0][1]*DMDX_ARRAY[0][1] + METRIC_ARRAY[0][0]*DMDX_ARRAY[1][1]);
    double dDetdy0 = 1./3. * (METRIC_ARRAY[1][1]*DMDY_ARRAY[0][0] - 2.*METRIC_ARRAY[0][1]*DMDY_ARRAY[0][1] + METRIC_ARRAY[0][0]*DMDY_ARRAY[1][1]);

    double d_areaMetric_dx0 = d_areaEuclidian_x0 * sq_det + euclideanArea / (2.*sq_det) * dDetdx0;
    double d_areaMetric_dy0 = d_areaEuclidian_y0 * sq_det + euclideanArea / (2.*sq_det) * dDetdy0;
    double logx = log(X);
    grad[0] += EDGE2TRI_RATIO * 2. * logx / X * d_areaMetric_dx0 / TARGET_AREA_IN_METRIC_SUBTRIANGLE;
    grad[1] += EDGE2TRI_RATIO * 2. * logx / X * d_areaMetric_dy0 / TARGET_AREA_IN_METRIC_SUBTRIANGLE;
  }
}

// Cost function with analytic gradient
class costFunctor final : public ceres::FirstOrderFunction {
  public:
    int _numParameters;
    int _numVertices;
    int _numVerticesToModify;
    std::vector<Vertex*> _vertices;
    std::vector<Vertex*> _verticesToModify;
    std::map<Vertex*, std::set<Edge, EdgeLessThan>> _vertNeighbours;
    std::vector<std::vector<Vertex*>> _elements;
    std::map<Vertex*, std::set<size_t>> _vert2tri;
    mutable std::vector<double> INITIALPOS;

  public:
    costFunctor(std::vector<Vertex*> &vertices,
      std::vector<Vertex*> &verticesToModify,
      std::map<Vertex*, std::set<Edge, EdgeLessThan>> &vertNeighbours,
      std::vector<std::vector<Vertex*>> &elements,
      std::map<Vertex*, std::set<size_t>> &vert2tri,
      int numVertices,
      int numVerticesToModify,
      int numParameters)
      : _vertices(vertices),
      _verticesToModify(verticesToModify),
      _vertNeighbours(vertNeighbours),
      _elements(elements),
      _vert2tri(vert2tri),
      _numVertices(numVertices),
      _numVerticesToModify(numVerticesToModify),
      _numParameters(numParameters)
      {
        INITIALPOS.resize(numParameters, 0.);
      };

    // The objective (cost) function and its gradient
    bool Evaluate(const double* parameters,
                  double* cost,
                  double* gradient) const override
    {
      // Save initial positions
      int cnt = 0;
      for(int i = 0; i < _numVerticesToModify; ++i) {
        Vertex *v = _verticesToModify[i];
        INITIALPOS[cnt] = v->x(); cnt++;
        INITIALPOS[cnt] = v->y(); cnt++;
      }

      // Update the positions
      cnt = 0;
      for(int i = 0; i < _numVerticesToModify; ++i) {
        Vertex *v = _verticesToModify[i];
        (*v)(0) = parameters[cnt]; cnt++;
        (*v)(1) = parameters[cnt]; cnt++;
      }

      // COMPUTE OBJECTIVE FUNCTION
      cost[0] = computeCostFunction(_vertices, _vertNeighbours, _elements);

      // COMPUTE GRADIENT
      if (gradient) {

        for(int i = 0; i < _numParameters; ++i)
          gradient[i] = 0.;

        cnt = 0;
        double vertexGradient[2] = {0., 0.};
        for(int i = 0; i < _numVerticesToModify; ++i) {
          computeCostGradient(i, _verticesToModify, _vertNeighbours, _elements, _vert2tri, vertexGradient);
          gradient[cnt] = vertexGradient[0]; cnt++;
          gradient[cnt] = vertexGradient[1]; cnt++;
        }
      }

      // Restore initial positions
      cnt = 0;
      for(int i = 0; i < _numVerticesToModify; ++i) {
        Vertex *v = _verticesToModify[i];
        (*v)(0) = INITIALPOS[cnt]; cnt++;
        (*v)(1) = INITIALPOS[cnt]; cnt++;
      }

      return true;
    }

    // When using numeric differentiation: only compute cost
    bool operator()(const double* parameters, double* cost) const 
    {
      return this->Evaluate(parameters, cost, nullptr);
    }

    int NumParameters() const override { return _numParameters; }
};

bool isInsideTriangle(std::vector<Vertex*> &tri_abc, Vertex *p) {
  std::vector<Vertex*> abp{ tri_abc[0], tri_abc[1], p};
  std::vector<Vertex*> bcp{ tri_abc[1], tri_abc[2], p};
  std::vector<Vertex*> cap{ tri_abc[2], tri_abc[0], p};
  double a1 = computeAreaTri(abp);
  double a2 = computeAreaTri(bcp);
  double a3 = computeAreaTri(cap);
  if (a1 >= 0 && a2 >= 0 && a3 >= 0) {
      return true;
  } else {
      return false;
  }
}

void getMetricFixedSubTriangle_fixedPoint(const double x0[2], const double R[2][2], std::vector<Vertex*> tri, MetricTensor &M, const double tol)
{
  FILE *myFile = fopen("fixedpoint.pos", "w"); fprintf(myFile, "View \" fixedpoint \"{\n");

  double xBary[2] = {x0[0], x0[1]}, xBaryPrev[2] = {x0[0] + 1., x0[1] + 1.};
  double dist;
  fprintf(myFile, "SP(%g,%g,0.){%g};\n", xBary[0], xBary[1], 1.);

  do {
    getMetric(M, xBary);
    MetricTensor Mmud = M.pow(-1./2.);
    // Compute M^(-1/2) * R
    double prod[2][2] = {{Mmud(0,0)*R[0][0] + Mmud(0,1)*R[1][0], Mmud(0,0)*R[0][1] + Mmud(0,1)*R[1][1]},
                         {Mmud(1,0)*R[0][0] + Mmud(1,1)*R[1][0], Mmud(1,0)*R[0][1] + Mmud(1,1)*R[1][1]}};

    double v0[2] = {tri[0]->x(), tri[1]->y()}, v1[2], v2[2];

    // Move tri according to M^(-1/2) and rotation R
    v1[0] = x0[0] + prod[0][0] * (tri[1]->x() - v0[0]) + prod[0][1] * (tri[1]->y() - v0[1]);
    v1[1] = x0[1] + prod[1][0] * (tri[1]->x() - v0[0]) + prod[1][1] * (tri[1]->y() - v0[1]);
    v2[0] = x0[0] + prod[0][0] * (tri[2]->x() - v0[0]) + prod[0][1] * (tri[2]->y() - v0[1]);
    v2[1] = x0[1] + prod[1][0] * (tri[2]->x() - v0[0]) + prod[1][1] * (tri[2]->y() - v0[1]);
    // Recompute barycenter
    xBary[0] = (x0[0] + v1[0] + v2[0])/3.;
    xBary[1] = (x0[1] + v1[1] + v2[1])/3.;

    fprintf(myFile, "SP(%g,%g,0.){%g};\n", xBary[0], xBary[1], 1.);
    fprintf(myFile, "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g){%d,%d,%d};\n",
              x0[0], x0[1], 0.,
              v1[0], v1[1], 0.,
              v2[0], v2[1], 0., 1, 1, 1);

    // Compute distance between barycenter iterations
    dist = sqrt(pow(xBary[0] - xBaryPrev[0], 2) + pow(xBary[1] - xBaryPrev[1], 2));
    xBaryPrev[0] = xBary[0];
    xBaryPrev[1] = xBary[1];
    feInfo("dist = %+-1.4e", dist);
  } while(dist > tol);
  fprintf(myFile, "};"); fclose(myFile);
}

void createInitialMesh(const int N,
  const double x0_ref[2],
  const double theta,
  std::vector<Vertex> &vertices,
  std::vector<Vertex*> &verticesToModify,
  std::map<Vertex*, std::set<Edge, EdgeLessThan>> &vertNeighbours,
  std::vector<std::vector<size_t>> &boundaryVertices,
  std::vector<std::vector<Vertex*>> &elements,
  std::map<Vertex*, std::set<size_t>> &vert2tri,
  size_t &iTriBarycenter)
{
  FILE *myFile = fopen("edges.pos", "w");
  fprintf(myFile, "View \" edges \"{\n");

  // Discretize reference rectangle triangle then map it to an equilateral
  vertices.resize((N+1)*(N+2)/2);
  verticesToModify.clear();
  boundaryVertices.resize(3);
  std::vector<std::vector<Vertex*>> verticesMat(N+1);
  vert2tri.clear();

  int cnt = 0;
  for(int i = 0; i <= N; ++i) {
    verticesMat[i].resize(N-i+1);
    for(int j = 0; j <= N-i; ++j) {
      Vertex v( (double)j/(double)N, (double)i/(double)N, 0., cnt);
      vertices[cnt] = v;
      verticesMat[i][j] = &vertices[cnt];

      // Indices of boundary vertices
      if(i == 0) {
        boundaryVertices[0].push_back(cnt);
      } 
      if(j == N-i) {
        boundaryVertices[1].push_back(cnt);
      } 
      if(j == 0) {
        boundaryVertices[2].push_back(cnt);
      }
      cnt++;
    }
  }

  int cnt_edge = 0;
  for(int i = 0; i <= N; ++i) {
    for(int j = 0; j <= N-i; ++j) {
      // Right neighbour
      if(j-1 >= 0) {
        Edge new_edge(verticesMat[i][j], verticesMat[i][j-1], cnt_edge++);
        vertNeighbours[verticesMat[i][j]].insert(new_edge);
        vertNeighbours[verticesMat[i][j-1]].insert(new_edge.flip());
        fprintf(myFile, "SL(%g,%g,0.,%g,%g,0.){1.,1.};\n",
          new_edge.getVertex(0)->x(), new_edge.getVertex(0)->y(),
          new_edge.getVertex(1)->x(), new_edge.getVertex(1)->y());
      }

      // Left neighbour
      if(j+1 <= N-i) {
        Edge new_edge(verticesMat[i][j], verticesMat[i][j+1], cnt_edge++);
        vertNeighbours[verticesMat[i][j]].insert(new_edge);
        vertNeighbours[verticesMat[i][j+1]].insert(new_edge.flip());
        fprintf(myFile, "SL(%g,%g,0.,%g,%g,0.){1.,1.};\n",
          new_edge.getVertex(0)->x(), new_edge.getVertex(0)->y(),
          new_edge.getVertex(1)->x(), new_edge.getVertex(1)->y());
      }

      // Bottom neighbour
      if(i-1 >= 0) {
        Edge new_edge(verticesMat[i][j], verticesMat[i-1][j], cnt_edge++);
        vertNeighbours[verticesMat[i][j]].insert(new_edge);
        vertNeighbours[verticesMat[i-1][j]].insert(new_edge.flip());
        fprintf(myFile, "SL(%g,%g,0.,%g,%g,0.){1.,1.};\n",
          new_edge.getVertex(0)->x(), new_edge.getVertex(0)->y(),
          new_edge.getVertex(1)->x(), new_edge.getVertex(1)->y());
      }

      // Top neighbour
      if(i+1 <= N && j <= N-(i+1)) {
        Edge new_edge(verticesMat[i][j], verticesMat[i+1][j], cnt_edge++);
        vertNeighbours[verticesMat[i][j]].insert(new_edge);
        vertNeighbours[verticesMat[i+1][j]].insert(new_edge.flip());
        fprintf(myFile, "SL(%g,%g,0.,%g,%g,0.){1.,1.};\n",
          new_edge.getVertex(0)->x(), new_edge.getVertex(0)->y(),
          new_edge.getVertex(1)->x(), new_edge.getVertex(1)->y());
      }

      // Bottom-right neighbour
      if(i-1 >= 0) {
        Edge new_edge(verticesMat[i][j], verticesMat[i-1][j+1], cnt_edge++);
        vertNeighbours[verticesMat[i][j]].insert(new_edge);
        vertNeighbours[verticesMat[i-1][j+1]].insert(new_edge.flip());
        fprintf(myFile, "SL(%g,%g,0.,%g,%g,0.){1.,1.};\n",
          new_edge.getVertex(0)->x(), new_edge.getVertex(0)->y(),
          new_edge.getVertex(1)->x(), new_edge.getVertex(1)->y());
      }

      // Top-left neighbour
      if(i+1 <= N && j-1 >= 0) {
        Edge new_edge(verticesMat[i][j], verticesMat[i+1][j-1], cnt_edge++);
        vertNeighbours[verticesMat[i][j]].insert(new_edge);
        vertNeighbours[verticesMat[i+1][j-1]].insert(new_edge.flip());
        fprintf(myFile, "SL(%g,%g,0.,%g,%g,0.){1.,1.};\n",
          new_edge.getVertex(0)->x(), new_edge.getVertex(0)->y(),
          new_edge.getVertex(1)->x(), new_edge.getVertex(1)->y());
      }
    }
  }
  fprintf(myFile, "};"); fclose(myFile);

  // Create the triangles
  int cnt_tri = 0;
  FILE *triFile = fopen("initialMeshRef.pos", "w");
  fprintf(triFile, "View \" initialMeshRef \"{\n");
  for(int i = 0; i <= N; ++i) {
    for(int j = 0; j <= N-i; ++j) {
      
      // Has a left neighbour (thus has a top neighbour)
      if(j+1 <= N-i) { 
        std::vector<Vertex*> tri(3);
        tri[0] = verticesMat[i][j];
        tri[1] = verticesMat[i][j+1];
        tri[2] = verticesMat[i+1][j];
        elements.push_back(tri);
        vert2tri[tri[0]].insert(cnt_tri);
        vert2tri[tri[1]].insert(cnt_tri);
        vert2tri[tri[2]].insert(cnt_tri);
        cnt_tri++;
        fprintf(triFile, "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g){%d,%d,%d};\n",
            tri[0]->x(), tri[0]->y(), 0.,
            tri[1]->x(), tri[1]->y(), 0.,
            tri[2]->x(), tri[2]->y(), 0., 1, 1, 1);
      }

      // Has a top and a right neighbour
      if(i+1 <= N && j <= N-(i+1) && j-1 >= 0) { 
        std::vector<Vertex*> tri(3);
        tri[0] = verticesMat[i][j];
        tri[1] = verticesMat[i+1][j];
        tri[2] = verticesMat[i+1][j-1];
        elements.push_back(tri);
        vert2tri[tri[0]].insert(cnt_tri);
        vert2tri[tri[1]].insert(cnt_tri);
        vert2tri[tri[2]].insert(cnt_tri);
        cnt_tri++;
        fprintf(triFile, "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g){%d,%d,%d};\n",
             tri[0]->x(), tri[0]->y(), 0.,
             tri[1]->x(), tri[1]->y(), 0.,
             tri[2]->x(), tri[2]->y(), 0., 1, 1, 1);
      }
      cnt++;
    }
  }
  feInfo("There are %d triangles", elements.size());
  fprintf(triFile, "};"); fclose(triFile);

  // Find the triangle that contains the barycenter
  // Vertex barycenter((0. + 1. + 0.5)/3., (0. + 0. + sqrt(3.)/2.)/3., 0., -1);
  Vertex barycenter(1./3., 1./3., 0., -1);
  for(size_t i = 0; i < elements.size(); ++i) {
    if(isInsideTriangle(elements[i], &barycenter)) {
      iTriBarycenter = i;
      FILE *myFile = fopen("baryRef.pos", "w");
      fprintf(myFile, "View \" baryRef \"{\n");
      fprintf(myFile, "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g){%d,%d,%d};\n",
             elements[i][0]->x(), elements[i][0]->y(), 0.,
             elements[i][1]->x(), elements[i][1]->y(), 0.,
             elements[i][2]->x(), elements[i][2]->y(), 0., 1, 1, 1);
      fprintf(myFile, "};"); fclose(myFile);
      break;
    }
  }
  for(size_t i = 0; i < vertices.size(); ++i) {
    if(vert2tri.at(&vertices[i]).find(iTriBarycenter) == vert2tri.at(&vertices[i]).end()) {
      verticesToModify.push_back(&vertices[i]);
    }
  }

  // Fix the triangle that contains the barycenter and move the initial mesh
  // according to its transformation M^(-1/2) * R(theta).
  // Determine the metric at the barycenter in the physical space using a fixed point.
  double R[2][2] = {{cos(theta), -sin(theta)}, {sin(theta), cos(theta)}};
  getMetricFixedSubTriangle_fixedPoint(x0_ref, R, elements[iTriBarycenter], METRIC, 1e-6);

  // Compute M^(-1/2) * R * J0
  MetricTensor Mmud = METRIC.pow(-1./2.);
  double J0[2][2] = {{1., 1./2.},{0, sqrt(3.)/2.}};
  double MmudR[2][2] = {{Mmud(0,0)*R[0][0] + Mmud(0,1)*R[1][0], Mmud(0,0)*R[0][1] + Mmud(0,1)*R[1][1]},
                        {Mmud(1,0)*R[0][0] + Mmud(1,1)*R[1][0], Mmud(1,0)*R[0][1] + Mmud(1,1)*R[1][1]}};
  double J[2][2] = {{MmudR[0][0]*J0[0][0] + MmudR[0][1]*J0[1][0], MmudR[0][0]*J0[0][1] + MmudR[0][1]*J0[1][1]},
                    {MmudR[1][0]*J0[0][0] + MmudR[1][1]*J0[1][0], MmudR[1][0]*J0[0][1] + MmudR[1][1]*J0[1][1]}};                  

  for(auto &p : vertices) {
    double new_x = x0_ref[0] + J[0][0] * p(0) + J[0][1] * p(1);
    double new_y = x0_ref[1] + J[1][0] * p(0) + J[1][1] * p(1);
    p(0) = new_x;
    p(1) = new_y;
  }

  triFile = fopen("initialMesh.pos", "w");
  fprintf(triFile, "View \" initialMesh \"{\n");
  for(auto tri : elements) {
    fprintf(triFile, "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g){%d,%d,%d};\n",
        tri[0]->x(), tri[0]->y(), 0.,
        tri[1]->x(), tri[1]->y(), 0.,
        tri[2]->x(), tri[2]->y(), 0., 1, 1, 1);
  }
  fprintf(triFile, "};"); fclose(triFile);

  triFile = fopen("bary.pos", "w");
  fprintf(triFile, "View \" bary \"{\n");
  fprintf(triFile, "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g){%d,%d,%d};\n",
         elements[iTriBarycenter][0]->x(), elements[iTriBarycenter][0]->y(), 0.,
         elements[iTriBarycenter][1]->x(), elements[iTriBarycenter][1]->y(), 0.,
         elements[iTriBarycenter][2]->x(), elements[iTriBarycenter][2]->y(), 0., 1, 1, 1);
  fprintf(triFile, "};"); fclose(triFile);
}

double funBisection(double x0[2], double x1[2], double c)
{
  double xMid[2] = {(x0[0]+x1[0])/2., (x0[1]+x1[1])/2.};
  getMetric(METRIC, DMDX, DMDY, xMid);
  double xc[2] = {x0[0] + c * (x1[0] - x0[0]), x0[1] + c * (x1[1] - x0[1])};
  double lengthInMetricSquared = getEdgeLengthInMetricSquared(METRIC, x0, xc);
  return lengthInMetricSquared - TARGET_LENGTH_IN_METRIC_SUBTRIANGLE*TARGET_LENGTH_IN_METRIC_SUBTRIANGLE;
}

// Set c with a bisection such that x0 + c * (x1-x0) has length in metric = TARGET_LENGTH_IN_METRIC_SUBTRIANGLE
double relocateFixedEdge_bisection(double x0[2], double x1[2], double cMin, double cMax, int maxIter, double tol)
{
  double c, ca = cMin, cb = cMax;
  double fa = funBisection(x0,x1,ca);
  double fb = funBisection(x0,x1,cb);
  if(fa*fb > 0) {
    feErrorMsg(FE_STATUS_ERROR, "Bisection bracket does not work: f(%+-1.3e) = %+-1.3e and f(%+-1.3e) = %+-1.3e",
      ca, fa, cb, fb);
  }
  
  int iter = 0;
  while(iter++ <= maxIter) {
    c = (ca+cb)/2.;
    double fc = funBisection(x0,x1,c);
    // feInfo("Bisection for c: c = %+-1.3e - f(c) = %+-1.3e", c, fc);
    if(fabs(c) < tol || (cb-ca)/2. < tol) {
      return c;
    }
    if(fa*fc < 0) {
      cb = c;
      fb = fc;
    } else {
      ca = c;
      fa = fc;
    }
  }
  feErrorMsg(FE_STATUS_ERROR, "Bisection for c did not converge :/");
  return c;
}

int main(int argc, char** argv)
{
  double x0_ref[2] = {3., 1.};
  double theta = 0.; //M_PI/4.;

  #define NCONV 1
  double rate[NCONV], Etot[NCONV], h[NCONV];

  int N_SUBTRIANGLES = 10;

  for(int iConv = 0; iConv < NCONV; ++iConv, N_SUBTRIANGLES *= 2)
  {
    SET_RATIO = true;

    // Step size
    h[iConv] = 1./N_SUBTRIANGLES;

    TARGET_LENGTH_IN_METRIC_SUBTRIANGLE = 1. / N_SUBTRIANGLES;
    TARGET_AREA_IN_METRIC_SUBTRIANGLE = A0 / (N_SUBTRIANGLES*N_SUBTRIANGLES);

    std::vector<Vertex> vertices;
    std::vector<Vertex*> verticesToModify;
    std::map<Vertex*, std::set<Edge, EdgeLessThan>> vertNeighbours;
    std::vector<std::vector<Vertex*>> elements;
    std::map<Vertex*, std::set<size_t>> vert2tri;
    size_t iTriBarycenter;
    std::vector<std::vector<size_t>> boundaryVertices;

    // Create the initial mesh:
    // - discretize the reference rectangle triangle
    // - find the subtriangle that contains the barycenter
    // - compute the metric at the barycenter of its physical coordinates
    //   using a fixed point
    // - map it and the whole mesh using this constant metric to form the initial mesh
    createInitialMesh(N_SUBTRIANGLES, x0_ref, theta, vertices, verticesToModify, 
      vertNeighbours, boundaryVertices, elements, vert2tri, iTriBarycenter);

    size_t numVertices = vertices.size();
    size_t numVerticesToModify = verticesToModify.size();
    size_t numParameters = 2*numVerticesToModify;
    size_t numTri = elements.size();

    std::vector<Vertex*> verticesPtr(numVertices);
    for(int i = 0; i < numVertices; ++i) {
      verticesPtr[i] = &vertices[i];
    }

    std::vector<double> optimalPos(numParameters, 0.);

    // Choice 1 (bad, rotation is free): Fix the first vertex and send the whole mesh to the physical space
    // for(auto *v : verticesPtr) {
    //   (*v)(0) += x0[0];
    //   (*v)(1) += x0[1];
    //   feInfo("V pos = %f - %f", (*v)(0), (*v)(1));
    // }  

    // Choice 2: Fix the subtriangle containing the barycenter via a fixed vertex and a rotation
    // // Rotate around 0
    // for(auto *v : verticesPtr) {
    //   double vx = (*v)(0);
    //   double vy = (*v)(1);
    //   (*v)(0) = R[0][0] * vx + R[0][1] * vy;
    //   (*v)(1) = R[1][0] * vx + R[1][1] * vy;
    // }
    // double v0x = elements[iTriBarycenter][0]->x();
    // double v0y = elements[iTriBarycenter][0]->y();
    // // Send v0 of iTriBarycenter to x0_ref
    // for(auto *v : verticesPtr) {
    //   (*v)(0) += x0_ref[0] - v0x;
    //   (*v)(1) += x0_ref[1] - v0y;
    // }
    // // Move the other two vertices of the subtriangle containing the barycenter to have the required lengths
    // for(size_t i = 1; i < 3; ++i) {
    //   Vertex *v = elements[iTriBarycenter][i];
    //   double vx = (*v)(0);
    //   double vy = (*v)(1);

    //   // Set c with a bisection such that x0_ref + c * (v-x0_ref) has length in metric = TARGET_LENGTH_IN_METRIC_SUBTRIANGLE
    //   // double c = TARGET_LENGTH_SUBTRIANGLE / sqrt(( pow(vx - x0_ref[0],2) + pow(vy - x0_ref[1],2) ));
    //   double x1[2] = {vx, vy};
    //   double c = relocateFixedEdge_bisection(x0_ref, x1, 0., 10., 50, 1e-3);
    //   (*v)(0) = x0_ref[0] + c * (vx - x0_ref[0]);
    //   (*v)(1) = x0_ref[1] + c * (vy - x0_ref[1]);
    // }

    // // Untangle the linear mesh
    // double tangled;
    // do {
    //   tangled = false;
    //   for(size_t i = 0; i < elements.size(); ++i) {
    //     if(computeAreaTri(elements[i]) <= 0.)
    //       tangled = true;
    //   }
    //   if(tangled) {
    //     // Inflate the whole mesh from the barycenter
    //     double xBary[2] = {0., 0.};
    //     for(size_t i = 0; i < 3; ++i) {
    //       xBary[0] += elements[iTriBarycenter][i]->x();
    //       xBary[1] += elements[iTriBarycenter][i]->y();
    //     }
    //     xBary[0] /= 3.; xBary[1] /= 3.;

    //     for(auto *v : verticesToModify) {
    //       double vx = v->x();
    //       double vy = v->y();
    //       (*v)(0) = xBary[0] + 1.1 * (vx - xBary[0]);
    //       (*v)(1) = xBary[1] + 1.1 * (vy - xBary[1]);
    //     }
    //   }
    // } while(tangled);

    // // Move whole mesh by the jacobian of the fixed triangle
    // for(auto *v : verticesPtr) {
    //   double vx = (*v)(0);
    //   double vy = (*v)(1);
    //   (*v)(0) = R[0][0] * vx + R[0][1] * vy;
    //   (*v)(1) = R[1][0] * vx + R[1][1] * vy;
    // }

    // Initial positions
    int cnt = 0;
    for(int i = 0; i < numVerticesToModify; ++i) {
      Vertex *v = verticesToModify[i];
      optimalPos[cnt] = (*v)(0); cnt++;
      optimalPos[cnt] = (*v)(1); cnt++;
    }

    // //////////////////////////////////////////////////////
    // To test the gradient implementation with finite differences:
    // double hFD = 1e-8;
    // double f0 = computeCostFunction(verticesPtr, vertNeighbours, elements);
    // feInfo("Initial cost is %f", f0);
    // // Move vertices and check gradient
    // for(size_t i = 0; i < numVerticesToModify; ++i) {

    //   Vertex *v = verticesToModify[i];
    //   double x0 = (*v)(0);
    //   double y0 = (*v)(1);

    //   (*v)(0) += hFD;
    //   double dx = computeCostFunction(verticesPtr, vertNeighbours, elements);
    //   (*v)(0) = x0;

    //   (*v)(1) += hFD;
    //   double dy = computeCostFunction(verticesPtr, vertNeighbours, elements);
    //   (*v)(1) = y0;

    //   double grad[2] = {0., 0.};
    //   computeCostGradient(i, verticesToModify, vertNeighbours, elements, vert2tri, grad);

    //   double gx = (dx-f0)/hFD;
    //   double gy = (dy-f0)/hFD;
    //   double errorAbs = fmax(fabs(grad[0] - gx), fabs(grad[1] - gy));
    //   double errorRel = fmax(fabs(grad[0] - gx)/fmax(1e-12, fabs(gx)), fabs(grad[1] - gy)/fmax(1e-12, fabs(gy)));

    //   feInfo("FD    grad = %+-1.10e - %+-1.10e", gx, gy);
    //   feInfo("Exact grad = %+-1.10e - %+-1.10e - errorAbs = %1.6e - errorRel = %1.6e", grad[0], grad[1], errorAbs, errorRel);
    // }
    //   exit(-1);
    // //////////////////////////////////////////////////////

    if(true) {
      ceres::GradientProblemSolver::Summary summary;

      costFunctor *costFun = new costFunctor(verticesPtr, verticesToModify, vertNeighbours, elements, vert2tri, numVertices, numVerticesToModify, numParameters);
      
      // Analytic diff
      ceres::GradientProblem problem(costFun);

      // Numeric diff (only to test)
      // auto *foo = new ceres::NumericDiffFirstOrderFunction<costFunctor, ceres::CENTRAL>(costFun, numParameters);
      // ceres::GradientProblem problem(foo);

      ceres::GradientProblemSolver::Options options;
      // options.function_tolerance = 1e-12;
      // options.gradient_tolerance = 1e-10;
      // options.parameter_tolerance = 1e-12;
      // // options.min_line_search_step_size = 1e-11;
      options.max_num_iterations = 50000;
      // options.max_num_line_search_step_size_iterations = 1000;

      // options.line_search_direction_type = ceres::LBFGS;
      // options.line_search_type = ceres::WOLFE;
      // options.line_search_interpolation_type = ceres::CUBIC;
      // options.max_lbfgs_rank = 5000;

      options.minimizer_progress_to_stdout = true;

      ceres::Solve(options, problem, optimalPos.data(), &summary);

      std::cout << summary.FullReport() << "\n";

      if(summary.termination_type == ceres::TerminationType::CONVERGENCE || 
         summary.termination_type == ceres::TerminationType::NO_CONVERGENCE)
      {
        // Converged or ran out of iterations, but did not fail
        // Continue with the last computed candidate
        // return summary.final_cost;
      } else {
        // Failed
        feErrorMsg(FE_STATUS_ERROR, "Minimization did not converge");
      }
    }

    // Apply new position
    for(int i = 0; i < numVerticesToModify; ++i) {
      Vertex *v = verticesToModify[i];
      (*v)(0) = optimalPos[2*i];
      (*v)(1) = optimalPos[2*i+1];
    }

    // Print solution
    std::string name = "afterMoving_N" + std::to_string(N_SUBTRIANGLES) + ".pos";
    FILE *myFile = fopen(name.data(), "w");
    fprintf(myFile, "View \" %s \"{\n", name.data());
    for(int i = 0; i < numVertices; ++i) {
      Vertex *v = verticesPtr[i];
      fprintf(myFile, "SP(%g,%g,%g){%d};\n", (*v)(0), (*v)(1), (*v)(2), v->getTag());
      for(auto e : vertNeighbours[v]) {
      fprintf(myFile, "SL(%g,%g,0.,%g,%g,0.){%d,%d};\n",
                e.getVertex(0)->x(), e.getVertex(0)->y(),
                e.getVertex(1)->x(), e.getVertex(1)->y(), v->getTag(), v->getTag());
      }
    }
    fprintf(myFile, "};"); fclose(myFile);

    if(iConv == NCONV-1) {
      // Compute geodesics
      myFile = fopen("geodesics.pos", "w");
      fprintf(myFile, "View \" geodesics \"{\n");
      double x0[2] = {vertices[0].x(), vertices[0].y() };
      double x1[2] = {vertices[N_SUBTRIANGLES].x(), vertices[N_SUBTRIANGLES].y() };
      double x2[2] = {vertices[vertices.size()-1].x(), vertices[vertices.size()-1].y() };
      fprintf(myFile, "SP(%g,%g,0.){%d};\n", x0[0], x0[1], 1);
      fprintf(myFile, "SP(%g,%g,0.){%d};\n", x1[0], x1[1], 1);
      fprintf(myFile, "SP(%g,%g,0.){%d};\n", x2[0], x2[1], 1);

      int maxIter = 100;
      double tol = 1e-5;
      double ds = tol/2.;
      double length, error;

      std::vector<std::vector<SPoint2>> geodesics(3);
      geodesics[0] = geodesicBetweenTwoPoints_arrayPtr(x0, x1, getMetric, maxIter, tol, ds, length, error);
      feInfo("Length of geodesic = %+-1.3e", length);
      feInfo("Distance to target = %+-1.3e", error);
      geodesics[1] = geodesicBetweenTwoPoints_arrayPtr(x1, x2, getMetric, maxIter, tol, ds, length, error);
      feInfo("Length of geodesic = %+-1.3e", length);
      feInfo("Distance to target = %+-1.3e", error);
      geodesics[2] = geodesicBetweenTwoPoints_arrayPtr(x2, x0, getMetric, maxIter, tol, ds, length, error);
      feInfo("Length of geodesic = %+-1.3e", length);
      feInfo("Distance to target = %+-1.3e", error);

      for(size_t i = 0; i < 3; ++i) {
        for(SPoint2 p : geodesics[i]) { fprintf(myFile, "SP(%+-1.10e,%+-1.10e,0.){%d};\n", p[0], p[1], 1); }
      }
      fprintf(myFile, "};"); fclose(myFile);

      // Estimate the error betweem the transformation and the geodesics
      // Measure the distance from each edge point to the geodesic
      double E[3] = {0., 0., 0.};
      Etot[iConv] = 0.;
      for(size_t i = 0; i < 3; ++i) {
        for(auto ind : boundaryVertices[i]) {
          double x[2] = {vertices[ind].x(), vertices[ind].y()};
          double dist = distGeodesicToPoint(x, geodesics[i]);
          E[i] += fabs(dist);
          Etot[iConv] = fmax(fabs(dist), Etot[iConv]);
        }
      }
      feInfo("E = %+-1.4e - %+-1.4e - %+-1.4e - %+-1.4e", E[0], E[1], E[2], Etot[iConv]);
    }
  }


  // Compute the convergence rate
  for(int i = 1; i < NCONV; ++i) {
    rate[i] = -log(Etot[i] / Etot[i-1]) / log(h[i] / h[i-1]);
  }
  printf("%12s \t %12s \t %12s \n", "h", "E", "rate");
  for(int i = 0; i < NCONV; ++i) {
    printf("%12.6e \t %12.6e \t %12.6e\n", h[i], Etot[i], rate[i]);
  }

  return 0;
}