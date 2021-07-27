#include "feCurvedAdapt.h"

#include <iostream>
#include <poly34.h>
#include "SBoundingBox3d.h"

#ifdef HAVE_GMSH
#include "gmsh.h"

static int U_ = 0;
static int U_X = 1;
static int U_Y = 2;
static int U_XX = 3;
static int U_XY = 4;
static int U_YX = 5;
static int U_YY = 6;
static int U_XXX = 7;
static int U_XXY = 8;
static int U_XYX = 9;
static int U_XYY = 10;
static int U_YXX = 11;
static int U_YXY = 12;
static int U_YYX = 13;
static int U_YYY = 14;

static void PROBE(int TAG, double x, double y, double z, std::vector<double> &val,
                  const int step = -1, const int numComp = -1) {
  val.clear();
  gmsh::view::probe(TAG, x, y, z, val, step, numComp, 0, 1.e-3);
  if(val.empty()) val.push_back(0.);
}

static double f(double x, double y) {
  std::vector<double> u;
  PROBE(U_, x, y, 0.0, u);
  return u[0];

  // gmsh::model::setCurrent(nameModel);
  // std::vector<double> u,v;
  // PROBE(U_,x,y,0.0,u);
  // PROBE(V_,x,y,0.0,v);
  // gmsh::model::setCurrent("rect");
  // return u[0]*u[0]+v[0]*v[0];
}

static double fxxx(double x, double y) {
  std::vector<double> u;
  PROBE(U_XXX, x, y, 0.0, u);
  return u[0];
  // gmsh::model::setCurrent(nameModel);
  // std::vector<double> u,v;
  // PROBE(U_,x,y,0.0,u);
  // PROBE(V_,x,y,0.0,v);
  // std::vector<double> ux, vx;
  // PROBE(U_X,x,y,0.0,ux);
  // PROBE(V_X,x,y,0.0,vx);
  // std::vector<double> uxx, vxx;
  // PROBE(U_XX,x,y,0.0,uxx,0,1);
  // PROBE(V_XX,x,y,0.0,vxx,0,1);

  // std::vector<double> uxxx, vxxx;
  // PROBE(U_XXX,x,y,0.0,uxxx);
  // PROBE(V_XXX,x,y,0.0,vxxx);
  // gmsh::model::setCurrent("rect");
  // return 6*(uxx[0]*ux[0] +vxx[0]*vx[0]) + 2* (u[0]*uxxx[0]+v[0]*vxxx[0]); ;
}

static double fyyy(double x, double y) {
  std::vector<double> u;
  PROBE(U_YYY, x, y, 0.0, u);
  return u[0];
  // gmsh::model::setCurrent(nameModel);
  // std::vector<double> u,v;
  // PROBE(U_,x,y,0.0,u);
  // PROBE(V_,x,y,0.0,v);
  // std::vector<double> uy, vy;
  // PROBE(U_Y,x,y,0.0,uy);
  // PROBE(V_Y,x,y,0.0,vy);
  // std::vector<double> uyy, vyy;
  // PROBE(U_YY,x,y,0.0,uyy,0,1);
  // PROBE(V_YY,x,y,0.0,vyy,0,1);
  // std::vector<double> uyyy, vyyy;
  // PROBE(U_YYY,x,y,0.0,uyyy);
  // PROBE(V_YYY,x,y,0.0,vyyy);
  // gmsh::model::setCurrent("rect");
  // return 6*(uyy[0]*uy[0] +vyy[0]*vy[0]) + 2* (u[0]*uyyy[0]+v[0]*vyyy[0]); ;
}

static double fxyy(double x, double y) {
  std::vector<double> u;
  PROBE(U_XYY, x, y, 0.0, u);
  return u[0];

  // gmsh::model::setCurrent(nameModel);
  // std::vector<double> u,v;
  // PROBE(U_,x,y,0.0,u);
  // PROBE(V_,x,y,0.0,v);
  // std::vector<double> uy, vy, ux, vx;
  // PROBE(U_Y,x,y,0.0,uy);
  // PROBE(V_Y,x,y,0.0,vy);
  // PROBE(U_X,x,y,0.0,ux);
  // PROBE(V_X,x,y,0.0,vx);
  // std::vector<double> uyx, vyx, uyy;
  // std::vector<double> uxy, vxy, vyy;
  // PROBE(U_YY,x,y,0.0,uyy);
  // PROBE(U_XY,x,y,0.0,uxy);
  // PROBE(U_YX,x,y,0.0,uyx);
  // PROBE(V_YY,x,y,0.0,vyy);
  // PROBE(V_XY,x,y,0.0,vxy);
  // PROBE(V_YX,x,y,0.0,vyx);
  // double UXY = 0.5*(uxy[0]+uyx[0]);
  // double VXY = 0.5*(vxy[0]+vyx[0]);
  // std::vector<double> uxyy, uyxy,uyyx ;
  // std::vector<double> vxyy, vyxy,vyyx ;
  // PROBE(U_YYX,x,y,0.0,uyyx);
  // PROBE(U_XYY,x,y,0.0,uxyy);
  // PROBE(U_YXY,x,y,0.0,uyxy);
  // PROBE(V_YYX,x,y,0.0,vyyx);
  // PROBE(V_XYY,x,y,0.0,vxyy);
  // PROBE(V_YXY,x,y,0.0,vyxy);
  // double UXYY = (uxyy[0]+uyxy[0]+uyyx[0])/3.0;
  // double VXYY = (vxyy[0]+vyxy[0]+vyyx[0])/3.0;
  // gmsh::model::setCurrent("rect");
  // return 4 * (UXY*uy[0] + VXY*vy[0]) + 2 * (uy[0]*uyy[0]+vy[0]*vyy[0]) +
  //   2 * (u[0]*UXYY+v[0]*VXYY);
}

static double fxxy(double x, double y) {
  std::vector<double> u;
  PROBE(U_XXY, x, y, 0.0, u);
  return u[0];

  // gmsh::model::setCurrent(nameModel);
  // std::vector<double> u,v;
  // PROBE(U_,x,y,0.0,u);
  // PROBE(V_,x,y,0.0,v);
  // std::vector<double> uy, vy, ux, vx;
  // PROBE(U_Y,x,y,0.0,uy);
  // PROBE(V_Y,x,y,0.0,vy);
  // PROBE(U_X,x,y,0.0,ux);
  // PROBE(V_X,x,y,0.0,vx);
  // std::vector<double> uyx, vyx, uxx;
  // std::vector<double> uxy, vxy, vxx;
  // PROBE(U_XX,x,y,0.0,uxx);
  // PROBE(U_XY,x,y,0.0,uxy);
  // PROBE(U_YX,x,y,0.0,uyx);
  // PROBE(V_XX,x,y,0.0,vxx);
  // PROBE(V_XY,x,y,0.0,vxy);
  // PROBE(V_YX,x,y,0.0,vyx);
  // double UXY = 0.5*(uxy[0]+uyx[0]);
  // double VXY = 0.5*(vxy[0]+vyx[0]);
  // std::vector<double> uxxy, uyxx,uxyx ;
  // std::vector<double> vxxy, vyxx,vxyx ;
  // PROBE(U_XXY,x,y,0.0,uxxy);
  // PROBE(U_YXX,x,y,0.0,uyxx);
  // PROBE(U_XYX,x,y,0.0,uxyx);
  // PROBE(V_XXY,x,y,0.0,vxxy);
  // PROBE(V_YXX,x,y,0.0,vyxx);
  // PROBE(V_XYX,x,y,0.0,vxyx);
  // gmsh::model::setCurrent("rect");
  // double UXXY = (uxxy[0]+uyxx[0]+uxyx[0])/3.0;
  // double VXXY = (vxxy[0]+vyxx[0]+vxyx[0])/3.0;
  // gmsh::model::setCurrent("rect");
  // return 4 * (UXY*ux[0] + VXY*vx[0]) + 2 * (uy[0]*uxx[0]+vy[0]*vxx[0]) +
  //   2 * (u[0]*UXXY+v[0]*VXXY);
}

static double dttt(const double x, const double y, double C, double S) {
  const double c111 = fxxx(x, y);
  const double c222 = fyyy(x, y);
  const double c112 = fxxy(x, y);
  const double c122 = fxyy(x, y);
  return C * C * C * c111 + S * S * S * c222 + 3 * C * C * S * c112 + 3 * C * S * S * c122;
}

void feCurvedAdapt::computeMetricP2(double x, double y, double lMin, double lMax, double eps,
                                    double &g00, double &g01, double &g11, FILE *F, double &C,
                                    double &S) {
  // void feCurvedAdapt::computeMetricP2 (int node, double lMin, double lMax, double eps, double
  // &g00, double &g01, double &g11, FILE *F, double &C, double &S){
  double c111 = fxxx(x, y);
  const double c222 = fyyy(x, y);
  const double c112 = fxxy(x, y);
  const double c122 = fxyy(x, y);

  // std::cout<<x<<" "<<y<<" : "<<c111<<" "<<c222<<" "<<c112<<" "<<c122<<std::endl;

  if(c111 == 0) {
    printf("coucou\n");
    c111 = .0001;
  }

  double c[3];
  int nRoots = SolveP3(c, c112 / c111, c122 / c111, c222 / c111);

  int myRoot = 0;

  if(nRoots == 3) {
    double maxdttt_ = 0;
    for(int i = 0; i < 3; i++) {
      C = c[i];
      S = 1;
      double L = sqrt(C * C + S * S);
      C /= L;
      S /= L;
      double dttt_ = fabs(dttt(x, y, -S, C));
      if(dttt_ > maxdttt_) {
        myRoot = i;
        maxdttt_ = dttt_;
      }
    }
  }

  C = c[myRoot];
  S = 1;

  const double L = sqrt(C * C + S * S);
  C /= L;
  S /= L;
  return;
}

void smoothDirections(std::map<size_t, double> &C, std::map<size_t, double> &S) {
  std::vector<int> elementTypes;
  std::vector<std::vector<std::size_t> > elementTags;
  std::vector<std::vector<std::size_t> > nodeTags;
  gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, 2);
  std::multimap<size_t, size_t> graph;
  std::set<size_t> nodes_to_treat;
  std::set<size_t> nodes;

  for(size_t i = 0; i < elementTags[0].size(); i++) {
    size_t n0 = nodeTags[0][3 * i + 0];
    size_t n1 = nodeTags[0][3 * i + 1];
    size_t n2 = nodeTags[0][3 * i + 2];
    graph.insert(std::make_pair(n0, n1));
    graph.insert(std::make_pair(n1, n2));
    graph.insert(std::make_pair(n2, n0));
    nodes.insert(n0);
    nodes.insert(n1);
    nodes.insert(n2);
  }

  double threshold = 0.1;

  for(auto n : nodes) {
    double c = C[n];
    double s = S[n];
    for(auto it = graph.lower_bound(n); it != graph.upper_bound(n); ++it) {
      size_t neigh = it->second;
      double cn = C[neigh];
      double sn = S[neigh];
      double pv = fabs(cn * c + sn * s);
      if(pv > threshold && pv < 1. - threshold) nodes_to_treat.insert(n);
    }
  }
  int iter = 0;
  while(iter++ < 100) {
    for(auto n : nodes_to_treat) {
      double c4 = 0;
      double s4 = 0;
      for(auto it = graph.lower_bound(n); it != graph.upper_bound(n); ++it) {
        size_t neigh = it->second;
        double cn = C[neigh];
        double sn = S[neigh];
        double theta = atan2(sn, cn);
        c4 += cos(4 * theta);
        s4 += sin(4 * theta);
      }
      double theta = 0.25 * atan2(s4, c4);
      C[n] = cos(theta);
      S[n] = sin(theta);
    }
  }

  FILE *f = fopen("dirs.pos", "w");
  fprintf(f, "View\"Dirs\"{\n");

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

    fprintf(f, "VP(%g,%g,0){%g,%g,0};", coord[0], coord[1], c, s);
  }

  fprintf(f, "};\n");
  fclose(f);
}

double ERROR_SQUARED_P2(double *xa, double *xb, double *xc, double *xab, double *xbc, double *xca) {
  int triangleP2 = gmsh::model::mesh::getElementType("Triangle", 2);
  std::vector<double> localCoord;
  std::vector<double> weights;
  gmsh::model::mesh::getIntegrationPoints(triangleP2, "Gauss12", localCoord, weights);
  std::vector<double> basisFunctions;
  int numComponents, numOrientations;
  gmsh::model::mesh::getBasisFunctions(triangleP2, localCoord, "Lagrange", numComponents,
                                       basisFunctions, numOrientations);

  std::vector<double> gradBasisFunctions;
  gmsh::model::mesh::getBasisFunctions(triangleP2, localCoord, "GradLagrange", numComponents,
                                       gradBasisFunctions, numOrientations);

  double F[6] = {f(xa[0], xa[1]),   f(xb[0], xb[1]),   f(xc[0], xc[1]),
                 f(xab[0], xab[1]), f(xbc[0], xbc[1]), f(xca[0], xca[1])};
  double X[6] = {xa[0], xb[0], xc[0], xab[0], xbc[0], xca[0]};
  double Y[6] = {xa[1], xb[1], xc[1], xab[1], xbc[1], xca[1]};
  double e2 = 0;
  for(size_t i = 0; i < weights.size(); i++) {
    double interpolated = 0, x = 0, y = 0, dxdu = 0, dxdv = 0, dydu = 0, dydv = 0;

    for(size_t j = 0; j < 6; j++) {
      x += basisFunctions[6 * i + j] * X[j];
      y += basisFunctions[6 * i + j] * Y[j];
      dxdu += gradBasisFunctions[3 * (6 * i + j) + 0] * X[j];
      dxdv += gradBasisFunctions[3 * (6 * i + j) + 1] * X[j];
      dydu += gradBasisFunctions[3 * (6 * i + j) + 0] * Y[j];
      dydv += gradBasisFunctions[3 * (6 * i + j) + 1] * Y[j];
      interpolated += basisFunctions[6 * i + j] * F[j];
    }
    double detJ = fabs(dxdu * dydv - dxdv * dydu);
    double exact = f(x, y);
    double diff = exact - interpolated;
    e2 += weights[i] * diff * diff * detJ;
  }
  return e2;
}

// static bool rtreeCallback(int id, void *ctx){
//   std::vector<int> *vec = reinterpret_cast<std::vector<int> *>(ctx);
//   vec->push_back(id);
//   return true;
// }
#endif

feCurvedAdapt::feCurvedAdapt(feMesh *mesh, std::vector<feRecovery *> &recovery) : _rec(recovery) {
#ifdef HAVE_GMSH

  // std::vector<Vertex> &vertices = mesh->getVertices();

  // // Add all interior elements to an rtree
  // for(feCncGeo *cnc : mesh->getCncGeo()){
  //   if(cnc->getDim() == mesh->getDim()){
  //     int nElm = cnc->getNbElm();
  //     int nNodePerElm = cnc->getNbNodePerElem(); // All nodes are used to create the bounding
  //     box, not just the vertices for(int iElm = 0; iElm < nElm; ++iElm){
  //       SBoundingBox3d bbox;
  //       printf("element %2d : %2d - %2d - %2d\n", iElm, cnc->getNodeConnectivity(iElm, 0),
  //       cnc->getNodeConnectivity(iElm, 1), cnc->getNodeConnectivity(iElm, 2)); for(int iNode = 0;
  //       iNode < nNodePerElm; ++iNode){
  //         int node = cnc->getNodeConnectivity(iElm, iNode);
  //         SPoint3 pt(vertices[node].x(), vertices[node].y(), vertices[node].z());
  //         bbox += pt;
  //         _rtree.Insert((double*)(bbox.min()), (double*)(bbox.max()), iElm);
  //       }
  //     }
  //   }
  // }

  // std::vector<int> &vertices = recovery[0]->getVertices();

  // for(int v : vertices){
  // computeMetricP2();
  // }

  gmsh::initialize();
  std::string metricMeshName = "../../data/square2Msh2Adapt.msh";
  gmsh::open(metricMeshName);
  std::string nameModel;
  gmsh::model::getCurrent(nameModel);
  gmsh::model::add("rect");
  gmsh::model::setCurrent("rect");
  gmsh::merge(metricMeshName);

  // get the nodes
  std::vector<std::size_t> nodeTags;
  std::vector<double> coord;
  std::vector<double> parametricCoord;
  int tag = 1;
  gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, -1, tag, true, false);

  double modelSize = 1.;
  const double lMin = modelSize / 1000.;
  const double lMax = modelSize / 6.;

  int viewTag = gmsh::view::add(":metricP2_straight");
  int viewTagF = gmsh::view::add("F");
  std::vector<std::vector<double> > data;
  std::vector<std::vector<double> > dataF;

  double g00, g01, g11;

  FILE *F = fopen("metric.pos", "w");
  fprintf(F, "View \" \"{\n");

  std::map<size_t, double> COS;
  std::map<size_t, double> SIN;

  double eps = 1e-1;

  for(size_t i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];

    double C, S;
    computeMetricP2(x, y, lMin, lMax, eps, g00, g01, g11, F, C, S);

    COS[nodeTags[i]] = C;
    SIN[nodeTags[i]] = S;
  }

  smoothDirections(COS, SIN);

  for(size_t i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];

    double C = COS[nodeTags[i]];
    double S = SIN[nodeTags[i]];

    fprintf(F, "VP(%g,%g,%g){%g,%g,%g};\n", x, y, 0., C, S, 0.);
    fprintf(F, "VP(%g,%g,%g){%g,%g,%g};\n", x, y, 0., -S, C, 0.);

    const double dttt0_ = fabs(dttt(x, y, C, S));
    const double dttt1_ = fabs(dttt(x, y, -S, C));

    double l0 = dttt0_ != 0. ? pow(6 * eps / dttt0_, 0.3333) : lMax;
    double l1 = dttt1_ != 0. ? pow(6 * eps / dttt1_, 0.3333) : lMax;

    l0 = std::min(l0, lMax);
    l0 = std::max(l0, lMin);
    l1 = std::min(l1, lMax);
    l1 = std::max(l1, lMin);

    double h0 = 1. / (l0 * l0);
    double h1 = 1. / (l1 * l1);

    double g00 = C * C * h0 + S * S * h1;
    double g11 = S * S * h0 + C * C * h1;
    double g01 = S * C * (h1 - h0);

    std::vector<double> v(9);
    std::vector<double> vF(1);

    v[0] = g00;
    v[1] = -g01;
    v[2] = 0;

    v[3] = -g01;
    v[4] = g11;
    v[5] = 0;

    v[6] = 0;
    v[7] = 0;
    v[8] = f(x, y); // export f as well.
    vF[0] = f(x, y);
    data.push_back(v);
    dataF.push_back(vF);
  }
  fprintf(F, "};\n");
  fclose(F);

  gmsh::view::addModelData(viewTag, 0, "rect", "NodeData", nodeTags, data);
  gmsh::view::addModelData(viewTagF, 0, "rect", "NodeData", nodeTags, dataF);

  printf("dataf.size = %lu\n", dataF.size());

  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);

  //  gmsh::write("metric.msh");
  gmsh::view::write(viewTagF, "sol.msh");
  gmsh::view::write(viewTag, "metric.msh");
  //  gmsh::view::write(viewTagF,"metric.msh",true);

  system("mmg2d -hgrad 3 metric.msh");

  gmsh::option::setNumber("Mesh.MshFileVersion", 4);
  gmsh::model::add("rect_adapt");
  gmsh::model::setCurrent("rect_adapt");
  gmsh::merge("metric.o.msh");

  int computePointsUsingScaledCrossFieldPlanarP2(
    const char *modelForMetric, const char *modelForMesh, int VIEW_TAG, int faceTag,
    std::vector<double> &pts,
    double er(double *, double *, double *, double *, double *, double *));

  bool vazy = true;

  std::vector<double> pts;
  if(vazy)
    computePointsUsingScaledCrossFieldPlanarP2("rect", "rect_adapt", viewTag, 4, pts,
                                               ERROR_SQUARED_P2);
  gmsh::model::setCurrent("rect_adapt");

  if(!vazy) gmsh::model::mesh::setOrder(2);
  gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, -1, -1, true);

  int viewTagF_adapt = gmsh::view::add("F_adapt");
  std::vector<std::vector<double> > dataF_adapt;

  for(size_t i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];
    std::vector<double> vF(1);
    vF[0] = f(x, y);
    dataF_adapt.push_back(vF);
  }
  gmsh::model::setCurrent("rect_adapt");
  gmsh::view::addModelData(viewTagF_adapt, 0, "rect_adapt", "NodeData", nodeTags, dataF_adapt);
  gmsh::view::write(viewTagF_adapt, "sol_adapt.msh");

  std::vector<int> elementTypes;

  gmsh::model::setCurrent("rect_adapt");

  {
    FILE *ERR = fopen("ERROR.pos", "w");
    fprintf(ERR, "View\"ERROR\"{\n");
    int triangleP2 = gmsh::model::mesh::getElementType("Triangle", 2);
    std::vector<double> localCoord;
    std::vector<double> weights;
    gmsh::model::mesh::getIntegrationPoints(triangleP2, "Gauss12", localCoord, weights);
    std::vector<double> basisFunctions;
    int numComponents, numOrientations;
    gmsh::model::mesh::getBasisFunctions(triangleP2, localCoord, "Lagrange", numComponents,
                                         basisFunctions, numOrientations);
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t> > elementTags;
    std::vector<std::vector<std::size_t> > nodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, 2);
    double Error = 0.0;

    for(size_t e = 0; e < elementTags[0].size(); e++) {
      int eDim, eTag;
      std::vector<double> xa, xb, xc, xab, xbc, xca, pc;
      gmsh::model::mesh::getNode(nodeTags[0][6 * e + 0], xa, pc, eDim, eTag);
      gmsh::model::mesh::getNode(nodeTags[0][6 * e + 1], xb, pc, eDim, eTag);
      gmsh::model::mesh::getNode(nodeTags[0][6 * e + 2], xc, pc, eDim, eTag);
      gmsh::model::mesh::getNode(nodeTags[0][6 * e + 3], xab, pc, eDim, eTag);
      gmsh::model::mesh::getNode(nodeTags[0][6 * e + 4], xbc, pc, eDim, eTag);
      gmsh::model::mesh::getNode(nodeTags[0][6 * e + 5], xca, pc, eDim, eTag);
      double ErrorElement = ERROR_SQUARED_P2(&xa[0], &xb[0], &xc[0], &xab[0], &xbc[0], &xca[0]);
      fprintf(ERR,
              "ST2(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g){%g,%g,%g,%g,%g,%g};\n",
              xa[0], xa[1], 0.0, xb[0], xb[1], 0.0, xc[0], xc[1], 0.0, xab[0], xab[1], 0.0, xbc[0],
              xbc[1], 0.0, xca[0], xca[1], 0.0, ErrorElement, ErrorElement, ErrorElement,
              ErrorElement, ErrorElement, ErrorElement);
      Error += ErrorElement;
    }
    fprintf(ERR, "};\n");
    fclose(ERR);
    printf("Error = %12.5E\n", sqrt(Error));
  }
#else
  printf("In feCurvedAdapt : Error : Gmsh is required to compute curved meshes.\n");
#endif
}

feCurvedAdapt::~feCurvedAdapt() {
  // for(feGridFunction *gF : _derivatives)
  // delete gF;
}