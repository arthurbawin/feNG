#include "feMetricTools.h"
#include "feSimplex.h"

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
// #include "eigen3/Eigen/Eigen"
// #include "eigen3/unsupported/Eigen/MatrixFunctions"
#include "../contrib/Eigen/Eigen"
#include "../contrib/unsupported/Eigen/MatrixFunctions"

double evaluateFieldFromRecovery(int indexDerivative, feRecovery *rec, double *x)
{
  return rec->evalDerivative(indexDerivative, x);
}

double evaluateFieldFromRecoveryCallback(int indexDerivative, void *recUserPtr, double *x)
{
  return ((feRecovery*) recUserPtr)->evalDerivative(indexDerivative, x);
}

double f(feRecovery *rec, double *x) { return evaluateFieldFromRecovery(0, rec, x); }
double fx(feRecovery *rec, double *x) { return evaluateFieldFromRecovery(1, rec, x); }
double fy(feRecovery *rec, double *x) { return evaluateFieldFromRecovery(2, rec, x); }
double fxx(feRecovery *rec, double *x) { return evaluateFieldFromRecovery(3, rec, x); }
double fxy(feRecovery *rec, double *x) { return evaluateFieldFromRecovery(4, rec, x); }
double fyx(feRecovery *rec, double *x) { return evaluateFieldFromRecovery(5, rec, x); }
double fyy(feRecovery *rec, double *x) { return evaluateFieldFromRecovery(6, rec, x); }
double fxxx(feRecovery *rec, double *x) { return evaluateFieldFromRecovery(7, rec, x); }
double fxxy(feRecovery *rec, double *x) { return evaluateFieldFromRecovery(8, rec, x); }
double fxyx(feRecovery *rec, double *x) { return evaluateFieldFromRecovery(9, rec, x); }
double fxyy(feRecovery *rec, double *x)
{
  return evaluateFieldFromRecovery(10, rec, x);
}
double fyxx(feRecovery *rec, double *x)
{
  return evaluateFieldFromRecovery(11, rec, x);
}
double fyxy(feRecovery *rec, double *x)
{
  return evaluateFieldFromRecovery(12, rec, x);
}
double fyyx(feRecovery *rec, double *x)
{
  return evaluateFieldFromRecovery(13, rec, x);
}
double fyyy(feRecovery *rec, double *x)
{
  return evaluateFieldFromRecovery(14, rec, x);
}

double dtt(double *x, double C, double S, feRecovery *rec)
{
  const double c11 = fxx(rec, x);
  const double c12 = (fxy(rec, x) + fyx(rec, x)) / 2.;
  const double c22 = fyy(rec, x);
  return C * C * c11 + S * S * c22 + 2. * C * S * c12;
}

double dttt(double *x, double C, double S, feRecovery *rec, int direction)
{
  const double c111 = fxxx(rec, x);
  const double c112 = (fxxy(rec, x) + fxyx(rec, x) + fyxx(rec, x)) / 3.;
  const double c122 = (fxyy(rec, x) + fyyx(rec, x) + fyxy(rec, x)) / 3.;
  const double c222 = fyyy(rec, x);

  const double c11 = fxx(rec, x);
  const double c12 = (fxy(rec, x) + fyx(rec, x)) / 2.;
  const double c22 = fyy(rec, x);

  const double c1 = fx(rec, x);
  const double c2 = fy(rec, x);

  double kappa1 =
    (-c2 * c2 * c11 + 2.0 * c1 * c2 * c12 - c1 * c1 * c22) / (pow(c1 * c1 + c2 * c2, 3.0 / 2.0));
  double kappa2 =
    (c1 * c2 * (c22 - c11) + (c1 * c1 - c2 * c2) * c12) / (pow(c1 * c1 + c2 * c2, 3.0 / 2.0));

  double g11 = C;
  double g12 = S;
  double g21 = -S;
  double g22 = C;

  if(direction == 0) {
    // return C * C * C * c111 + S * S * S * c222 + 3. * C * C * S * c112 + 3. * C * S * S * c122
    // + 3.0 * kappa1 * (g11*g21*c11 + g11*g22*c12 + g12*g21*c12 + g12*g22*c22);
    return g11 * g11 * g11 * c111 + g12 * g12 * g12 * c222 + 3. * g11 * g11 * g12 * c112 +
           3. * g11 * g12 * g12 * c122 +
           3.0 * kappa1 * (g11 * g21 * c11 + g11 * g22 * c12 + g12 * g21 * c12 + g12 * g22 * c22);
  } else if(direction == 1) {
    // return C * C * C * c111 + S * S * S * c222 + 3. * C * C * S * c112 + 3. * C * S * S * c122
    // + 3.0 * kappa2 * (g11*g21*c11 + g21*g12*c12 + g22*g11*c12 + g12*g22*c22);
    return g21 * g21 * g21 * c111 + g22 * g22 * g22 * c222 + 3. * g21 * g21 * g22 * c112 +
           3. * g21 * g22 * g22 * c122 +
           3.0 * kappa2 * (g11 * g21 * c11 + g21 * g12 * c12 + g22 * g11 * c12 + g12 * g22 * c22);
  } else {
    // Compute size along a straight edge
    return C * C * C * c111 + S * S * S * c222 + 3. * C * C * S * c112 + 3. * C * S * S * c122;
  }
}

double sech(double x) { return 1. / cosh(x); }

// For the solution u = (1. + tanh(a*(r-r0)))/2.
double dtttAnalytical(double x, double y, double C, double S, feRecovery *rec,
                      int direction)
{
  double a = 10.0;
  double x0 = 0.2;
  double y0 = 0.4;
  double r0 = 0.2;

  double r = sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0));
  double drdx = (x - x0) / r;
  double drdy = (y - y0) / r;

  double dudr = a / 2. * sech(a * (r - r0)) * sech(a * (r - r0));
  double d2udr2 = -a * a * tanh(a * (r - r0)) * sech(a * (r - r0)) * sech(a * (r - r0));
  double d3udr3 =
    -a * a * a * sech(a * (r - r0)) * sech(a * (r - r0)) *
    (sech(a * (r - r0)) * sech(a * (r - r0)) - 2. * tanh(a * (r - r0)) * tanh(a * (r - r0)));

  const double c1 = dudr * drdx;
  const double c2 = dudr * drdy;

  const double c11 = d2udr2 * drdx * drdx;
  const double c12 = d2udr2 * drdx * drdy;
  const double c22 = d2udr2 * drdy * drdy;

  const double c111 = d3udr3 * drdx * drdx * drdx;
  const double c222 = d3udr3 * drdy * drdy * drdy;
  const double c112 = d3udr3 * drdx * drdx * drdy;
  const double c122 = d3udr3 * drdx * drdy * drdy;

  double kappa1 =
    (-c2 * c2 * c11 + 2.0 * c1 * c2 * c12 - c1 * c1 * c22) / (pow(c1 * c1 + c2 * c2, 3.0 / 2.0));
  double kappa2 =
    (c1 * c2 * (c22 - c11) + (c1 * c1 - c2 * c2) * c12) / (pow(c1 * c1 + c2 * c2, 3.0 / 2.0));

  // double A, b;
  // A = c1;
  // b = c2;

  // double normGrad = sqrt(A * A + b * b);

  // if(normGrad > 1e-5) {
  //   double theta1 = atan2(b, A);
  //   C = cos(theta1);
  //   S = sin(theta1);
  // } else {
  //   // Gradient norm is too small : directions will be smoothed
  //   C = 0.;
  //   S = 0.;
  // }

  double g11 = C;
  double g12 = S;
  double g21 = -S;
  double g22 = C;

  // if(x == 0.25 && y == 0.25){
  //   feInfo("g11 = %4.4f - g11 = %4.4f - g21 = %4.4f - g22 = %4.4f", g11, g12, g21, g22);
  //   feInfo("k1 = %8.8e - k2 = %8.8e", kappa1, kappa2);
  //   feInfo("C1 = %8.8e - H1 = %8.8e", g11 * g11 * g11 * c111 + g12 * g12 * g12 * c222 + 3. * g11
  //   * g11 * g12 * c112 +
  //          3. * g11 * g12 * g12 * c122, g11 * g21 * c11 + g11 * g22 * c12 + g12 * g21 * c12 + g12
  //          * g22 * c22);
  //   feInfo("C2 = %8.8e - H2 = %8.8e", g21 * g21 * g21 * c111 + g22 * g22 * g22 * c222 + 3. * g21
  //   * g21 * g22 * c112 +
  //          3. * g21 * g22 * g22 * c122, g11 * g21 * c11 + g21 * g12 * c12 + g22 * g11 * c12 + g12
  //          * g22 * c22);
  // }

  if(direction == 0) {
    // return C * C * C * c111 + S * S * S * c222 + 3. * C * C * S * c112 + 3. * C * S * S * c122
    // + 3.0 * kappa1 * (g11*g21*c11 + g11*g22*c12 + g12*g21*c12 + g12*g22*c22);
    return g11 * g11 * g11 * c111 + g12 * g12 * g12 * c222 + 3. * g11 * g11 * g12 * c112 +
           3. * g11 * g12 * g12 * c122 +
           3.0 * kappa1 * (g11 * g21 * c11 + g11 * g22 * c12 + g12 * g21 * c12 + g12 * g22 * c22);
  } else if(direction == 1) {
    // return C * C * C * c111 + S * S * S * c222 + 3. * C * C * S * c112 + 3. * C * S * S * c122
    // + 3.0 * kappa2 * (g11*g21*c11 + g21*g12*c12 + g22*g11*c12 + g12*g22*c22);
    return g21 * g21 * g21 * c111 + g22 * g22 * g22 * c222 + 3. * g21 * g21 * g22 * c112 +
           3. * g21 * g22 * g22 * c122 +
           3.0 * kappa2 * (g11 * g21 * c11 + g21 * g12 * c12 + g22 * g11 * c12 + g12 * g22 * c22);
  } else {
    // Compute size along a straight edge
    return C * C * C * c111 + S * S * S * c222 + 3. * C * C * S * c112 + 3. * C * S * S * c122;
  }
}

void computeDirectionFieldFromGradient(double *x, double &C, double &S, double tol,
                                       feRecovery *rec, FILE *F_grad, FILE *F_iso)
{
  double a, b;
  a = fx(rec, x);
  b = fy(rec, x);

  double normGrad = sqrt(a * a + b * b);

  if(normGrad > tol) {
    double theta1 = atan2(b, a);
    C = cos(theta1);
    S = sin(theta1);
    if(F_grad != nullptr) { fprintf(F_grad, "VP(%g,%g,%g){%g,%g,%g};\n", x[0], x[1], 0.,  C, S, 0.); }
    if(F_iso  != nullptr) { fprintf(F_iso , "VP(%g,%g,%g){%g,%g,%g};\n", x[0], x[1], 0., -S, C, 0.); }
  } else {
    // Gradient norm is too small : directions will be smoothed
    C = 1.;
    S = 0.;
  }
}

void computeDirectionFieldFromHessian(double *x, double &C, double &S, double tol,
                                      feRecovery *rec, FILE *F)
{
  double a, b, c;
  a = fxx(rec, x);
  b = (fxy(rec, x) + fyx(rec, x)) / 2.;
  c = fyy(rec, x);

  // Eigenvalues and spectral radius of the hessian
  double l1 = (a + c - sqrt(a * a - 2. * a * c + 4. * b * b + c * c)) / 2.;
  double l2 = (a + c + sqrt(a * a - 2. * a * c + 4. * b * b + c * c)) / 2.;

  double spectralRadius = fmax(fabs(l1), fabs(l2));

  if(spectralRadius > tol) {
    // First eigenvector
    double v11 = a - c + sqrt(a * a + 4. * b * b + c * c - 2 * a * c);
    double v12 = 2. * b;
    C = v11;
    S = v12;
    const double L = fmax(1e-10, sqrt(C * C + S * S));
    C /= L;
    S /= L;
    if(F != nullptr) {
      fprintf(F, "VP(%g,%g,%g){%g,%g,%g};\n", x[0], x[1], 0., C, S, 0.);
    }
  } else {
    // Hessian's eigenvalues are too small : directions will be smoothed
    C = 0.;
    S = 0.;
  }
}

void smoothDirections(std::map<size_t, double> &C, std::map<size_t, double> &S, int nIter, double tol)
{
#if defined(HAVE_GMSH)
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
        c4 += cos(4. * theta);
        s4 += sin(4. * theta);
      }
      double theta = 0.25 * atan2(s4, c4);
      C[n] = cos(theta);
      S[n] = sin(theta);
    }
  }

  FILE *fIso = fopen("isoAfterSmoothing.pos", "w"); fprintf(fIso, "View\"isoAfterSmoothing\"{\n");
  FILE *fGra = fopen("graAfterSmoothing.pos", "w"); fprintf(fGra, "View\"graAfterSmoothing\"{\n");

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

    fprintf(fGra, "VP(%g,%g,0){%g,%g,0};", coord[0], coord[1],  c, s);
    fprintf(fIso, "VP(%g,%g,0){%g,%g,0};", coord[0], coord[1], -s, c);
  }

  fprintf(fIso, "};"); fclose(fIso);
  fprintf(fGra, "};"); fclose(fGra);
#endif
}

void smoothSizes(std::map<size_t, double> &L1, std::map<size_t, double> &L2, int nIter,
                      double tol)
{
#if defined(HAVE_GMSH)
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

  double threshold = tol;

  for(auto n : nodes) {
    double l1 = L1[n];
    double l2 = L2[n];
    for(auto it = graph.lower_bound(n); it != graph.upper_bound(n); ++it) {
      size_t neigh = it->second;
      double l1n = L1[neigh];
      double l2n = L2[neigh];
      double grad1 = fabs(l1n - l1)/l1;
      double grad2 = fabs(l2n - l2)/l2;
      if(grad1 > threshold) nodes_to_treat.insert(n);
      if(grad2 > threshold) nodes_to_treat.insert(n);
    }
  }
  int iter = 0;
  while(iter++ < nIter) {
    for(auto n : nodes_to_treat) {
      double l1 = 0;
      double l2 = 0;
      int cnt = 0;
      for(auto it = graph.lower_bound(n); it != graph.upper_bound(n); ++it) {
        size_t neigh = it->second;
        l1 += L1[neigh];
        l2 += L2[neigh];
        cnt++;
      }
      L1[n] = l1/cnt;
      L2[n] = l2/cnt;
    }
  }
#endif
}

/* Compute metric using a brute-force method to solve optimization problem
   Original code from R. Kuate */
void metricHechtKuate(int nbpoints, double *x, double *y, double &A, double &B, double &C,
                      double epsilon, double *xNew, double *yNew)
{
  // C = 0.0;

  // int bool_assert = 1;

  // double epsilon0 = 1e-5, precision = 1e-18, delta = 1e-10;
  // double inf = 1e100;

  // double Rmin = 1e20, Rmax = 0;

  // int indiceX0 = 0;

  // // double* xNew = new double[nbpoints];
  // // double* yNew = new double[nbpoints];
  // // double *PPoint = new double[nbpoints];

  // for(int i = 0; i < nbpoints; i++) {
  //   Rmax = fmax(Rmax, sqrt(x[i] * x[i] + y[i] * y[i]));

  //   //---déplacement des points situées sur les axes--------------
  //   if(abs(x[i]) <= precision) {
  //     if(y[i] < 0) {
  //       x[i] = -delta;
  //       y[i] = -sqrt(pow(y[i], 2) - pow(x[i], 2));
  //     }
  //     if(y[i] > 0) {
  //       x[i] = delta;
  //       y[i] = sqrt(pow(y[i], 2) - pow(x[i], 2));
  //     }
  //   }

  //   if(abs(y[i]) <= precision) {
  //     if(x[i] < 0) {
  //       y[i] = -delta;
  //       x[i] = -sqrt(pow(x[i], 2) - pow(y[i], 2));
  //     }
  //     if(x[i] > 0) {
  //       y[i] = delta;
  //       x[i] = sqrt(pow(x[i], 2) - pow(y[i], 2));
  //     }
  //   }
  //   //-----------------------------------------------------------
  //   if(bool_assert) {
  //     // std::cout<<"foo : "<<abs(x[i]*y[i])<<std::endl;
  //     // std::cout<<"foo : "<<pow(precision,2)<<std::endl;
  //     assert(abs(x[i] * y[i]) >= pow(precision, 2));
  //   }

  //   if(Rmin > sqrt(x[i] * x[i] + y[i] * y[i])) {
  //     indiceX0 = i;
  //     Rmin = sqrt(x[i] * x[i] + y[i] * y[i]);
  //   }
  // }

  // //-------permutation des indices de la liste des points :
  // // ranger la liste en commençant par le point X0-------
  // for(int k = 0; k < nbpoints - indiceX0; k++) {
  //   xNew[k] = x[k + indiceX0];
  //   yNew[k] = y[k + indiceX0];
  // }
  // for(int k = nbpoints - indiceX0; k < nbpoints; k++) {
  //   xNew[k] = x[k - nbpoints + indiceX0];
  //   yNew[k] = y[k - nbpoints + indiceX0];
  // }
  // for(int i = 0; i < nbpoints; i++) {
  //   x[i] = xNew[i];
  //   y[i] = yNew[i];
  // }

  // //----------------------------------------------------------------

  // int test = -1;

  // double X0, Y0;
  // double bmin = 0., bmax = inf, b1, b2, aik = 0., bik = 0., cik = 0.;
  // double Xk = 0., Yk = 0., Ck = 0., Bk = 0., Xi = 0., Yi = 0., ri, detXY = 0., Ri, R0, r0;

  // X0 = x[0];
  // Y0 = y[0];
  // r0 = sqrt(x[0] * x[0] + y[0] * y[0]);
  // if(bool_assert) {
  //   // std::cout<<"Assert ON"<<std::endl;
  //   assert(r0 == Rmin);
  // }

  // // std::cout<<" Rmin = "<<Rmin<<" Rmax =  "<<Rmax<<std::endl;

  // double EPS = 0.0; // pour recuperer la valeur de epsilon0 optimale

  // int condition = -1;

  // if(r0 <= epsilon0) epsilon0 = r0 * epsilon0;

  // A = 1. / ((r0 - epsilon0) * (r0 - epsilon0));
  // B = A;

  // double epsilon0min = epsilon0;

  // // if(abs(Rmin-Rmax) > 1e-5){
  // if(abs(Rmin - Rmax) > precision) {
  //   // for(int ee=0; ee<neps-1; ee++){ //boucle sur epsilon0---------------
  //   // epsilon0= Tabepsilon[ee];
  //   if(r0 <= epsilon0) epsilon0 = r0 * epsilon0;
  //   if(bool_assert) {
  //     assert(r0 > epsilon0);
  //   }
  //   R0 = r0 / (r0 - epsilon0);

  //   for(int i = 1; i < nbpoints; i++) { // boucle sur chaque noeud
  //     Xi = x[i];
  //     Yi = y[i];
  //     ri = sqrt(x[i] * x[i] + y[i] * y[i]);

  //     if(ri <= epsilon) epsilon = ri * epsilon;

  //     if(bool_assert) {
  //       assert(ri > epsilon);
  //     }

  //     Ri = ri / (ri - epsilon);

  //     detXY = Xi * Y0 - Yi * X0;

  //     //------deplacement des points alignés avec l'origine et X0-----------
  //     if(abs(detXY) <= precision) {
  //       printf("Point %d - \t x = %10.15f - \t y =%10.15f - \t x0 = %10.15f - \t y0 =%10.15f - \t "
  //              "detXY = %10.15f \n",
  //              i, x[i], y[i], X0, Y0, detXY);
  //       Xi += delta;

  //       if(Yi < 0)
  //         Yi = -sqrt(pow(ri, 2) - pow(Xi, 2));
  //       else
  //         Yi = sqrt(pow(ri, 2) - pow(Xi, 2));
  //       x[i] = Xi;
  //       y[i] = Yi;

  //       std::cout << ri << std::endl;
  //       ri = sqrt(x[i] * x[i] + y[i] * y[i]);

  //       if(ri <= epsilon) epsilon = ri * epsilon;

  //       std::cout << "xi =" << x[i] << " yi =" << y[i] << " ri = " << ri << " epsilon = " << epsilon
  //                 << std::endl;
  //       if(bool_assert) {
  //         assert(ri > epsilon);
  //       }
  //       Ri = ri / (ri - epsilon);
  //     }

  //     detXY = Xi * Y0 - Yi * X0;

  //     if(bool_assert) {
  //       assert(abs(detXY) >= precision);
  //     }

  //     //-----racines du polynome en b à minimiser----------------------------
  //     double bb1 =
  //       (1. / pow(detXY, 2)) *
  //       (pow(X0 * Ri, 2) + pow(Xi * R0, 2) -
  //        2. * abs(Xi * X0) * sqrt(pow(R0 * Ri, 2) - pow(detXY / (Rmax * (r0 - epsilon0)), 2)));
  //     double bb2 =
  //       (1. / pow(detXY, 2)) *
  //       (pow(X0 * Ri, 2) + pow(Xi * R0, 2) +
  //        2. * abs(Xi * X0) * sqrt(pow(R0 * Ri, 2) - pow(detXY / (Rmax * (r0 - epsilon0)), 2)));
  //     //--fin----racines du polynome en b à minimiser--------------------

  //     bmax = fmin(bb2, pow(Rmax / pow((r0), 2), 2));
  //     bmin = fmax(1. / (Rmax * Rmax), bb1); // minoration de b
  //     double Cte = fmax(1e-9, (bmax - bmin) * 1e-9);
  //     bmin = bmin * (1. + Cte);
  //     bmax = bmax * (1. - Cte);

  //     // bornes de b-----------------------------------------------------------

  //     // cas:  majoration de c --------------------------------------------
  //     double Li = X0 * Xi * (pow(Rmax / pow(r0 - epsilon0min, 2), 2) - 1. / pow(Rmax, 2)) +
  //                 (pow(Ri * X0, 2) - pow(R0 * Xi, 2)) / detXY;
  //     double LiXY = Xi * Y0 + Yi * X0;

  //     if(abs(LiXY) >= precision) {
  //       condition = 1;

  //       if(Xi * X0 > 0) {
  //         if(LiXY > 0)
  //           bmin = fmax(bmin, -Li / LiXY);
  //         else
  //           bmax = fmin(bmax, -Li / LiXY);
  //       } else {
  //         if(LiXY < 0)
  //           bmin = fmax(bmin, -Li / LiXY);
  //         else
  //           bmax = fmin(bmax, -Li / LiXY);
  //       }
  //     } else {
  //       if(Li < 0)
  //         condition = 0;
  //       else
  //         condition = 1;
  //     }

  //     // cas  minoration de c --------------------------------------------
  //     Li = X0 * Xi * (-pow(Rmax / pow(r0 - epsilon0min, 2), 2) + 1. / pow(Rmax, 2)) +
  //          (pow(Ri * X0, 2) - pow(R0 * Xi, 2)) / detXY;
  //     LiXY = Xi * Y0 + Yi * X0;

  //     if(abs(LiXY) >= precision) {
  //       condition = 1;
  //       if(Xi * X0 > 0) {
  //         if(LiXY < 0)
  //           bmin = fmax(bmin, -Li / LiXY);
  //         else
  //           bmax = fmin(bmax, -Li / LiXY);
  //       } else {
  //         if(LiXY > 0)
  //           bmin = fmax(bmin, -Li / LiXY);
  //         else
  //           bmax = fmin(bmax, -Li / LiXY);
  //       }
  //     } else {
  //       if(Li > 0)
  //         condition = 0;
  //       else
  //         condition = 1;
  //     }

  //     if(condition) {
  //       //--cas : minoration de a-----------------------------------------------

  //       double Gi =
  //         ((Xi * Yi * R0 * R0 - X0 * Y0 * Ri * Ri) / detXY + Xi * X0 / (Rmax * Rmax)) / (Yi * Y0);

  //       if(Xi * X0 > 0) {
  //         if(Yi * Y0 > 0)
  //           bmin = fmax(bmin, Gi);
  //         else
  //           bmax = fmin(bmax, Gi);
  //       } else {
  //         if(Yi * Y0 < 0)
  //           bmin = fmax(bmin, Gi);
  //         else
  //           bmax = fmin(bmax, Gi);
  //       }

  //       // cas :majoration de a------------------------------------------------
  //       double Hi = (Xi * X0 * Rmax * Rmax / pow((r0 - epsilon0min), 4) +
  //                    (Xi * Yi * R0 * R0 - X0 * Y0 * Ri * Ri) / detXY) /
  //                   (Yi * Y0);
  //       if(Xi * X0 > 0) {
  //         if(Yi * Y0 > 0)
  //           bmax = fmin(bmax, Hi);
  //         else
  //           bmin = fmax(bmin, Hi);
  //       } else {
  //         if(Yi * Y0 < 0)
  //           bmax = fmin(bmax, Hi);
  //         else
  //           bmin = fmax(bmin, Hi);
  //       }
  //       //------fin bornes de b------------------------------------------------
  //       b2 = bmax;
  //       b1 = bmin;

  //       for(int k = 1; k < nbpoints; k++) { // on balaye les contraintes
  //         Xk = x[k];
  //         Yk = y[k];
  //         Bk = (Yk * Yk * Xi * X0 + Xk * (Xk * Yi * Y0 - Yk * (Yi * X0 + Xi * Y0))) / (Xi * X0);
  //         Ck = (X0 * Xi * detXY -
  //               Xk * (Xi * R0 * R0 * (Yk * Xi - Yi * Xk) + X0 * Ri * Ri * (-Yk * X0 + Y0 * Xk))) /
  //              (Xi * X0 * detXY);

  //         if(bool_assert) {
  //           assert(abs(Xi * X0 * Y0 * Yi * Xk * Yk) >= pow(precision, 5));
  //         }
  //         if(abs(Bk) > precision) { // non nul

  //           if(Bk <= 0)
  //             bmax = fmin(bmax, Ck / Bk);
  //           else
  //             bmin = fmax(bmin, Ck / Bk);

  //           if((bmax < b1) || (bmin > b2) || (bmin > bmax)) {
  //             test = 0;
  //             break;
  //           }

  //           else
  //             test = 1;
  //         } else {
  //           if(Ck > precision) {
  //             test = 0;
  //             break;
  //           } else // Ck<=0
  //             test = -1; // 1 peut etre
  //         }
  //       }

  //       if(test) {
  //         double a0 = -pow((detXY / (Xi * X0)), 2);
  //         double a1 = 2. * (pow(Ri / Xi, 2) + pow(R0 / X0, 2));
  //         if(((a0 * bmax + a1) * bmax) < ((a0 * bmin + a1) * bmin))
  //           bik = bmax;
  //         else
  //           bik = bmin;

  //         aik = (Ri * Ri * Y0 * X0 - R0 * R0 * Yi * Xi + bik * Yi * Y0 * detXY) / (detXY * Xi * X0);
  //         cik = (-Ri * Ri * X0 * X0 + R0 * R0 * Xi * Xi - bik * (Yi * X0 + Y0 * Xi) * detXY) /
  //               (detXY * Xi * X0);

  //         if(bool_assert) {
  //           assert((4. * aik * bik - cik * cik) >= 0.); // aire positive
  //           assert(abs((4. * aik * bik - cik * cik) - pow(2. / (Rmax * (r0 - epsilon0)), 2)) >
  //                  0); // aire positive
  //         }
  //         if((4. * aik * bik - cik * cik) <= (4. * A * B - C * C)) {
  //           A = aik;
  //           B = bik;
  //           C = cik;
  //           EPS = epsilon0;
  //         }
  //       } // if(test)
  //     } // if(condition)
  //   } // for(int i=1;i<nbpoints;i++)
  //   // } // for(int ee=0; ee<neps-1; ee++)
  // } // if(abs(Rmin-Rmax)>1e-5)
  // else {
  //   A = 1. / (Rmin * Rmin);
  //   B = A;
  //   C = 0.;
  // }
}

void metriqueSimplexe2D(int nPhi, std::vector<double> phi, std::vector<double> erreur, double &A,
                        double &B, double &C, int max_iter, std::vector<int> &expX,
                        std::vector<int> &expY, int dimRecovery, int degreeSolution)
{
  // double xi, yi, xj, yj, Q11, Q12, Q21, Q22;
  // double err, normeXj;
  // double L1, L2, L3, expL11, expL12, expL21, expL22;
  // double Aprev, Bprev, Cprev;
  // // Variables LAPACK
  // // char jobz='v', uplo='l';
  // // const int n=2, lwork=18;
  // // int info;
  // // double w[2], work[18];

  // SVector3 v0, v1, v2;

  // // std::vector<int> &expX = _recovery->getXExponentsRecovery();
  // // std::vector<int> &expY = _recovery->getYExponentsRecovery();

  // int nSimplex = 6, mSimplex = nPhi;

  // std::vector<double> AConstraint(nPhi * nSimplex), bConstraint(nPhi), cObjective(nSimplex);
  // cObjective[0] = -1.0;
  // cObjective[1] = 1.0;
  // cObjective[2] = 0.0;
  // cObjective[3] = 0.0;
  // cObjective[4] = -1.0;
  // cObjective[5] = 1.0;

  // A = 1.0; // Candidat initial : Q = I
  // C = 0.0;
  // B = 1.0;

  // std::pair<std::vector<double>, double> retSimplex;

  // int i = 0;
  // double residu = 1.0;
  // while(i < max_iter && residu > 1e-8) {
  //   Aprev = A;
  //   Bprev = B;
  //   Cprev = C;
  //   //==================================================================================
  //   // Calcul de Q^(-1/2)
  //   SMetric3 Q(1.0);
  //   Q.set_m11(A);
  //   Q.set_m21(C);
  //   Q.set_m22(B);
  //   fullMatrix<double> V(3, 3);
  //   fullVector<double> S(3);
  //   Q.eig(V, S, false);

  //   // S(0) = pow(S(0), -0.5);
  //   // S(1) = pow(S(1), -0.5);
  //   // S(2) = pow(S(2), -0.5);

  //   v0 = SVector3(V(0, 0), V(0, 1), V(0, 2)); // Attention c'est peut-être la transposée
  //   v1 = SVector3(V(1, 0), V(1, 1), V(1, 2));
  //   v2 = SVector3(V(2, 0), V(2, 1), V(2, 2));

  //   Q = SMetric3(pow(S(0), -0.5), pow(S(1), -0.5), pow(S(2), -0.5), v0, v1, v2);

  //   // double aTest[4] = {A, C, C, B};
  //   // dsyev(&jobz, &uplo, &n, aTest, &n, w, work, &lwork, &info);

  //   // l1 = w[0];  v11 = aTest[0];  v21 = aTest[2];
  //   // l2 = w[1];  v12 = aTest[1];  v22 = aTest[3];

  //   // Q11 = v11*v11*pow(l1, -0.5) + v21*v21*pow(l2, -0.5);
  //   // Q12 = v12*v11*pow(l1, -0.5) + v21*v22*pow(l2, -0.5);
  //   // Q21 = v11*v12*pow(l1, -0.5) + v21*v22*pow(l2, -0.5);
  //   // Q22 = v12*v12*pow(l1, -0.5) + v22*v22*pow(l2, -0.5);
  //   //==================================================================================
  //   for(int iPhi = 0; iPhi < nPhi; ++iPhi) {
  //     xi = cos(phi[iPhi]);
  //     yi = sin(phi[iPhi]);

  //     xj = Q(0, 0) * xi + Q(0, 1) * yi;
  //     yj = Q(1, 0) * xi + Q(1, 1) * yi;

  //     // err = erreur[0] * xj*xj + erreur[1] * xj*yj + erreur[2] * yj*yj;
  //     err = 0.0;
  //     int indice = 0;
  //     for(int iDeg = 0; iDeg < dimRecovery; ++iDeg) {
  //       if((expX[iDeg] + expY[iDeg]) == degreeSolution + 1) {
  //         err += erreur[indice] * pow(xj, expX[iDeg]) * pow(yj, expY[iDeg]);
  //         ++indice;
  //       }
  //     }

  //     xj = xi / fmax(pow(fabs(err), 1.0 / (degreeSolution + 1)), 1e-20);
  //     yj = yi / fmax(pow(fabs(err), 1.0 / (degreeSolution + 1)), 1e-20);

  //     // ==============================================================================
  //     // Contraintes lineaires
  //     AConstraint[iPhi * nSimplex] = xj * xj;
  //     AConstraint[iPhi * nSimplex + 1] = -xj * xj;
  //     AConstraint[iPhi * nSimplex + 2] = 2.0 * xj * yj;
  //     AConstraint[iPhi * nSimplex + 3] = -2.0 * xj * yj;
  //     AConstraint[iPhi * nSimplex + 4] = yj * yj;
  //     AConstraint[iPhi * nSimplex + 5] = -yj * yj;

  //     normeXj = xj * xj + yj * yj;

  //     bConstraint[iPhi] = -(-normeXj * log(normeXj));
  //   }

  //   retSimplex = simplex(nSimplex, mSimplex, AConstraint, bConstraint, cObjective, 0.0);

  //   if(isinf(retSimplex.second)) {
  //     if(retSimplex.first[0] == -1)
  //       printf("Objective function unbounded!\n");
  //     else if(retSimplex.first[0] == -2)
  //       printf("Linear program infeasible!\n");
  //   } else {
  //     // printf("Solution: (");
  //     // for (int i=0;i<nSimplex+mSimplex;i++) printf("%lf%s", retSimplex.first[i], (i < nSimplex +
  //     // mSimplex - 1) ? ", " : ")\n"); for (int i=0;i<nSimplex;i++) printf("%lf%s",
  //     // retSimplex.first[i], (i < nSimplex - 1) ? ", " : ")\n"); printf("Optimal objective value:
  //     // %lf\n", retSimplex.second);
  //     L1 = retSimplex.first[0] - retSimplex.first[1];
  //     L2 = retSimplex.first[2] - retSimplex.first[3];
  //     L3 = retSimplex.first[4] - retSimplex.first[5];
  //   }

  //   //==================================================================================
  //   // Calcul de Q^(1/2)
  //   // Q11 = v11*v11*pow(l1, 0.5) + v21*v21*pow(l2, 0.5);
  //   // Q12 = v12*v11*pow(l1, 0.5) + v21*v22*pow(l2, 0.5);
  //   // Q21 = v11*v12*pow(l1, 0.5) + v21*v22*pow(l2, 0.5);
  //   // Q22 = v12*v12*pow(l1, 0.5) + v22*v22*pow(l2, 0.5);
  //   Q = SMetric3(pow(S(0), 0.5), pow(S(1), 0.5), pow(S(2), 0.5), v0, v1, v2);
  //   Q11 = Q(0, 0);
  //   Q12 = Q(0, 1);
  //   Q21 = Q(1, 0);
  //   Q22 = Q(1, 1);
  //   //==================================================================================
  //   // Calcul de expm(L)

  //   SMetric3 L(1.0);
  //   L.set_m11(L1);
  //   L.set_m21(L2);
  //   L.set_m22(L3);
  //   // fullMatrix<double> V(3, 3);
  //   // fullVector<double> S(3);
  //   L.eig(V, S, false);

  //   // S(0) = pow(S(0), -0.5);
  //   // S(1) = pow(S(1), -0.5);
  //   // S(2) = pow(S(2), -0.5);

  //   v0 = SVector3(V(0, 0), V(0, 1), V(0, 2)); // Attention c'est peut-être la transposée
  //   v1 = SVector3(V(1, 0), V(1, 1), V(1, 2));
  //   v2 = SVector3(V(2, 0), V(2, 1), V(2, 2));

  //   L = SMetric3(exp(S(0)), exp(S(1)), exp(S(2)), v0, v1, v2);

  //   expL11 = L(0, 0);
  //   expL12 = L(0, 1);
  //   expL21 = L(1, 0);
  //   expL22 = L(1, 1);

  //   // double aTest2[4] = {L1, L2, L2, L3};
  //   // dsyev(&jobz, &uplo, &n, aTest2, &n, w, work, &lwork, &info);

  //   // l1 = w[0];  v11 = aTest2[0];  v21 = aTest2[2];
  //   // l2 = w[1];  v12 = aTest2[1];  v22 = aTest2[3];

  //   // expL11 = v11*v11*exp(l1) + v21*v21*exp(l2);
  //   // expL12 = v12*v11*exp(l1) + v21*v22*exp(l2);
  //   // expL21 = v11*v12*exp(l1) + v21*v22*exp(l2);
  //   // expL22 = v12*v12*exp(l1) + v22*v22*exp(l2);
  //   //==================================================================================
  //   // Calcul de Q = Q^(1/2) * expm(L) * Q^(1/2)
  //   A = Q11 * (Q11 * expL11 + Q12 * expL21) + Q21 * (Q11 * expL12 + Q12 * expL22);
  //   C = Q12 * (Q11 * expL11 + Q12 * expL21) + Q22 * (Q11 * expL12 + Q12 * expL22);
  //   B = Q12 * (Q21 * expL11 + Q22 * expL21) + Q22 * (Q21 * expL12 + Q22 * expL22);

  //   residu =
  //     sqrt((A - Aprev) * (A - Aprev) + 2 * (C - Cprev) * (C - Cprev) + (B - Bprev) * (B - Bprev));
  //   // cout<<"======== ITERATION "<<i+1<<" - RESIDU = "<<residu<<" ==========="<<endl;
  //   ++i;
  // }
}

static void contractDerivativesCurved(double xdot[2], double xdotdot[2], const double hess[2][2], 
  const double cijk[2][2][2], double &C, double &H)
{
  H = 0;
  C = 0;
  for(int i = 0; i < 2; ++i){
    for(int j = 0; j < 2; ++j){
      H += hess[i][j] * xdot[i] * xdotdot[j];
      for(int k = 0; k < 2; ++k){
        C += cijk[i][j][k] * xdot[i] * xdot[j] * xdot[k];
      }
    }
  }
}

static double bissectionFun(const double e, const double sBar, double xdotdot[2], const double v[2], const double w[2], 
  const double k, const double hess[2][2], const double cijk[2][2][2])
{
  double I = 0.;
  double ds = 1e-2;
  double s = 0.;
  double xdot[2]; 
  double C, H, f, fprev = 0.;

  while(s < sBar){
    xdot[0] = v[0] + k*w[0]*s;
    xdot[1] = v[1] + k*w[1]*s;
    contractDerivativesCurved(xdot, xdotdot, hess, cijk, C, H);
    f = (sBar-s)*(sBar-s) * fabs(C + 3.*H);
    I += (f+fprev)/2. * ds;
    fprev = f;
    s += ds;
  }

  return I - e;
}

// Turn around the vertex with a quadratic parametrization, construct error curve
// and fit an ellipse inside (brute-force for now).
void computeWorstMetric(int nTheta, int nIncr, double e, double *x, double cG, double sG, 
  double &h1Opt, double &h2Opt, feRecovery *rec, SMetric3 &M, double lMin, double lMax){

  double theta, k, C, H, s;
  double v[2], w[2], g1[2], g2[2], xdot[2], xdotdot[2];

  std::vector<double> xE(4*nTheta, 0.);
  std::vector<double> yE(4*nTheta, 0.);

  g1[0] = cG;
  g1[1] = sG;
  g2[0] = -sG;
  g2[1] = cG;

  double c1 = fx(rec, x);
  double c2 = fy(rec, x);

  double c11 = fxx(rec, x);
  double c12 = (fxy(rec, x) + fyx(rec, x)) / 2.;
  double c22 = fyy(rec, x);

  double c111 = fxxx(rec, x);
  double c222 = fyyy(rec, x);
  double c112 = (fxxy(rec, x) + fxyx(rec, x) + fyxx(rec, x)) / 3.;
  double c122 = (fxyy(rec, x) + fyyx(rec, x) + fyxy(rec, x)) / 3.;

  // double grad[2] = {c1, c2};
  double hess[2][2] = {{c11, c12}, {c12, c22}};
  double cijk[2][2][2] = {{{c111, c112}, {c112, c122}}, {{c112, c122}, {c122, c222}}};
  
  double k1 = (-c2 * c2 * c11 + 2.0 * c1 * c2 * c12 - c1 * c1 * c22) / (pow(c1 * c1 + c2 * c2, 3.0 / 2.0));
  double k2 = (c1 * c2 * (c22 - c11) + (c1 * c1 - c2 * c2) * c12) / (pow(c1 * c1 + c2 * c2, 3.0 / 2.0));

  for(int iDir = 0; iDir < 4; ++iDir){
    for(int i = 0; i < nTheta; ++i){
      theta = i * M_PI/2. / nTheta;

      if(iDir == 0 || iDir == 1){
        v[0] = g1[0] * cos(theta) + g2[0] * sin(theta);
        v[1] = g1[1] * cos(theta) + g2[1] * sin(theta);
      } else{
        v[0] = -g2[0] * cos(theta) + g1[0] * sin(theta);
        v[1] = -g2[1] * cos(theta) + g1[1] * sin(theta);
      }

      if(iDir == 1 || iDir == 3){
        v[0] = -v[0];
        v[1] = -v[1];
      }

      // w = cross(ez, v)
      w[0] = -v[1];
      w[1] = v[0];

      switch(iDir){
        case 0 : k = + k1 * cos(theta) + k2*sin(theta); break;
        case 1 : k = - k1 * cos(theta) - k2*sin(theta); break;
        case 2 : k = - k2 * cos(theta) + k1*sin(theta); break;
        case 3 : k = + k2 * cos(theta) - k1*sin(theta); break;
      }

      xdot[0] = v[0];
      xdot[1] = v[1];
      xdotdot[0] = k*w[0];
      xdotdot[1] = k*w[1];

      contractDerivativesCurved(xdot, xdotdot, hess, cijk, C, H);

      s = pow( (6. * e / fabs(C + 3.*H)) , 1./3. );
      xE[iDir*nTheta+i] = x[0] + v[0]*s + k*w[0]*s*s/2.;
      yE[iDir*nTheta+i] = x[1] + v[1]*s + k*w[1]*s*s/2.;

      // // // // ////////////////////////////////////////////////

      // // Bissection :

      // double sA = 0.01;
      // double sB = 2.0*s;

      // double fA = bissectionFun(e, sA, xdotdot, v, w, k, hess, cijk);
      // double fB = bissectionFun(e, sB, xdotdot, v, w, k, hess, cijk);
      // if(fA*fB > 0){
      //   feWarning("Bad first interval for bissection with starting s = %1.4e", s);
      //   feWarning("Initial guesses : sA = %1.4e - fA = %+-1.4e - sB = %1.4e fB = %+-1.4e", 
      //     sA, fA, sB, fB);
      //   exit(-1);
      //   break;
      // }

      // double sMid, fMid, delta = (sB-sA)/2.0, tol = 1e-2;
      // int iter = 0;
      // while(iter < 50 && fabs(delta) > tol){
      //   iter++;
      //   delta = (sB-sA)/2.0;
      //   sMid = sA + delta;
      //   fMid = bissectionFun(e, sMid, xdotdot, v, w, k, hess, cijk);

      //   if(fA*fMid > 0){
      //     sA = sMid;
      //     fA = fMid;
      //   } else{
      //     sB = sMid;
      //     fB = fMid;
      //   }
      // }

      // // feInfo("Bissection converged in %d iterations : s = %1.4e changed to sMid = %1.4e", iter, s, sMid);

      // xE[iDir*nTheta+i] = x + v[0]*sMid + k*w[0]*sMid*sMid/2.;
      // yE[iDir*nTheta+i] = y + v[1]*sMid + k*w[1]*sMid*sMid/2.;

    }
  }

  // Fit an ellipse with the prescribed principal directions
  double h1, h2, H1 = h1Opt, H2 = h2Opt;
  double area, areaMax = 0., d, r;
  double thetaG = atan2(sG,cG);
  bool OK;
  for(int i = 0; i < nIncr; ++i){
    h1 = i * H1 / (double) nIncr;
    for(int j = 0; j < nIncr; ++j){
      h2 = j * H2 / (double) nIncr;
      area = h1*h2;
      if(area > areaMax){
        // Check the constraints (error curve)
        OK = true;
        for(int k = 0; k < 4*nTheta; ++k){
          d = sqrt( (xE[k] - x[0])*(xE[k] - x[0]) + (yE[k] - x[1])*(yE[k] - x[1]) );
          theta = atan2(yE[k] - x[1], xE[k] - x[0]);
          r = h1*h2/sqrt( (h1*sin(theta-thetaG))*(h1*sin(theta-thetaG)) + (h2*cos(theta-thetaG))*(h2*cos(theta-thetaG)) );
          if(r > d){
            OK = false;
            break;
          }
        }
        if(OK){
          areaMax = area;
          h1Opt = h1;
          h2Opt = h2;
        }
      }
    }
  }

  // double AA, BB, CC;
  // double *xNew = (double *) malloc(sizeof(double) * 4*nTheta);
  // double *yNew = (double *) malloc(sizeof(double) * 4*nTheta);
  // metricHechtKuate(4*nTheta, xE.data(), yE.data(), AA, BB, CC, 1e-5, xNew, yNew);
  // free(xNew);
  // free(yNew);

  // SMetric3 m;
  // m.set_m11(AA);
  // m.set_m21(CC/2.0);
  // m.set_m22(BB);
  // m.set_m33(1.0);

  // fullMatrix<double> V(3, 3);
  // fullVector<double> S(3);
  // m.eig(V, S, false);
  // SVector3 v0(V(0, 0), V(0, 1), V(0, 2));
  // SVector3 v1(V(1, 0), V(1, 1), V(1, 2));
  // SVector3 v2(V(2, 0), V(2, 1), V(2, 2));
  // m = SMetric3(fmin(lMax, fmax(lMin, S(0))), fmin(lMax, fmax(lMin, S(1))),
  //                fmin(lMax, fmax(lMin, S(2))), v0, v1, v2);

  // M = m;

  // double tol = 1e-1;
  // if(fabs(x - 0.158414) < tol && fabs(y - 0.453362) < tol){
  //   FILE *myf = fopen("errorCurveBissection.pos", "w");
  //   fprintf(myf, "View \" \"{\n");

  //   // Print the error level
  //   for(int i = 0; i < xE.size(); ++i) {
  //     fprintf(myf, "SP(%.16g,%.16g,%.16g){%f};\n", xE[i], yE[i], 0., 1.);
  //   }
  //   // Print the kuate ellipse
  //   int nt = 30;
  //   std::vector<double> xP(nt, 0.);
  //   std::vector<double> yP(nt, 0.);
  //   getEllipsePoints(M(0, 0), 2.0 * M(0, 1), M(1, 1), x, y, xP, yP);
  //   for(int j = 0; j < nt; ++j) {
  //     if(j != nt - 1) {
  //       fprintf(myf, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
  //               0., xP[j + 1], yP[j + 1], 0., 1, 1);
  //     } else {
  //       fprintf(myf, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
  //               0., xP[0], yP[0], 0., 1, 1);
  //     }
  //   }

  //   fprintf(myf, "};");
  //   fclose(myf);
  // }
}

SVector3 E_X(1., 0., 0.);
SVector3 E_Y(0., 1., 0.);
SVector3 E_Z(0., 0., 1.);
double tolerance_gradation = 1e-1; // Metrics are modified if || M - M_grown || < tol
double tolerance_diagonal = 1e-4;  // If fabs(non-diag) is below this tol, matrix is diag

static inline double matNorm2(const SMetric3 &m){
  double sqr = 0;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      sqr += m(i, j)*m(i, j);
    }
  }
  return sqrt(sqr);
}

static inline double matNorm2(const SMetric3 &m1, const SMetric3 &m2){
  double sqr = 0;
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      sqr += (m1(i, j) - m2(i, j))*(m1(i, j) - m2(i, j));
    }
  }
  return sqrt(sqr);
}

static inline SMetric3 intersectionReductionSimultaneeExplicite(const SMetric3 &m1, const SMetric3 &m2)
{
  double a1 = m1(0,0), b1 = m1(0,1), c1 = m1(1,1);
  double a2 = m2(0,0), b2 = m2(0,1), c2 = m2(1,1);
  double v00, v01, vTInv00, vTInv01, vTInv10, vTInv11;
  // Metrics are identical
  double relDiff = matNorm2(m1,m2)/matNorm2(m1);
  if(matNorm2(m1,m2)/matNorm2(m1) < 1e-2){
      // std::cout<<"Identical matrix : returning m1"<<std::endl;
    return m2;
  }

  if(fabs(b1) < tolerance_diagonal && fabs(b2) < tolerance_diagonal){ 

    // Both matrices are diagonal
    return SMetric3(fmax(a1,a2), fmax(c1,c2), 1., E_X, E_Y, E_Z);

  } else if(fabs(b1) < tolerance_diagonal){ 

    // m1 is diagonal
    v00 = (a1*c2 + a2*c1 - sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 + 4*a1*b2*b2*c1 + a2*a2*c1*c1))/(2*a1*b2) - c2/b2;
    v01 = (a1*c2 + a2*c1 + sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 + 4*a1*b2*b2*c1 + a2*a2*c1*c1))/(2*a1*b2) - c2/b2;
    SVector3 v0(v00, 1., 0.); //v0.normalize();
    SVector3 v1(v01, 1., 0.); //v1.normalize();
    double l0 = fmax(dot(v0, m1, v0), dot(v0, m2, v0));
    double l1 = fmax(dot(v1, m1, v1), dot(v1, m2, v1));
    vTInv00 = -(a1*b2)/sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 + 4*a1*b2*b2*c1 + a2*a2*c1*c1);
    vTInv01 = (sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 + 4*a1*b2*b2*c1 + a2*a2*c1*c1) - a1*c2 + a2*c1)/(2*sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 + 4*a1*b2*b2*c1 + a2*a2*c1*c1));
    vTInv10 =  (a1*b2)/sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 + 4*a1*b2*b2*c1 + a2*a2*c1*c1);
    vTInv11 = (sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 + 4*a1*b2*b2*c1 + a2*a2*c1*c1) + a1*c2 - a2*c1)/(2*sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 + 4*a1*b2*b2*c1 + a2*a2*c1*c1));
    if(isnan(vTInv00) || isnan(vTInv01) || isnan(vTInv10) || isnan(vTInv11) || isinf(vTInv00) || isinf(vTInv01) || isinf(vTInv10) || isinf(vTInv11)){
      std::cout<<relDiff<<std::endl;
      m1.print("1 - m1");
      m2.print("1 - m2");
      std::cout<<"NAN alert : "<<vTInv00<<" "<<vTInv01<<" "<<vTInv10<<" "<<vTInv11<<std::endl;
      exit(-1);
    }
    v0 = SVector3(vTInv00, vTInv01, 0.); //v0.normalize();
    v1 = SVector3(vTInv10, vTInv11, 0.); //v1.normalize();
    return SMetric3(l0, l1, 1., v0, v1, E_Z);

  } else if(fabs(b2) < tolerance_diagonal){ 

    // m2 is diagonal
    v00 = (a1*c2)/(a2*b1) - (a1*c2 + a2*c1 + sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 + a2*a2*c1*c1 + 4*a2*b1*b1*c2))/(2*a2*b1);
    v01 = (a1*c2)/(a2*b1) - (a1*c2 + a2*c1 - sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 + a2*a2*c1*c1 + 4*a2*b1*b1*c2))/(2*a2*b1);
    SVector3 v0(v00, 1., 0.); //v0.normalize();
    SVector3 v1(v01, 1., 0.); //v1.normalize();
    double l0 = fmax(dot(v0, m1, v0), dot(v0, m2, v0));
    double l1 = fmax(dot(v1, m1, v1), dot(v1, m2, v1));
    vTInv00 = -(a2*b1)/sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 + a2*a2*c1*c1 + 4*a2*b1*b1*c2);
    vTInv01 = (sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 + a2*a2*c1*c1 + 4*a2*b1*b1*c2) + a1*c2 - a2*c1)/(2*sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 + a2*a2*c1*c1 + 4*a2*b1*b1*c2));
    vTInv10 =  (a2*b1)/sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 + a2*a2*c1*c1 + 4*a2*b1*b1*c2);
    vTInv11 = (sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 + a2*a2*c1*c1 + 4*a2*b1*b1*c2) - a1*c2 + a2*c1)/(2*sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 + a2*a2*c1*c1 + 4*a2*b1*b1*c2));
    if(isnan(vTInv00) || isnan(vTInv01) || isnan(vTInv10) || isnan(vTInv11) || isinf(vTInv00) || isinf(vTInv01) || isinf(vTInv10) || isinf(vTInv11)){
      std::cout<<relDiff<<std::endl;
      m1.print("2 - m1");
      m2.print("2 - m2");
      std::cout<<"NAN alert : "<<vTInv00<<" "<<vTInv01<<" "<<vTInv10<<" "<<vTInv11<<std::endl;
      exit(-1);
    }
    v0 = SVector3(vTInv00, vTInv01, 0.); //v0.normalize();
    v1 = SVector3(vTInv10, vTInv11, 0.); //v1.normalize();
    return SMetric3(l0, l1, 1., v0, v1, E_Z);

  } else{

    // Check if metrics are multiple of one another
    if( fabs(a2/a1 - c2/c1) < 1e-3 && fabs(a2/a1 - b2/b1) < 1e-3 ){
      return (a2/a1 <= 1) ? m1 : m2;
    } else{
      v00 = (sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 - 4*a1*b1*b2*c2 + 4*a1*b2*b2*c1 + a2*a2*c1*c1 + 4*a2*b1*b1*c2 - 4*a2*b1*b2*c1) + a1*c2 + a2*c1 - 2*b1*b2)/(2*(a1*b2 - a2*b1)) - (a1*c2 - b1*b2)/(a1*b2 - a2*b1);
      v01 = - (a1*c2 - b1*b2)/(a1*b2 - a2*b1) - (sqrt(a1*a1*c2*c2 - 2.*a1*a2*c1*c2 - 4.*a1*b1*b2*c2 + 4.*a1*b2*b2*c1 + a2*a2*c1*c1 + 4.*a2*b1*b1*c2 - 4.*a2*b1*b2*c1) - a1*c2 - a2*c1 + 2.*b1*b2)/(2.*(a1*b2 - a2*b1));
      SVector3 v0(v00, 1., 0.); //v0.normalize();
      SVector3 v1(v01, 1., 0.); //v1.normalize();
      double l0 = fmax(dot(v0, m1, v0), dot(v0, m2, v0));
      double l1 = fmax(dot(v1, m1, v1), dot(v1, m2, v1));
      vTInv00 =  (a1*b2 - a2*b1)/sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 - 4*a1*b1*b2*c2 + 4*a1*b2*b2*c1 + a2*a2*c1*c1 + 4*a2*b1*b1*c2 - 4*a2*b1*b2*c1);
      vTInv01 = (sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 - 4*a1*b1*b2*c2 + 4*a1*b2*b2*c1 + a2*a2*c1*c1 + 4*a2*b1*b1*c2 - 4*a2*b1*b2*c1) + a1*c2 - a2*c1)/(2*sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 - 4*a1*b1*b2*c2 + 4*a1*b2*b2*c1 + a2*a2*c1*c1 + 4*a2*b1*b1*c2 - 4*a2*b1*b2*c1));
      vTInv10 = -(a1*b2 - a2*b1)/sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 - 4*a1*b1*b2*c2 + 4*a1*b2*b2*c1 + a2*a2*c1*c1 + 4*a2*b1*b1*c2 - 4*a2*b1*b2*c1);
      vTInv11 = (sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 - 4*a1*b1*b2*c2 + 4*a1*b2*b2*c1 + a2*a2*c1*c1 + 4*a2*b1*b1*c2 - 4*a2*b1*b2*c1) - a1*c2 + a2*c1)/(2*sqrt(a1*a1*c2*c2 - 2*a1*a2*c1*c2 - 4*a1*b1*b2*c2 + 4*a1*b2*b2*c1 + a2*a2*c1*c1 + 4*a2*b1*b1*c2 - 4*a2*b1*b2*c1));  
      if(isnan(vTInv00) || isnan(vTInv01) || isnan(vTInv10) || isnan(vTInv11) || isinf(vTInv00) || isinf(vTInv01) || isinf(vTInv10) || isinf(vTInv11)){
        std::cout<<relDiff<<std::endl;
        m1.print("3 - m1");
        m2.print("3 - m2");
        std::cout<<"NAN alert : "<<vTInv00<<" "<<vTInv01<<" "<<vTInv10<<" "<<vTInv11<<std::endl;
        exit(-1);
      }
      v0 = SVector3(vTInv00, vTInv01, 0.); //v0.normalize();
      v1 = SVector3(vTInv10, vTInv11, 0.); //v1.normalize();
      return SMetric3(l0, l1, 1., v0, v1, E_Z);
    }
  }
}

static inline SMetric3 gradationAlauzet(double grad, const SMetric3 &m, double x1, double y1, double x2, double y2){
  SVector3 pq(x2-x1, y2-y1, 0.);
  double eta = (1. + sqrt(dot(pq,m,pq)) * log(grad));
  eta = 1./(eta*eta);
  SMetric3 tmp(m);
  tmp *= eta;
  tmp.set_m33(1.0);
  return tmp;
}

void gradationMetriques(double gradation, int nmax, std::vector<double> &coord1, std::map<int, SMetric3> &metricsOnGmshModel)
{
#if defined(HAVE_GMSH)
  std::vector<int> elementTypes;
  std::vector<std::vector<std::size_t> > elementTags;
  std::vector<std::vector<std::size_t> > elemNodeTags;
  gmsh::model::mesh::getElements(elementTypes, elementTags, elemNodeTags, 2);

  // std::vector<std::pair<size_t,size_t>> edges;
  std::set<std::pair<size_t,size_t>, gmshEdgeLessThan> edges;

  std::vector<std::size_t> nodeTags;
  std::vector<double> coord;
  std::vector<double> nodecoord1;
  std::vector<double> nodecoord2;
  std::vector<double> nodecoord3;
  std::vector<double> parametricCoord;
  std::vector<double> nodeparametricCoord;
  gmsh::model::mesh::getNodesByElementType(elementTypes[0], nodeTags,
                                          coord,
                                          parametricCoord);

  for(auto val : elementTypes)
    feInfo("Elements classified on dim 2 are of type %d", val);
  feInfo("There are %d elements of dim 2 and type [0]", elementTags[0].size());
  feInfo("There are %d nodes total in elem [0]", nodeTags.size());

  int nodedim, nodetag;

  FILE* fff = fopen("CHECK.pos","w");
  fprintf(fff, "View \" CHECK \"{\n");

  size_t maxtag = -1;
  std::pair<std::set<std::pair<size_t,size_t>, gmshEdgeLessThan>::iterator, bool> ret;

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
  feInfo("done here");

  fprintf(fff, "};");
  fclose(fff);

  int nt = 30;
  std::vector<double> xP(nt, 0.);
  std::vector<double> yP(nt, 0.);

  SMetric3 mp, mq, mpq, mqp;
  double xp, yp, xq, yq;
  bool correction = true;
  int iter = 0;
  while(correction && iter < nmax){

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

      gmsh::model::mesh::getNode(n    , nodecoord1, nodeparametricCoord, nodedim, nodetag);
      gmsh::model::mesh::getNode(neigh, nodecoord2, nodeparametricCoord, nodedim, nodetag);

      xp = nodecoord1[0];
      yp = nodecoord1[1];
      xq = nodecoord2[0];
      yq = nodecoord2[1];      

      // fprintf(fff, "SL(%.16g,%.16g,0.,%.16g,%.16g,0.){%u, %u};\n", xp, yp, xq, yq, 1, 1);
      // fprintf(fff, "SP(%.16g,%.16g,0.){%u};\n", xp, yp, n);
      // Span mp to neighbouring vertex
      mpq = gradationAlauzet(gradation, mp, xp, yp, xq, yq);
      if(fabs(mpq(2,2) - 1.0) > 1e-4){
        std::cout<<"NONZERO z entry in mpq";
          exit(-1);
      }
      mpq = intersectionReductionSimultaneeExplicite(mq, mpq);
      if(fabs(mpq(2,2) - 1.0) > 1e-4){
        std::cout<<"NONZERO z entry in mpq apres gradation";
          exit(-1);
      }

      // mpq = intersection_reductionSimultanee(mq, mpq);
      // metricsOnGmshModel[neigh] = mpq;

      if(matNorm2(metricsOnGmshModel[neigh], mpq)/matNorm2(mpq) > tolerance_gradation){
        metricsOnGmshModel[neigh] = mpq;
        correction = true;
      }
      // Span mq to n
      mqp = gradationAlauzet(gradation, mq, xq, yq, xp, yp);
      if(fabs(mpq(2,2) - 1.0) > 1e-4){
        std::cout<<"NONZERO z entry in mqp";
          exit(-1);
      }
      mqp = intersectionReductionSimultaneeExplicite(mp, mqp);
      if(fabs(mpq(2,2) - 1.0) > 1e-4){
        std::cout<<"NONZERO z entry in mqp apres gradation";
          exit(-1);
      }

      // mqp = intersection_reductionSimultanee(mp, mqp);
      // metricsOnGmshModel[n] = mqp;

      if(matNorm2(metricsOnGmshModel[n], mqp)/matNorm2(mqp) > tolerance_gradation){
        metricsOnGmshModel[n] = mqp;
        correction = true;
      }

      // mp = gradationAlauzet(gradation, mref, X0, Y0, xp, yp);

      // double factor = 1000.;
      // getEllipsePoints(factor * mp(0, 0), factor * 2.0 * mp(0, 1), factor * mp(1, 1), xp, yp, xP, yP);
      // for(int j = 0; j < nt; ++j) {
      //   if(j != nt - 1) {
      //     fprintf(fff, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
      //             0., xP[j + 1], yP[j + 1], 0., 1, 1);
      //   } else {
      //     fprintf(fff, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
      //             0., xP[0], yP[0], 0., 1, 1);
      //   }
      // }
    }

    feInfo("Passe %d", iter);
    // correction = false;

    // fprintf(fff, "};");
    // fclose(fff);

  }
  feInfo("%d coord size", coord.size());
  feInfo("%d edges size", edges.size());
#endif
}

static Eigen::VectorXd CMINUS(7);
static Eigen::VectorXd CPLUS(7);
// static Eigen::VectorXd CMINUS(8);
// static Eigen::VectorXd CPLUS(8);

static double solveErrorFunction(double k, double *v, double *w, double H[2][2], double C[2][2][2], std::set<double> &roots)
{
  double a1 = k*w[0];
  double a2 = k*w[1];
  // double a1 = 0.0;
  // double a2 = 0.0;
  double v1 = v[0];
  double v2 = v[1];

  double H11 = H[0][0];
  double H12 = H[0][1];
  double H21 = H[1][0];
  double H22 = H[1][1];

  double C111 =  C[0][0][0];
  double C112 = (C[0][0][1] + C[0][1][0] + C[1][0][0])/3.;
  double C122 = (C[0][1][1] + C[1][0][1] + C[1][1][0])/3.;
  double C222 =  C[1][1][1];

  double cs6 = (C111*a1*a1*a1 + 3.*C112*a1*a1*a2 + 3.*C122*a1*a2*a2 + C222*a2*a2*a2)/120.;
  double cs5 = C111*a1*a1*v1/20. + C112*a1*a1*v2/20. + C122*a2*a2*v1/20. + C222*a2*a2*v2/20. + C112*a1*a2*v1/10. + C122*a1*a2*v2/10.;
  double cs4 = C111*a1*v1*v1/8.  + C112*a2*v1*v1/8.  + C122*a1*v2*v2/8.  + C222*a2*v2*v2/8.  + C112*a1*v1*v2/4.  + C122*a2*v1*v2/4.;
  double cs3 = C111*v1*v1*v1/6.  + C112*v1*v1*v2/2.  + C122*v1*v2*v2/2.  + C222*v2*v2*v2/6. ;
  cs4 += (H11*a1*a1 + H12*a1*a2 + H21*a1*a2 + H22*a2*a2)/8.;
  cs3 += (H11*a1*v1 + H12*a2*v1 + H21*a1*v2 + H22*a2*v2)/2.;

  // Ajout du développement limité au premier ordre de Hij = Hij(x0) + s*Cijk*vk
  cs5 += (a1*a1/20.) * (C111*v1 + C112*v2);
  cs4 += (v1*a1/8. ) * (C111*v1 + C112*v2);

  cs5 += (a1*a2/20.) * (C112*v1 + C122*v2); // C121*v1 + C122*v2
  cs4 += (v1*a2/8. ) * (C112*v1 + C122*v2);

  cs5 += (a1*a2/20.) * (C112*v1 + C122*v2); // C211*v1 + C212*v2
  cs4 += (v2*a1/8. ) * (C112*v1 + C122*v2);

  cs5 += (a2*a2/20.) * (C122*v1 + C222*v2); // C221*v1 + C222*v2
  cs4 += (v2*a2/8. ) * (C122*v1 + C222*v2);

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
  for (auto it = roots.begin(); it != roots.end(); it++)
  {
    s = fmin(s, *it);
  }

  roots = RootFinder::solvePolynomial(CPLUS,  0., INFINITY, 1e-8);
  for (auto it = roots.begin(); it != roots.end(); it++)
  {
      s = fmin(s, *it);
  }

  // feInfo("lowest root is %1.16e", s);
  return s;
}

#if defined(HAVE_SOPLEX)
static void computeLvl1(double K1, double K2, double *g1, double *g2, 
  double Hij[2][2], double Cijk[2][2][2], int nPhi, linearProblem &myLP)
{
  double v1[2], v2[2], k1, k2, theta, k, v[2], w[2], s, x0[2] = {0., 0.};

  std::set<double> roots;

  for(int iDir = 0; iDir < 4; ++iDir){
    switch(iDir){
      case 0 : 
        v1[0] = g1[0];
        v1[1] = g1[1];
        v2[0] = g2[0];
        v2[1] = g2[1];
        k1 = K1;
        k2 = K2;
        break;
      case 1 :
        v1[0] = -g1[0];
        v1[1] = -g1[1];
        v2[0] = -g2[0];
        v2[1] = -g2[1];
        k1 = -K1;
        k2 = -K2;
        break;
      case 2 :
        v1[0] = -g2[0];
        v1[1] = -g2[1];
        v2[0] =  g1[0];
        v2[1] =  g1[1];
        k1 = -K2;
        k2 =  K1;
        break;
      case 3 :
        v1[0] =  g2[0];
        v1[1] =  g2[1];
        v2[0] = -g1[0];
        v2[1] = -g1[1];
        k1 =  K2;
        k2 = -K1;
        break;
    }

    for(int i = 0; i < nPhi; ++i){
      theta = M_PI/2. * (double) i/ (double) nPhi;
      k = k1*cos(theta)*cos(theta) + k2*sin(theta)*sin(theta);
      v[0] = v1[0]*cos(theta) + v2[0]*sin(theta);
      v[1] = v1[1]*cos(theta) + v2[1]*sin(theta);
      w[0] = -v[1];
      w[1] =  v[0];

      s = solveErrorFunction(k, v, w, Hij, Cijk, roots);

      // Shoot :
      myLP.lvl1[iDir*2*nPhi + 2*i + 0] = x0[0] + s*v[0] + s*s/2. * k*w[0];
      myLP.lvl1[iDir*2*nPhi + 2*i + 1] = x0[1] + s*v[1] + s*s/2. * k*w[1];
    }
  }
}

// Not used if we use Eigen's Q.sqrt()
static void powerOneHalf(Eigen::Matrix2d &Q, double Q12[3])
{
  double a = Q(0,0);
  double b = Q(0,1);
  double c = Q(1,1);
  // Symbolic expression from Matlab for Q12 = V * sqrt(Lambda) * V^(-1)
  double f8 = sqrt(a*a - 2.*a*c + 4.*b*b + c*c);
  double rPlus  = sqrt(a + c + f8);
  double rMinus = sqrt(a + c - f8);
  double r2 = sqrt(2.);
  double f7 = r2 * rPlus  * f8;
  double f6 = r2 * rMinus * f8;
  double f5 = r2 * a * rPlus;
  double f4 = r2 * a * rMinus;
  double f3 = r2 * c * rPlus;
  double f2 = r2 * c * rMinus;
  double f1 = r2 * b * (rPlus - rMinus) / (2. * f8);
  // We only store the 3 coefficients a, b, c in Q12
  Q12[0] = (f5 - f2 + f7 - f4 + f2 + f6) / (4. * f8);
  Q12[1] = f1;
  Q12[2] = (f3 - f5 + f7 + f4 - f2 + f6) / (4. * f8);
}

static bool solveLP(linearProblem &myLP, int nPhi, Eigen::Matrix2d &L)
{
  // DSVector row(3);
  // LPRow lprow(3);
  // LPRowSet lpRowSet(4*nPhi,3);
  // lprow.setRhs(infinity);

  double xi, yi, normSquared, lhs;
  for(int i = 0; i < 4*nPhi; ++i){
    xi = myLP.constraints[2*i];
    yi = myLP.constraints[2*i+1];
    normSquared = xi*xi + yi*yi;
    lhs = -normSquared * log(normSquared);
    myLP.row.add(0,    xi*xi);
    myLP.row.add(1, 2.*xi*yi);
    myLP.row.add(2,    yi*yi);
    myLP.lprow.setLhs(lhs);
    myLP.lprow.setRowVector(myLP.row);
    myLP.lprowset.add(myLP.lprow);
    // problem.changeRowReal(i, lprow);
    // problem.changeRowReal(i, LPRow(lhs, row, infinity));
    myLP.row.clear();
  }
  myLP.problem.removeRowRangeReal(0, 4*nPhi-1);
  myLP.problem.addRowsReal(myLP.lprowset);
  myLP.lprowset.clear();
  /* solve LP */
  // SPxSolver::Status stat;
  myLP.stat = myLP.problem.optimize();

  if(myLP.stat == SPxSolver::OPTIMAL) {
    myLP.problem.getPrimal(myLP.prim);
    L(0,0) = myLP.prim[0];
    L(0,1) = myLP.prim[1];
    L(1,0) = myLP.prim[1];
    L(1,1) = myLP.prim[2];
    return true;
  } else {
    L(0,0) = 1.0;
    L(0,1) = 0.0;
    L(1,0) = 0.0;
    L(1,1) = 1.0;
    return false;
  }
}

static bool logSimplexCurved(double *x, double k1, double k2, double g1[2], double g2[2], 
  double Hij[2][2], double Cijk[2][2][2], int maxIter, int nThetaPerQuadrant, 
  double tol, Eigen::Matrix2d &Q, int &numIter, linearProblem &myLP)
{
  int nPhi = nThetaPerQuadrant;
  // int size = 2 * 4 * nPhi; // 2 coordonnees * 4 quadrants * nPhi angles/quadrant
  // std::vector<double> lvl1(size, 0.), constraints(size, 0.);
  computeLvl1(k1, k2, g1, g2, Hij, Cijk, nPhi, myLP);

  Eigen::Matrix2d Qprev, Q12, L, diff;

  Q = Eigen::Matrix2d::Identity();

  double xi, yi;
  bool res;

  for(int iter = 0; iter < maxIter; ++iter){
    Qprev = Q;
    Q12 = Q.sqrt();

    // Apply transformation Q12 to the initial constraints
    for(int i = 0; i < 4*nPhi; ++i){
      xi = myLP.lvl1[2*i];
      yi = myLP.lvl1[2*i+1];
      myLP.constraints[2*i]   = Q12(0,0) * xi + Q12(0,1) * yi;
      myLP.constraints[2*i+1] = Q12(1,0) * xi + Q12(1,1) * yi;
    } 

    // Solve the linear optimization problem for L
    res = solveLP(myLP, nPhi, L);

    if(res){
      // Recover Q from L
      Q = Q12 * L.exp() * Q12;

      diff = Q - Qprev;
      if(diff.norm() < tol){
        numIter = iter;
        break;
      }
    } else{
      // LP solver returned an error : return the last valid Q
      if(iter == 0){
        feWarning("In logSimplexCurved : LP solver returned an error at iteration %d. Returning identity matrix.", iter);
        FILE *f = fopen("lvl1_check.pos", "w");
        fprintf(f, "View\"lvl1\"{\n");
        for(int i = 0; i < 4*nPhi; ++i){
          // feInfo("%f - %f", myLP.lvl1[2*i], myLP.lvl1[2*i+1]);
          fprintf(f, "SP(%g,%g,0){%g,%g,0};", x[0] + myLP.lvl1[2*i], x[1] + myLP.lvl1[2*i+1], 1., 1.);
        }
        int nt = 30;
        std::vector<double> xP(nt, 0.);
        std::vector<double> yP(nt, 0.);
        double factor = 1.;
        getEllipsePoints(factor * Q(0, 0), factor * 2.0 * Q(0, 1), factor * Q(1, 1), x[0], x[1], xP, yP);
        for(int j = 0; j < nt; ++j) {
          if(j != nt - 1) {
            fprintf(f, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
                    0., xP[j + 1], yP[j + 1], 0., 1, 1);
          } else {
            fprintf(f, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
                    0., xP[0], yP[0], 0., 1, 1);
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

bool computeMetricLogSimplexCurved(double *x, double cG, double sG, feRecovery *rec, 
  Eigen::Matrix2d &Q, int maxIter, int nThetaPerQuadrant, double tol, int &numIter, linearProblem &myLP)
{
  double g1[2], g2[2];
  g1[0] = cG;
  g1[1] = sG;
  g2[0] = -sG;
  g2[1] = cG;

  // Derivatives from recoveries
  double c1 = fx(rec, x);
  double c2 = fy(rec, x);

  double c11 = fxx(rec, x);
  double c12 = (fxy(rec, x) + fyx(rec, x)) / 2.;
  double c22 = fyy(rec, x);

  double c111 = fxxx(rec, x);
  double c222 = fyyy(rec, x);
  double c112 = (fxxy(rec, x) + fxyx(rec, x) + fyxx(rec, x)) / 3.;
  double c122 = (fxyy(rec, x) + fyyx(rec, x) + fyxy(rec, x)) / 3.;

  // Analytical derivatives
  // double X = x[0];
  // double Y = x[1];

  // tanh
  // double a = 10.;
  // double b = 1.5;
  // double T = pow(tanh(a*(X/2. - sin(M_PI*b*Y)/4.)),2);

  // double c1 = -(a*(T - 1))/4.;
  // double c2 = (a*b*M_PI*cos(M_PI*b*Y)*(T - 1.))/8.;

  // double c11 = (a*a*tanh(a*(X/2. - sin(M_PI*b*Y)/4.))*(T - 1))/4;
  // double c12 = -(a*a*b*M_PI*tanh(a*(X/2. - sin(M_PI*b*Y)/4.))*cos(M_PI*b*Y)*(T - 1.))/8.;
  // double c22 = (a*a*b*b*M_PI*M_PI*tanh(a*(X/2. - sin(M_PI*b*Y)/4.))*pow(cos(M_PI*b*Y),2)*(T - 1.))/16. - (a*b*b*M_PI*M_PI*sin(M_PI*b*Y)*(T - 1.))/8.;

  // double y = Y;
  // double pi = M_PI;
  // double c111 = - (a*a*a*pow((T - 1),2))/8 - (a*a*a*T*(T - 1))/4;
  // double c112 = (a*a*a*b*pi*cos(pi*b*y)*pow((T - 1),2))/16 + (a*a*a*b*pi*T*cos(pi*b*y)*(T - 1))/8;
  // double c122 = (a*a*b*b*M_PI*M_PI*tanh(a*(X/2 - sin(pi*b*y)/4))*sin(pi*b*y)*(T - 1))/8 - (a*a*a*b*b*M_PI*M_PI*T*pow(cos(pi*b*y),2)*(T - 1))/16 - (a*a*a*b*b*M_PI*M_PI*pow(cos(pi*b*y),2)*pow((T - 1),2))/32;
  // double c222 = (a*a*a*b*b*b*M_PI*M_PI*M_PI*pow(cos(pi*b*y),3)*pow((T - 1),2))/64 - (a*b*b*b*M_PI*M_PI*M_PI*cos(pi*b*y)*(T - 1))/8 + (a*a*a*b*b*b*M_PI*M_PI*M_PI*T*pow(cos(pi*b*y),3)*(T - 1))/32 - (3*a*a*b*b*b*M_PI*M_PI*M_PI*tanh(a*(X/2 - sin(pi*b*y)/4))*cos(pi*b*y)*sin(pi*b*y)*(T - 1))/16;

  // Analytique : r4 = (x^2+y^2)^2 = x^4 + y^4 + 2*x^2 y^2
  // double c1 = 4.*X*X*X + 4. * X*Y*Y;
  // double c2 = 4.*Y*Y*Y + 4. * X*X*Y;

  // double c11 = 12.*X*X + 4.*Y*Y;
  // double c12 = 8.*X*Y;
  // double c22 = 12.*Y*Y + 4.*X*X;

  // double c111 = 24.*X;
  // double c112 = 8.*Y;
  // double c122 = 8.*X;
  // double c222 = 24.*Y;

  // Analytique : x^3 + 20*y^3
  // double c1 =3*X*X;
  // double c2 = 60.*Y*Y;

  // double c11 = 6.0*X;
  // double c12 = 0.;
  // double c22 = 120.*Y;

  // double c111 = 6.0;
  // double c112 = 0.;
  // double c122 = 120.; /// ???
  // double c222 = 0.;

  // Analytique papier IMR : atan(10*(sin(3*pi*y/2) - 2x))
  // double X = x[0];
  // double Y = x[1];
  // double pi = M_PI;
  // double T = pow(20.*X - 10.*sin((3.*pi*Y)/2.),2);
  // double TT = (T + 1.)*(T + 1.);
  // double TTT = (T + 1.)*(T + 1.)*(T + 1.);

  // double R = 800*X - 400*sin((3*pi*Y)/2);
  // double S = cos((3*pi*Y)/2);

  // double c1 = -20/(T + 1);
  // double c2 = (15*pi*S)/(T + 1);
  
  // double c11 = (20*R)/TT;
  // double c12 = -(600*pi*S*(20*X - 10*sin((3*pi*Y)/2)))/TT;
  // double c22 = (450*pi*pi*S*S*(20*X - 10*sin((3*pi*Y)/2)))/TT - (45*pi*pi*sin((3*pi*Y)/2))/(2*(T + 1));
  
  // double c111 = 16000/TT - (40*R*R)/TTT;
  // double c112 = (1200*pi*S*(20*X - 10*sin((3*pi*Y)/2))*R)/TTT - (12000*pi*S)/TT;
  // double c122 = (9000*pi*pi*S*S)/TT + (900*pi*pi*sin((3*pi*Y)/2)*(20*X - 10*sin((3*pi*Y)/2)))/TT - (36000*pi*pi*S*S*T)/TTT;
  // double c222 = (27000*pi*pi*pi*S*S*S*T)/TTT - (135*pi*pi*pi*S)/(4*(T + 1)) - (6750*pi*pi*pi*S*S*S)/TT - (2025*pi*pi*pi*S*sin((3*pi*Y)/2)*(20*X - 10*sin((3*pi*Y)/2)))/TT;

  double Hij[2][2] = {{c11, c12}, {c12, c22}};
  double Cijk[2][2][2] = {{{c111, c112}, {c112, c122}}, {{c112, c122}, {c122, c222}}};
  
  double k1 = (-c2 * c2 * c11 + 2.0 * c1 * c2 * c12 - c1 * c1 * c22) / (pow(c1 * c1 + c2 * c2, 3.0 / 2.0));
  // k1 = -k1;
  double k2 = (c1 * c2 * (c22 - c11) + (c1 * c1 - c2 * c2) * c12) / (pow(c1 * c1 + c2 * c2, 3.0 / 2.0));

  bool res = logSimplexCurved(x, k1, k2, g1, g2, Hij, Cijk, maxIter, nThetaPerQuadrant, tol, Q, numIter, myLP);
  return res;
}
#endif