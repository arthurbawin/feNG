#include "feMetricTools.h"
#include "feSimplex.h"

#if defined(HAVE_GMSH)
#include "gmsh.h"
#endif
#include "fullMatrix.h"
#include "ellipseToolbox.h"

double evaluateFieldFromRecovery(int indexDerivative, feRecovery *rec, double x, double y)
{
  std::vector<double> pos = {x, y, 0.0};
  return rec->evalDerivative(indexDerivative, pos);
}

double f(feRecovery *rec, double x, double y) { return evaluateFieldFromRecovery(0, rec, x, y); }
double fx(feRecovery *rec, double x, double y) { return evaluateFieldFromRecovery(1, rec, x, y); }
double fy(feRecovery *rec, double x, double y) { return evaluateFieldFromRecovery(2, rec, x, y); }
double fxx(feRecovery *rec, double x, double y) { return evaluateFieldFromRecovery(3, rec, x, y); }
double fxy(feRecovery *rec, double x, double y) { return evaluateFieldFromRecovery(4, rec, x, y); }
double fyx(feRecovery *rec, double x, double y) { return evaluateFieldFromRecovery(5, rec, x, y); }
double fyy(feRecovery *rec, double x, double y) { return evaluateFieldFromRecovery(6, rec, x, y); }
double fxxx(feRecovery *rec, double x, double y) { return evaluateFieldFromRecovery(7, rec, x, y); }
double fxxy(feRecovery *rec, double x, double y) { return evaluateFieldFromRecovery(8, rec, x, y); }
double fxyx(feRecovery *rec, double x, double y) { return evaluateFieldFromRecovery(9, rec, x, y); }
double fxyy(feRecovery *rec, double x, double y)
{
  return evaluateFieldFromRecovery(10, rec, x, y);
}
double fyxx(feRecovery *rec, double x, double y)
{
  return evaluateFieldFromRecovery(11, rec, x, y);
}
double fyxy(feRecovery *rec, double x, double y)
{
  return evaluateFieldFromRecovery(12, rec, x, y);
}
double fyyx(feRecovery *rec, double x, double y)
{
  return evaluateFieldFromRecovery(13, rec, x, y);
}
double fyyy(feRecovery *rec, double x, double y)
{
  return evaluateFieldFromRecovery(14, rec, x, y);
}

double dtt(const double x, const double y, double C, double S, feRecovery *rec)
{
  const double c11 = fxx(rec, x, y);
  const double c12 = (fxy(rec, x, y) + fyx(rec, x, y)) / 2.;
  const double c22 = fyy(rec, x, y);
  return C * C * c11 + S * S * c22 + 2. * C * S * c12;
}

double dttt(const double x, const double y, double C, double S, feRecovery *rec, int direction)
{
  const double c111 = fxxx(rec, x, y);
  const double c222 = fyyy(rec, x, y);
  const double c112 = (fxxy(rec, x, y) + fxyx(rec, x, y) + fyxx(rec, x, y)) / 3.;
  const double c122 = (fxyy(rec, x, y) + fyyx(rec, x, y) + fyxy(rec, x, y)) / 3.;

  const double c11 = fxx(rec, x, y);
  const double c12 = (fxy(rec, x, y) + fyx(rec, x, y)) / 2.;
  const double c22 = fyy(rec, x, y);

  const double c1 = fx(rec, x, y);
  const double c2 = fy(rec, x, y);

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
double dtttAnalytical(const double x, const double y, double C, double S, feRecovery *rec,
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

void computeDirectionFieldFromGradient(double x, double y, double &C, double &S, double tol,
                                       feRecovery *rec, FILE *F)
{
  double a, b;
  a = fx(rec, x, y);
  b = fy(rec, x, y);

  double normGrad = sqrt(a * a + b * b);

  if(normGrad > tol) {
    double theta1 = atan2(b, a);
    C = cos(theta1);
    S = sin(theta1);
    if(F != nullptr) {
      fprintf(F, "VP(%g,%g,%g){%g,%g,%g};\n", x, y, 0., C, S, 0.);
    }
  } else {
    // Gradient norm is too small : directions will be smoothed
    C = 0.;
    S = 0.;
  }
}

void computeDirectionFieldFromHessian(double x, double y, double &C, double &S, double tol,
                                      feRecovery *rec, FILE *F)
{
  double a, b, c;
  a = fxx(rec, x, y);
  b = (fxy(rec, x, y) + fyx(rec, x, y)) / 2.;
  c = fyy(rec, x, y);

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
      fprintf(F, "VP(%g,%g,%g){%g,%g,%g};\n", x, y, 0., C, S, 0.);
    }
  } else {
    // Hessian's eigenvalues are too small : directions will be smoothed
    C = 0.;
    S = 0.;
  }
}

#if defined(HAVE_GMSH)
void smoothDirections(std::map<size_t, double> &C, std::map<size_t, double> &S, FILE *F, int nIter,
                      double tol)
{
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

  // FILE *f = fopen("dirs.pos", "w");
  // fprintf(f, "View\"Dirs\"{\n");

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

    if(F != nullptr) fprintf(F, "VP(%g,%g,0){%g,%g,0};", coord[0], coord[1], c, s);
  }

  // fprintf(f, "};\n");
  // fclose(f);
}

void smoothSizes(std::map<size_t, double> &L1, std::map<size_t, double> &L2, int nIter,
                      double tol)
{
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
}
#endif

/* Compute metric using a brute-force method to solve optimization problem
   Original code from R. Kuate */
void metricHechtKuate(int nbpoints, double *x, double *y, double &A, double &B, double &C,
                      double epsilon, double *xNew, double *yNew)
{
  C = 0.0;

  int bool_assert = 1;

  double epsilon0 = 1e-5, precision = 1e-18, delta = 1e-10;
  double inf = 1e100;

  double Rmin = 1e20, Rmax = 0;

  int indiceX0 = 0;

  // double* xNew = new double[nbpoints];
  // double* yNew = new double[nbpoints];
  // double *PPoint = new double[nbpoints];

  for(int i = 0; i < nbpoints; i++) {
    Rmax = fmax(Rmax, sqrt(x[i] * x[i] + y[i] * y[i]));

    //---déplacement des points situées sur les axes--------------
    if(abs(x[i]) <= precision) {
      if(y[i] < 0) {
        x[i] = -delta;
        y[i] = -sqrt(pow(y[i], 2) - pow(x[i], 2));
      }
      if(y[i] > 0) {
        x[i] = delta;
        y[i] = sqrt(pow(y[i], 2) - pow(x[i], 2));
      }
    }

    if(abs(y[i]) <= precision) {
      if(x[i] < 0) {
        y[i] = -delta;
        x[i] = -sqrt(pow(x[i], 2) - pow(y[i], 2));
      }
      if(x[i] > 0) {
        y[i] = delta;
        x[i] = sqrt(pow(x[i], 2) - pow(y[i], 2));
      }
    }
    //-----------------------------------------------------------
    if(bool_assert) {
      // std::cout<<"foo : "<<abs(x[i]*y[i])<<std::endl;
      // std::cout<<"foo : "<<pow(precision,2)<<std::endl;
      assert(abs(x[i] * y[i]) >= pow(precision, 2));
    }

    if(Rmin > sqrt(x[i] * x[i] + y[i] * y[i])) {
      indiceX0 = i;
      Rmin = sqrt(x[i] * x[i] + y[i] * y[i]);
    }
  }

  //-------permutation des indices de la liste des points :
  // ranger la liste en commençant par le point X0-------
  for(int k = 0; k < nbpoints - indiceX0; k++) {
    xNew[k] = x[k + indiceX0];
    yNew[k] = y[k + indiceX0];
  }
  for(int k = nbpoints - indiceX0; k < nbpoints; k++) {
    xNew[k] = x[k - nbpoints + indiceX0];
    yNew[k] = y[k - nbpoints + indiceX0];
  }
  for(int i = 0; i < nbpoints; i++) {
    x[i] = xNew[i];
    y[i] = yNew[i];
  }

  //----------------------------------------------------------------

  int test = -1;

  double X0, Y0;
  double bmin = 0., bmax = inf, b1, b2, aik = 0., bik = 0., cik = 0.;
  double Xk = 0., Yk = 0., Ck = 0., Bk = 0., Xi = 0., Yi = 0., ri, detXY = 0., Ri, R0, r0;

  X0 = x[0];
  Y0 = y[0];
  r0 = sqrt(x[0] * x[0] + y[0] * y[0]);
  if(bool_assert) {
    // std::cout<<"Assert ON"<<std::endl;
    assert(r0 == Rmin);
  }

  // std::cout<<" Rmin = "<<Rmin<<" Rmax =  "<<Rmax<<std::endl;

  double EPS = 0.0; // pour recuperer la valeur de epsilon0 optimale

  int condition = -1;

  if(r0 <= epsilon0) epsilon0 = r0 * epsilon0;

  A = 1. / ((r0 - epsilon0) * (r0 - epsilon0));
  B = A;

  double epsilon0min = epsilon0;

  // if(abs(Rmin-Rmax) > 1e-5){
  if(abs(Rmin - Rmax) > precision) {
    // for(int ee=0; ee<neps-1; ee++){ //boucle sur epsilon0---------------
    // epsilon0= Tabepsilon[ee];
    if(r0 <= epsilon0) epsilon0 = r0 * epsilon0;
    if(bool_assert) {
      assert(r0 > epsilon0);
    }
    R0 = r0 / (r0 - epsilon0);

    for(int i = 1; i < nbpoints; i++) { // boucle sur chaque noeud
      Xi = x[i];
      Yi = y[i];
      ri = sqrt(x[i] * x[i] + y[i] * y[i]);

      if(ri <= epsilon) epsilon = ri * epsilon;

      if(bool_assert) {
        assert(ri > epsilon);
      }

      Ri = ri / (ri - epsilon);

      detXY = Xi * Y0 - Yi * X0;

      //------deplacement des points alignés avec l'origine et X0-----------
      if(abs(detXY) <= precision) {
        printf("Point %d - \t x = %10.15f - \t y =%10.15f - \t x0 = %10.15f - \t y0 =%10.15f - \t "
               "detXY = %10.15f \n",
               i, x[i], y[i], X0, Y0, detXY);
        Xi += delta;

        if(Yi < 0)
          Yi = -sqrt(pow(ri, 2) - pow(Xi, 2));
        else
          Yi = sqrt(pow(ri, 2) - pow(Xi, 2));
        x[i] = Xi;
        y[i] = Yi;

        std::cout << ri << std::endl;
        ri = sqrt(x[i] * x[i] + y[i] * y[i]);

        if(ri <= epsilon) epsilon = ri * epsilon;

        std::cout << "xi =" << x[i] << " yi =" << y[i] << " ri = " << ri << " epsilon = " << epsilon
                  << std::endl;
        if(bool_assert) {
          assert(ri > epsilon);
        }
        Ri = ri / (ri - epsilon);
      }

      detXY = Xi * Y0 - Yi * X0;

      if(bool_assert) {
        assert(abs(detXY) >= precision);
      }

      //-----racines du polynome en b à minimiser----------------------------
      double bb1 =
        (1. / pow(detXY, 2)) *
        (pow(X0 * Ri, 2) + pow(Xi * R0, 2) -
         2. * abs(Xi * X0) * sqrt(pow(R0 * Ri, 2) - pow(detXY / (Rmax * (r0 - epsilon0)), 2)));
      double bb2 =
        (1. / pow(detXY, 2)) *
        (pow(X0 * Ri, 2) + pow(Xi * R0, 2) +
         2. * abs(Xi * X0) * sqrt(pow(R0 * Ri, 2) - pow(detXY / (Rmax * (r0 - epsilon0)), 2)));
      //--fin----racines du polynome en b à minimiser--------------------

      bmax = fmin(bb2, pow(Rmax / pow((r0), 2), 2));
      bmin = fmax(1. / (Rmax * Rmax), bb1); // minoration de b
      double Cte = fmax(1e-9, (bmax - bmin) * 1e-9);
      bmin = bmin * (1. + Cte);
      bmax = bmax * (1. - Cte);

      // bornes de b-----------------------------------------------------------

      // cas:  majoration de c --------------------------------------------
      double Li = X0 * Xi * (pow(Rmax / pow(r0 - epsilon0min, 2), 2) - 1. / pow(Rmax, 2)) +
                  (pow(Ri * X0, 2) - pow(R0 * Xi, 2)) / detXY;
      double LiXY = Xi * Y0 + Yi * X0;

      if(abs(LiXY) >= precision) {
        condition = 1;

        if(Xi * X0 > 0) {
          if(LiXY > 0)
            bmin = fmax(bmin, -Li / LiXY);
          else
            bmax = fmin(bmax, -Li / LiXY);
        } else {
          if(LiXY < 0)
            bmin = fmax(bmin, -Li / LiXY);
          else
            bmax = fmin(bmax, -Li / LiXY);
        }
      } else {
        if(Li < 0)
          condition = 0;
        else
          condition = 1;
      }

      // cas  minoration de c --------------------------------------------
      Li = X0 * Xi * (-pow(Rmax / pow(r0 - epsilon0min, 2), 2) + 1. / pow(Rmax, 2)) +
           (pow(Ri * X0, 2) - pow(R0 * Xi, 2)) / detXY;
      LiXY = Xi * Y0 + Yi * X0;

      if(abs(LiXY) >= precision) {
        condition = 1;
        if(Xi * X0 > 0) {
          if(LiXY < 0)
            bmin = fmax(bmin, -Li / LiXY);
          else
            bmax = fmin(bmax, -Li / LiXY);
        } else {
          if(LiXY > 0)
            bmin = fmax(bmin, -Li / LiXY);
          else
            bmax = fmin(bmax, -Li / LiXY);
        }
      } else {
        if(Li > 0)
          condition = 0;
        else
          condition = 1;
      }

      if(condition) {
        //--cas : minoration de a-----------------------------------------------

        double Gi =
          ((Xi * Yi * R0 * R0 - X0 * Y0 * Ri * Ri) / detXY + Xi * X0 / (Rmax * Rmax)) / (Yi * Y0);

        if(Xi * X0 > 0) {
          if(Yi * Y0 > 0)
            bmin = fmax(bmin, Gi);
          else
            bmax = fmin(bmax, Gi);
        } else {
          if(Yi * Y0 < 0)
            bmin = fmax(bmin, Gi);
          else
            bmax = fmin(bmax, Gi);
        }

        // cas :majoration de a------------------------------------------------
        double Hi = (Xi * X0 * Rmax * Rmax / pow((r0 - epsilon0min), 4) +
                     (Xi * Yi * R0 * R0 - X0 * Y0 * Ri * Ri) / detXY) /
                    (Yi * Y0);
        if(Xi * X0 > 0) {
          if(Yi * Y0 > 0)
            bmax = fmin(bmax, Hi);
          else
            bmin = fmax(bmin, Hi);
        } else {
          if(Yi * Y0 < 0)
            bmax = fmin(bmax, Hi);
          else
            bmin = fmax(bmin, Hi);
        }
        //------fin bornes de b------------------------------------------------
        b2 = bmax;
        b1 = bmin;

        for(int k = 1; k < nbpoints; k++) { // on balaye les contraintes
          Xk = x[k];
          Yk = y[k];
          Bk = (Yk * Yk * Xi * X0 + Xk * (Xk * Yi * Y0 - Yk * (Yi * X0 + Xi * Y0))) / (Xi * X0);
          Ck = (X0 * Xi * detXY -
                Xk * (Xi * R0 * R0 * (Yk * Xi - Yi * Xk) + X0 * Ri * Ri * (-Yk * X0 + Y0 * Xk))) /
               (Xi * X0 * detXY);

          if(bool_assert) {
            assert(abs(Xi * X0 * Y0 * Yi * Xk * Yk) >= pow(precision, 5));
          }
          if(abs(Bk) > precision) { // non nul

            if(Bk <= 0)
              bmax = fmin(bmax, Ck / Bk);
            else
              bmin = fmax(bmin, Ck / Bk);

            if((bmax < b1) || (bmin > b2) || (bmin > bmax)) {
              test = 0;
              break;
            }

            else
              test = 1;
          } else {
            if(Ck > precision) {
              test = 0;
              break;
            } else // Ck<=0
              test = -1; // 1 peut etre
          }
        }

        if(test) {
          double a0 = -pow((detXY / (Xi * X0)), 2);
          double a1 = 2. * (pow(Ri / Xi, 2) + pow(R0 / X0, 2));
          if(((a0 * bmax + a1) * bmax) < ((a0 * bmin + a1) * bmin))
            bik = bmax;
          else
            bik = bmin;

          aik = (Ri * Ri * Y0 * X0 - R0 * R0 * Yi * Xi + bik * Yi * Y0 * detXY) / (detXY * Xi * X0);
          cik = (-Ri * Ri * X0 * X0 + R0 * R0 * Xi * Xi - bik * (Yi * X0 + Y0 * Xi) * detXY) /
                (detXY * Xi * X0);

          if(bool_assert) {
            assert((4. * aik * bik - cik * cik) >= 0.); // aire positive
            assert(abs((4. * aik * bik - cik * cik) - pow(2. / (Rmax * (r0 - epsilon0)), 2)) >
                   0); // aire positive
          }
          if((4. * aik * bik - cik * cik) <= (4. * A * B - C * C)) {
            A = aik;
            B = bik;
            C = cik;
            EPS = epsilon0;
          }
        } // if(test)
      } // if(condition)
    } // for(int i=1;i<nbpoints;i++)
    // } // for(int ee=0; ee<neps-1; ee++)
  } // if(abs(Rmin-Rmax)>1e-5)
  else {
    A = 1. / (Rmin * Rmin);
    B = A;
    C = 0.;
  }
}

void metriqueSimplexe2D(int nPhi, std::vector<double> phi, std::vector<double> erreur, double &A,
                        double &B, double &C, int max_iter, std::vector<int> &expX,
                        std::vector<int> &expY, int dimRecovery, int degreeSolution)
{
  double xi, yi, xj, yj, Q11, Q12, Q21, Q22;
  double err, normeXj;
  double L1, L2, L3, expL11, expL12, expL21, expL22;
  double Aprev, Bprev, Cprev;
  // Variables LAPACK
  // char jobz='v', uplo='l';
  // const int n=2, lwork=18;
  // int info;
  // double w[2], work[18];

  SVector3 v0, v1, v2;

  // std::vector<int> &expX = _recovery->getXExponentsRecovery();
  // std::vector<int> &expY = _recovery->getYExponentsRecovery();

  int nSimplex = 6, mSimplex = nPhi;

  std::vector<double> AConstraint(nPhi * nSimplex), bConstraint(nPhi), cObjective(nSimplex);
  cObjective[0] = -1.0;
  cObjective[1] = 1.0;
  cObjective[2] = 0.0;
  cObjective[3] = 0.0;
  cObjective[4] = -1.0;
  cObjective[5] = 1.0;

  A = 1.0; // Candidat initial : Q = I
  C = 0.0;
  B = 1.0;

  std::pair<std::vector<double>, double> retSimplex;

  int i = 0;
  double residu = 1.0;
  while(i < max_iter && residu > 1e-8) {
    Aprev = A;
    Bprev = B;
    Cprev = C;
    //==================================================================================
    // Calcul de Q^(-1/2)
    SMetric3 Q(1.0);
    Q.set_m11(A);
    Q.set_m21(C);
    Q.set_m22(B);
    fullMatrix<double> V(3, 3);
    fullVector<double> S(3);
    Q.eig(V, S, false);

    // S(0) = pow(S(0), -0.5);
    // S(1) = pow(S(1), -0.5);
    // S(2) = pow(S(2), -0.5);

    v0 = SVector3(V(0, 0), V(0, 1), V(0, 2)); // Attention c'est peut-être la transposée
    v1 = SVector3(V(1, 0), V(1, 1), V(1, 2));
    v2 = SVector3(V(2, 0), V(2, 1), V(2, 2));

    Q = SMetric3(pow(S(0), -0.5), pow(S(1), -0.5), pow(S(2), -0.5), v0, v1, v2);

    // double aTest[4] = {A, C, C, B};
    // dsyev(&jobz, &uplo, &n, aTest, &n, w, work, &lwork, &info);

    // l1 = w[0];  v11 = aTest[0];  v21 = aTest[2];
    // l2 = w[1];  v12 = aTest[1];  v22 = aTest[3];

    // Q11 = v11*v11*pow(l1, -0.5) + v21*v21*pow(l2, -0.5);
    // Q12 = v12*v11*pow(l1, -0.5) + v21*v22*pow(l2, -0.5);
    // Q21 = v11*v12*pow(l1, -0.5) + v21*v22*pow(l2, -0.5);
    // Q22 = v12*v12*pow(l1, -0.5) + v22*v22*pow(l2, -0.5);
    //==================================================================================
    for(int iPhi = 0; iPhi < nPhi; ++iPhi) {
      xi = cos(phi[iPhi]);
      yi = sin(phi[iPhi]);

      xj = Q(0, 0) * xi + Q(0, 1) * yi;
      yj = Q(1, 0) * xi + Q(1, 1) * yi;

      // err = erreur[0] * xj*xj + erreur[1] * xj*yj + erreur[2] * yj*yj;
      err = 0.0;
      int indice = 0;
      for(int iDeg = 0; iDeg < dimRecovery; ++iDeg) {
        if((expX[iDeg] + expY[iDeg]) == degreeSolution + 1) {
          err += erreur[indice] * pow(xj, expX[iDeg]) * pow(yj, expY[iDeg]);
          ++indice;
        }
      }

      xj = xi / fmax(pow(fabs(err), 1.0 / (degreeSolution + 1)), 1e-20);
      yj = yi / fmax(pow(fabs(err), 1.0 / (degreeSolution + 1)), 1e-20);

      // ==============================================================================
      // Contraintes lineaires
      AConstraint[iPhi * nSimplex] = xj * xj;
      AConstraint[iPhi * nSimplex + 1] = -xj * xj;
      AConstraint[iPhi * nSimplex + 2] = 2.0 * xj * yj;
      AConstraint[iPhi * nSimplex + 3] = -2.0 * xj * yj;
      AConstraint[iPhi * nSimplex + 4] = yj * yj;
      AConstraint[iPhi * nSimplex + 5] = -yj * yj;

      normeXj = xj * xj + yj * yj;

      bConstraint[iPhi] = -(-normeXj * log(normeXj));
    }

    retSimplex = simplex(nSimplex, mSimplex, AConstraint, bConstraint, cObjective, 0.0);

    if(isinf(retSimplex.second)) {
      if(retSimplex.first[0] == -1)
        printf("Objective function unbounded!\n");
      else if(retSimplex.first[0] == -2)
        printf("Linear program infeasible!\n");
    } else {
      // printf("Solution: (");
      // for (int i=0;i<nSimplex+mSimplex;i++) printf("%lf%s", retSimplex.first[i], (i < nSimplex +
      // mSimplex - 1) ? ", " : ")\n"); for (int i=0;i<nSimplex;i++) printf("%lf%s",
      // retSimplex.first[i], (i < nSimplex - 1) ? ", " : ")\n"); printf("Optimal objective value:
      // %lf\n", retSimplex.second);
      L1 = retSimplex.first[0] - retSimplex.first[1];
      L2 = retSimplex.first[2] - retSimplex.first[3];
      L3 = retSimplex.first[4] - retSimplex.first[5];
    }

    //==================================================================================
    // Calcul de Q^(1/2)
    // Q11 = v11*v11*pow(l1, 0.5) + v21*v21*pow(l2, 0.5);
    // Q12 = v12*v11*pow(l1, 0.5) + v21*v22*pow(l2, 0.5);
    // Q21 = v11*v12*pow(l1, 0.5) + v21*v22*pow(l2, 0.5);
    // Q22 = v12*v12*pow(l1, 0.5) + v22*v22*pow(l2, 0.5);
    Q = SMetric3(pow(S(0), 0.5), pow(S(1), 0.5), pow(S(2), 0.5), v0, v1, v2);
    Q11 = Q(0, 0);
    Q12 = Q(0, 1);
    Q21 = Q(1, 0);
    Q22 = Q(1, 1);
    //==================================================================================
    // Calcul de expm(L)

    SMetric3 L(1.0);
    L.set_m11(L1);
    L.set_m21(L2);
    L.set_m22(L3);
    // fullMatrix<double> V(3, 3);
    // fullVector<double> S(3);
    L.eig(V, S, false);

    // S(0) = pow(S(0), -0.5);
    // S(1) = pow(S(1), -0.5);
    // S(2) = pow(S(2), -0.5);

    v0 = SVector3(V(0, 0), V(0, 1), V(0, 2)); // Attention c'est peut-être la transposée
    v1 = SVector3(V(1, 0), V(1, 1), V(1, 2));
    v2 = SVector3(V(2, 0), V(2, 1), V(2, 2));

    L = SMetric3(exp(S(0)), exp(S(1)), exp(S(2)), v0, v1, v2);

    expL11 = L(0, 0);
    expL12 = L(0, 1);
    expL21 = L(1, 0);
    expL22 = L(1, 1);

    // double aTest2[4] = {L1, L2, L2, L3};
    // dsyev(&jobz, &uplo, &n, aTest2, &n, w, work, &lwork, &info);

    // l1 = w[0];  v11 = aTest2[0];  v21 = aTest2[2];
    // l2 = w[1];  v12 = aTest2[1];  v22 = aTest2[3];

    // expL11 = v11*v11*exp(l1) + v21*v21*exp(l2);
    // expL12 = v12*v11*exp(l1) + v21*v22*exp(l2);
    // expL21 = v11*v12*exp(l1) + v21*v22*exp(l2);
    // expL22 = v12*v12*exp(l1) + v22*v22*exp(l2);
    //==================================================================================
    // Calcul de Q = Q^(1/2) * expm(L) * Q^(1/2)
    A = Q11 * (Q11 * expL11 + Q12 * expL21) + Q21 * (Q11 * expL12 + Q12 * expL22);
    C = Q12 * (Q11 * expL11 + Q12 * expL21) + Q22 * (Q11 * expL12 + Q12 * expL22);
    B = Q12 * (Q21 * expL11 + Q22 * expL21) + Q22 * (Q21 * expL12 + Q22 * expL22);

    residu =
      sqrt((A - Aprev) * (A - Aprev) + 2 * (C - Cprev) * (C - Cprev) + (B - Bprev) * (B - Bprev));
    // cout<<"======== ITERATION "<<i+1<<" - RESIDU = "<<residu<<" ==========="<<endl;
    ++i;
  }
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
void computeWorstMetric(int nTheta, int nIncr, double e, double x, double y, double cG, double sG, 
  double &h1Opt, double &h2Opt, feRecovery *rec, SMetric3 &M, double lMin, double lMax){

  double theta, k, C, H, s;
  double v[2], w[2], g1[2], g2[2], xdot[2], xdotdot[2];

  std::vector<double> xE(4*nTheta, 0.);
  std::vector<double> yE(4*nTheta, 0.);

  g1[0] = cG;
  g1[1] = sG;
  g2[0] = -sG;
  g2[1] = cG;

  double c1 = fx(rec, x, y);
  double c2 = fy(rec, x, y);

  double c11 = fxx(rec, x, y);
  double c12 = (fxy(rec, x, y) + fyx(rec, x, y)) / 2.;
  double c22 = fyy(rec, x, y);

  double c111 = fxxx(rec, x, y);
  double c222 = fyyy(rec, x, y);
  double c112 = (fxxy(rec, x, y) + fxyx(rec, x, y) + fyxx(rec, x, y)) / 3.;
  double c122 = (fxyy(rec, x, y) + fyyx(rec, x, y) + fyxy(rec, x, y)) / 3.;

  double grad[2] = {c1, c2};
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
      xE[iDir*nTheta+i] = x + v[0]*s + k*w[0]*s*s/2.;
      yE[iDir*nTheta+i] = y + v[1]*s + k*w[1]*s*s/2.;

      // // // ////////////////////////////////////////////////

      // // Bissection :

      // double sA = 0.1;
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
          d = sqrt( (xE[k] - x)*(xE[k] - x) + (yE[k] - y)*(yE[k] - y) );
          theta = atan2(yE[k] - y, xE[k] - x);
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

  double tol = 1e-1;
  if(fabs(x - 0.158414) < tol && fabs(y - 0.453362) < tol){
    FILE *myf = fopen("errorCurveBissection.pos", "w");
    fprintf(myf, "View \" \"{\n");

    // Print the error level
    for(int i = 0; i < xE.size(); ++i) {
      fprintf(myf, "SP(%.16g,%.16g,%.16g){%f};\n", xE[i], yE[i], 0., 1.);
    }
    // Print the kuate ellipse
    int nt = 30;
    std::vector<double> xP(nt, 0.);
    std::vector<double> yP(nt, 0.);
    getEllipsePoints(M(0, 0), 2.0 * M(0, 1), M(1, 1), x, y, xP, yP);
    for(int j = 0; j < nt; ++j) {
      if(j != nt - 1) {
        fprintf(myf, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
                0., xP[j + 1], yP[j + 1], 0., 1, 1);
      } else {
        fprintf(myf, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
                0., xP[0], yP[0], 0., 1, 1);
      }
    }

    fprintf(myf, "};");
    fclose(myf);
  }
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

static inline SMetric3 intersectionReductionSimultaneeExplicite(const SMetric3 &m1, const SMetric3 &m2){
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
    if( fabs(a2/a1 - c2/c1) < 1e-6 && fabs(a2/a1 - b2/b1) < 1e-6 ){
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
  return tmp;
}

void gradationMetriques(double gradation, std::vector<double> &coord1, std::map<int, SMetric3> &metricsOnGmshModel)
{
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

  int nodedim, nodetag;

  FILE* fff = fopen("CHECK.pos","w");
  fprintf(fff, "View \" CHECK \"{\n");

  size_t maxtag= -1;
  std::pair<std::set<std::pair<size_t,size_t>, gmshEdgeLessThan>::iterator, bool> ret1;
  std::pair<std::set<std::pair<size_t,size_t>, gmshEdgeLessThan>::iterator, bool> ret2;
  std::pair<std::set<std::pair<size_t,size_t>, gmshEdgeLessThan>::iterator, bool> ret3;
  for(size_t i = 0; i < elementTags[0].size(); i++) {
  // for(size_t i = 0; i < 5; i++) {
    size_t n0 = elemNodeTags[0][3 * i + 0];
    size_t n1 = elemNodeTags[0][3 * i + 1];
    size_t n2 = elemNodeTags[0][3 * i + 2];
    ret1 = edges.insert(std::make_pair(n0, n1));
    ret2 = edges.insert(std::make_pair(n1, n2));
    ret3 = edges.insert(std::make_pair(n2, n0));
    if(!ret1.second){
      feInfo("Not added : elem %d edge 1 (%d,%d)", elementTags[0][i], n0, n1);
    }
    if(!ret2.second){
      feInfo("Not added : elem %d edge 2 (%d,%d)", elementTags[0][i], n1, n2);
    }
    if(!ret3.second){
      feInfo("Not added : elem %d edge 3 (%d,%d)", elementTags[0][i], n2, n0);
    }
    // feInfo("Added elem %d edge 1 (%d,%d)", elementTags[0][i], n0, n1);
    // feInfo("Added elem %d edge 2 (%d,%d)", elementTags[0][i], n1, n2);
    // feInfo("Added elem %d edge 3 (%d,%d)", elementTags[0][i], n2, n0);
    // double x1 = coord[3*n0+0];
    // double y1 = coord[3*n0+1];
    // double x2 = coord[3*n1+0];
    // double y2 = coord[3*n1+1];
    // double x3 = coord[3*n2+0];
    // double y3 = coord[3*n2+1];
    // fprintf(fff, "SL(%.16g,%.16g,0.,%.16g,%.16g,0.){%u, %u};\n", x1, y1, x2, y2, 1, 1);
    // fprintf(fff, "SL(%.16g,%.16g,0.,%.16g,%.16g,0.){%u, %u};\n", x2, y2, x3, y3, 1, 1);
    // fprintf(fff, "SL(%.16g,%.16g,0.,%.16g,%.16g,0.){%u, %u};\n", x3, y3, x1, y1, 1, 1);
    // x1 = coord[3*(n0-1)+0];
    // y1 = coord[3*(n0-1)+1];
    // x2 = coord[3*(n1-1)+0];
    // y2 = coord[3*(n1-1)+1];
    // x3 = coord[3*(n2-1)+0];
    // y3 = coord[3*(n2-1)+1];
    // fprintf(fff, "SL(%.16g,%.16g,0.,%.16g,%.16g,0.){%u, %u};\n", x1, y1, x2, y2, 1, 1);
    // fprintf(fff, "SL(%.16g,%.16g,0.,%.16g,%.16g,0.){%u, %u};\n", x2, y2, x3, y3, 1, 1);
    // fprintf(fff, "SL(%.16g,%.16g,0.,%.16g,%.16g,0.){%u, %u};\n", x3, y3, x1, y1, 1, 1);
    // x1 = coord[3*0];
    // y1 = coord[3*0+1];
    // fprintf(fff, "SP(%.16g,%.16g,0.){%u};\n", x1, y1, 1);
    // feInfo("nodeTags[0] = %d", nodeTags[0]);
    // maxtag = fmax(maxtag,n0);
    // maxtag = fmax(maxtag,n1);
    // maxtag = fmax(maxtag,n2);

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
  int iter = 0, nmax = 20;
  while(correction && iter < nmax){

    correction = false;
    iter++;

    FILE* fff = fopen("CHECK.pos","w");
    fprintf(fff, "View \" CHECK \"{\n");

    size_t n0 = elemNodeTags[0][400];
    gmsh::model::mesh::getNode(n0, nodecoord3, nodeparametricCoord, nodedim, nodetag);
    double X0 = nodecoord3[0];
    double Y0 = nodecoord3[1];
    SMetric3 mref = metricsOnGmshModel[n0];
    fprintf(fff, "SP(%.16g,%.16g,0.){%u};\n", X0, Y0, 1);

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
      mpq = intersectionReductionSimultaneeExplicite(mq, mpq);

      // mpq = intersection_reductionSimultanee(mq, mpq);
      // metricsOnGmshModel[neigh] = mpq;

      if(matNorm2(metricsOnGmshModel[neigh], mpq)/matNorm2(mpq) > tolerance_gradation){
        metricsOnGmshModel[neigh] = mpq;
        correction = true;
      }
      // Span mq to n
      mqp = gradationAlauzet(gradation, mq, xq, yq, xp, yp);
      mqp = intersectionReductionSimultaneeExplicite(mp, mqp);

      // mqp = intersection_reductionSimultanee(mp, mqp);
      // metricsOnGmshModel[n] = mqp;

      if(matNorm2(metricsOnGmshModel[n], mqp)/matNorm2(mqp) > tolerance_gradation){
        metricsOnGmshModel[n] = mqp;
        correction = true;
      }

      mp = gradationAlauzet(gradation, mref, X0, Y0, xp, yp);

      double factor = 1000.;
      getEllipsePoints(factor * mp(0, 0), factor * 2.0 * mp(0, 1), factor * mp(1, 1), xp, yp, xP, yP);
      for(int j = 0; j < nt; ++j) {
        if(j != nt - 1) {
          fprintf(fff, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
                  0., xP[j + 1], yP[j + 1], 0., 1, 1);
        } else {
          fprintf(fff, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
                  0., xP[0], yP[0], 0., 1, 1);
        }
      }
    }

    feInfo("Passe %d", iter);
    // correction = false;

    fprintf(fff, "};");
    fclose(fff);

  }
  feInfo("%d noeuds dans ", maxtag);
  feInfo("%d coord size", coord.size());
  feInfo("%d edges size", edges.size());
}