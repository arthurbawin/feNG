#include "feMetricTools.h"
#include "feSimplex.h"

#if defined(HAVE_GMSH)
#include "gmsh.h"
#endif
#include "fullMatrix.h"
#include "STensor3.h"

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

double sech(double x){
  return 1./cosh(x);
}

// For the solution u = (1. + tanh(a*(r-r0)))/2.
double dtttAnalytical(const double x, const double y, double C, double S, feRecovery *rec, int direction)
{
  double a = 10.0;
  double x0 = 0.2;
  double y0 = 0.4;
  double r0 = 0.2;

  double r = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));
  double drdx = (x-x0)/r;
  double drdy = (y-y0)/r;

  double dudr = a/2. * sech(a*(r-r0)) * sech(a*(r-r0));
  double d2udr2 = -a*a * tanh(a*(r-r0)) * sech(a*(r-r0)) * sech(a*(r-r0));
  double d3udr3 = -a*a*a * sech(a*(r-r0)) * sech(a*(r-r0)) * (sech(a*(r-r0)) * sech(a*(r-r0)) - 2. * tanh(a*(r-r0)) * tanh(a*(r-r0)));

  const double c1 = dudr*drdx;
  const double c2 = dudr*drdy;

  const double c11 = d2udr2*drdx*drdx;
  const double c12 = d2udr2*drdx*drdy;
  const double c22 = d2udr2*drdy*drdy;

  const double c111 = d3udr3*drdx*drdx*drdx;
  const double c222 = d3udr3*drdy*drdy*drdy;
  const double c112 = d3udr3*drdx*drdx*drdy;
  const double c122 = d3udr3*drdx*drdy*drdy;

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
  //   feInfo("C1 = %8.8e - H1 = %8.8e", g11 * g11 * g11 * c111 + g12 * g12 * g12 * c222 + 3. * g11 * g11 * g12 * c112 +
  //          3. * g11 * g12 * g12 * c122, g11 * g21 * c11 + g11 * g22 * c12 + g12 * g21 * c12 + g12 * g22 * c22);
  //   feInfo("C2 = %8.8e - H2 = %8.8e", g21 * g21 * g21 * c111 + g22 * g22 * g22 * c222 + 3. * g21 * g21 * g22 * c112 +
  //          3. * g21 * g22 * g22 * c122, g11 * g21 * c11 + g21 * g12 * c12 + g22 * g11 * c12 + g12 * g22 * c22);
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