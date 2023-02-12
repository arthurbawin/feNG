/*
    Example 1: 1D steady-state convection-diffusion equation.

    We solve the following PDE in the unknown u(x) on [xa, xb]:

    div(k*grad(u)) - div(c*u) + f = 0, 

    with
      source term f(x) = -1
      velocity field c = 1
      variable (but homogenous in space) diffusivity k

    With constant diffusivity k and a divergence-free velocity field c,
    this simplifies to:

    k*lap(u) - c dot grad(u) + f = 0,

    which, in 1D, becomes:

    k*u''(x) - c*u'(x) + f = 0.

    The boundary conditions are u(xa) = u(xb) = 0. In this case, the exact
    solution is given by:

            1   /     1 - exp(cx/k) \
    u(x) = --- | x - -------------   |
            c   \     1 - exp(c/k)  /
*/

#include "feAPI.h"
#include "glfem.h"

// The function handling the exact solution u(t,x).
// The arguments for user-defined functions must always be:
//
//  t   : the time at which to evaluate fun (unused here)
//  pos : the position vector, with pos[0] = x, pos[1] = y, pos[2] = z
//  par : the parameters vector (may be empty). These must be consistent with the parameters
//        given we creating the function in the main function. Here, the exact solution
//        has two parameters, c and k. They are given in the main in the order {c, k}
//        hence we must respect this order here.
double fSol(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  double x = pos[0];
  double x0 = par[0];
  double k = par[1];
  // return 1./2. * (1 - erf(x / sqrt(4. * k * t)));
  // return 1. / sqrt(4.*M_PI*k*t) * exp(-(x-x0)*(x-x0)/(4.*k*t));

  double L = 1.;
  double res = 0.;
  double res1 = 0.;
  double res2 = 0.;
  for(int n = 1; n < 100; n+=2){
    double l = n*M_PI/L;
    // double An = 2. * (L*L* ( cos(M_PI*n)*(-M_PI*M_PI*(L-1.)*n*n + 2.*L) + M_PI*(2.*L-1.)*n*sin(M_PI*n) - 2.*L)) / (M_PI*M_PI*M_PI * n*n*n); // u0 = x(1-x)
    double An = -8. * L*L*L / (M_PI*M_PI*M_PI * n*n*n); // u0 = x(1-x)
    double Bn = 4./(M_PI*n); // u0 = 1

    // res += An * sin(n*M_PI*x/L) * exp(-k * (n*M_PI/L)*(n*M_PI/L) *t );
    res1 += An * sin(l*x) * exp(-k*l*l*t) + Bn * sin(l*x) * (1 - exp(-k*l*l*t))/(k*l*l);
  }
  return res1;
}

double fIni(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  double x = pos[0];
  return -x*(1.-x);
  return sin(M_PI*x);
  return (x < 0) ? 1. : 0.;
}

// A zero function to initialize the values of the degrees of freedom
double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  return 0.;
}

double fConstant(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  return par[0];
}

int main(int argc, char **argv)
{
  petscInitialize(argc, argv);

  double x0 = 0.;
  double k = 0.1;
  feFunction *funInitial   = new feFunction(fIni,  {});
  feFunction *funSol       = new feFunction(fSol,  {x0, k});
  feFunction *funZero      = new feFunction(fZero, {});
  feFunction *kDiffusivity = new feFunction(fConstant, {-k});
  feFunction *fSource = new feFunction(fConstant, {-1.0});

  int deg = 2, degreeQuadrature = 10;
  double xa = 0.;
  double xb = 1.;
  int nElm = 50;
  feMesh1DP1 mesh(xa, xb, nElm, "BXA", "BXB", "Domaine");
  feMesh2DP1 mesh2;

  feSpace *uG, *uD, *uDomaine;
  feCheck(createFiniteElementSpace(uG, &mesh, elementType::LAGRANGE, deg, "U", "BXA", degreeQuadrature, funZero));
  feCheck(createFiniteElementSpace(uD, &mesh, elementType::LAGRANGE, deg, "U", "BXB", degreeQuadrature, funZero));
  feCheck(createFiniteElementSpace(uDomaine, &mesh, elementType::LAGRANGE, deg, "U", "Domaine", degreeQuadrature, funInitial));

  std::vector<feSpace*> spaces = {uG, uD, uDomaine};
  std::vector<feSpace*> essentialSpaces = {uG, uD};

  feMetaNumber numbering(&mesh, spaces, essentialSpaces);

  feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
  sol.initializeUnknowns(&mesh);

  feBilinearForm *mass, *diff, *source;
  feCheck(createBilinearForm( mass, {uDomaine}, new feSysElm_TransientMass(1.0))   );
  feCheck(createBilinearForm( diff, {uDomaine}, new feSysElm_Diffusion<1>(kDiffusivity)) );
  feCheck(createBilinearForm( source, {uDomaine}, new feSysElm_Source(fSource)) );

  feLinearSystem *system;
  feCheck(createLinearSystem(system, PETSC, {mass, diff, source}, numbering.getNbUnknowns(), argc, argv));

  feNorm normU(L2_ERROR, {uDomaine}, &sol, funSol);
  std::vector<feNorm *> norms = {&normU};

  feExportData exportData = {nullptr, 1, ""};

  TimeIntegrator *solver;
  feTolerances tol{1e-9, 1e-8, 50};
  double t0 = 0;
  double tEnd = 1;
  double nSteps = 500;
  feCheck(createTimeIntegrator(solver, BDF1, tol, system, &numbering, &sol, &mesh, norms, exportData, t0, tEnd, nSteps));

  int nInteriorPlotNodes = 10;
  feBasicViewer viewer("test", mesh.getNumInteriorElements(), nInteriorPlotNodes);

  double xLim[2] = {xa, xb};
  double yLim[2] = {0, 1.};
  viewer.setAxesLimits(xLim, yLim);

  int stepFrequency = 2;
  int nFrames = 0;

  do{
    viewer.reshapeWindowBox(1.2, 1.2, mesh, sol);

    if(nFrames++ % stepFrequency == 0 && sol.getCurrentTime() < tEnd){
      feCheck(solver->makeStep());
      // double umax = 0.;
      // for(auto val : sol.getSolutionReference()){
      //   umax = fmax(umax, fabs(val));
      //   feInfo("uMax = %f", umax);
      // }
    }

    viewer.draw1DCurve(mesh, numbering, sol, uDomaine, funSol);
    std::string message = "BDF1: iteration = " + std::to_string(solver->getCurrentStep()) + " - t = " + std::to_string(sol.getCurrentTime());
    viewer.drawMessage(message);
    viewer.windowUpdate();
  } while(!viewer.windowShouldClose());

  petscFinalize();
  return 0;
}