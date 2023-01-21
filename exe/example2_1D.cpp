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
  return 1. / sqrt(4.*M_PI*k*t) * exp(-(x-x0)*(x-x0)/(4.*k*t));
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
  feFunction *funSol       = new feFunction(fSol,  {x0, k});
  feFunction *funZero      = new feFunction(fZero, {});
  feFunction *kDiffusivity = new feFunction(fConstant, {k});

  int dim, deg = 2, degreeQuadrature = 10;
  double xa = -2.;
  double xb = 1.;
  int nElm = 10;
  feMesh1DP1 mesh(xa, xb, nElm, "BXA", "BXB", "Domaine");
  feMesh2DP1 mesh2;

  feSpace *uG, *uD, *uDomaine;
  feCheck(createFiniteElementSpace(uG, &mesh, dim = 0, POINT, LAGRANGE, deg, "U", "BXA", degreeQuadrature, funSol));
  feCheck(createFiniteElementSpace(uD, &mesh, dim = 0, POINT, LAGRANGE, deg, "U", "BXB", degreeQuadrature, funSol));
  feCheck(createFiniteElementSpace(uDomaine, &mesh, dim = 1, LINE, LAGRANGE, deg, "U", "Domaine", degreeQuadrature, funSol));

  std::vector<feSpace*> spaces = {uG, uD, uDomaine};
  std::vector<feSpace*> essentialSpaces = {uG, uD};

  feMetaNumber numbering(&mesh, spaces, essentialSpaces);

  feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
  sol.initializeUnknowns(&mesh);

  feBilinearForm *mass, *diff;
  feCheck(createBilinearForm( mass, {uDomaine}, new feSysElm_1D_Masse(1.0, nullptr))   );
  feCheck(createBilinearForm( diff, {uDomaine}, new feSysElm_1D_Diffusion(kDiffusivity)) );

  feLinearSystem *system;
  feCheck(createLinearSystem(system, PETSC, {mass, diff}, numbering.getNbUnknowns(), argc, argv));

  feNorm normU(L2_ERROR, {uDomaine}, &sol, funSol);
  std::vector<feNorm *> norms = {&normU};

  feExportData exportData = {nullptr, 1, ""};

  TimeIntegrator *solver;
  feTolerances tol{1e-9, 1e-8, 50};
  double t0 = 0.1;
  double tEnd = 5;
  double nSteps = 500;
  feCheck(createTimeIntegrator(solver, BDF1, tol, system, &numbering, &sol, &mesh, norms, exportData, t0, tEnd, nSteps));

  int nInteriorPlotNodes = 10;
  feBasicViewer viewer("test", mesh.getNbInteriorElems(), nInteriorPlotNodes);

  double xLim[2] = {xa, xb};
  double yLim[2] = {0, 1.};
  viewer.setAxesLimits(xLim, yLim);

  int stepFrequency = 2;
  int nFrames = 0;

  do{
    viewer.reshapeWindowBox(1.2, mesh, sol);

    if(nFrames++ % stepFrequency == 0 && sol.getCurrentTime() < tEnd){
      feCheck(solver->makeStep());
    }

    viewer.draw1DCurve(mesh, numbering, sol, uDomaine, funSol);
    std::string message = "BDF1: iteration = " + std::to_string(solver->getCurrentStep()) + " - t = " + std::to_string(sol.getCurrentTime());
    viewer.drawMessage(message);
    viewer.windowUpdate();
  } while(!viewer.windowShouldClose());

  petscFinalize();
  return 0;
}