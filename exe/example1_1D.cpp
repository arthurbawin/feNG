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
  double c = par[0];
  double k = par[1];
  return 1./c * (x - (1. - exp(c * x / k))/(1. - exp(c / k)));
}

// The source term f(x)
double fSource(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  return -1.;
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

  double c = 5.0;
  double k = 1.0;
  feFunction *funSol    = new feFunction(fSol,    {c, k});
  feFunction *funSource = new feFunction(fSource, {});
  feFunction *funZero   = new feFunction(fZero,   {});
  feFunction *cVelocity = new feFunction(fConstant, {c});

  int dim, deg = 1, degreeQuadrature = 10;
  double xa = 0.;
  double xb = 1.;
  int nElm = 500;
  feMesh1DP1 mesh(xa, xb, nElm, "BXA", "BXB", "Domaine");

  feSpace *uG, *uD, *uDomaine;
  feCheck(createFiniteElementSpace(uG, &mesh, dim = 0, POINT, LAGRANGE, deg, "U", "BXA", degreeQuadrature, funSol));
  feCheck(createFiniteElementSpace(uD, &mesh, dim = 0, POINT, LAGRANGE, deg, "U", "BXB", degreeQuadrature, funSol));
  feCheck(createFiniteElementSpace(uDomaine, &mesh, dim = 1, LINE, LAGRANGE, deg, "U", "Domaine", degreeQuadrature, funZero));

  std::vector<feSpace*> spaces = {uG, uD, uDomaine};
  std::vector<feSpace*> essentialSpaces = {uG, uD};

  feMetaNumber numbering(&mesh, spaces, essentialSpaces);

  feSolution sol(&mesh, spaces, essentialSpaces, &numbering);
  sol.initializeUnknowns(&mesh, &numbering);

  feBilinearForm *adv, *diff, *source;
  feCheck(createBilinearForm(    adv, {uDomaine}, &mesh, degreeQuadrature, new feSysElm_1D_Advection(1.0, cVelocity)) );
  feCheck(createBilinearForm(   diff, {uDomaine}, &mesh, degreeQuadrature, new feSysElm_1D_Diffusion(k, nullptr))     );
  feCheck(createBilinearForm( source, {uDomaine}, &mesh, degreeQuadrature, new feSysElm_1D_Source(1.0, funSource))    );

  feLinearSystem *system;
  feCheck(createLinearSystem(system, PETSC, spaces, {adv, diff, source}, &numbering, &mesh, argc, argv));

  feNorm normU(uDomaine, &mesh, degreeQuadrature, funSol);
  std::vector<feNorm *> norms = {&normU};

  feExporter *exporter;
  feCheck(createVisualizationExporter(exporter, VTK, &numbering, &sol, &mesh, spaces));
  int exportEveryNSteps = 1;
  std::string vtkFileRoot = "output";
  feExportData exportData = {nullptr, exportEveryNSteps, vtkFileRoot};

  TimeIntegrator *solver;
  feTolerances tol{1e-9, 1e-8, 50};
  double t0 = 0.;
  double tEnd = 0.1;
  double nSteps = 50;
  feCheck(createTimeIntegrator(solver, BDF1, tol, system, &numbering, &sol, &mesh, norms, exportData, t0, tEnd, nSteps));
  feCheck(solver->makeStep());

  int nInteriorPlotNodes = 40;
  feBasicViewer viewer("test", mesh.getNbInteriorElems(), nInteriorPlotNodes);

  double xLim[2] = {xa, xb};
  double yLim[2] = {0, 0.2};
  viewer.setAxesLimits(xLim, yLim);

  do{
    viewer.reshapeWindowBox(1.2, mesh, sol);
    viewer.draw1DCurve(mesh, numbering, sol, uDomaine, funSol);
    viewer.windowUpdate();
  } while(!viewer.windowShouldClose());

  petscFinalize();
  return 0;
}