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

void fVelocity(const double t, const std::vector<double> &pos, const std::vector<double> &par, std::vector<double> &res)
{
  // Constant x-velocity equal to the given parameter
  res[0] = par[0];
  res[1] = 0.;
}

double fConstant(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  return par[0];
}

int main(int argc, char **argv)
{
  petscInitialize(argc, argv);

  // int nElm = 10;

  double ceq = 1.0;
  double keq = 1.0;
  double feq = 1.0;

  feOptionsParser options(argc, argv);
  // options.addOption(&nElm, "-n", "--nElm", "Number of mesh elements");
  options.addOption(&ceq, "-c", "--velocity", "Constant velocity");
  options.addOption(&feq, "-f", "--source", "Constant source term");
  feCheck(options.parse());

  double ccode = ceq;
  double kcode = -keq;
  double fcode = -feq;

  feFunction *funSol      = new feFunction(fSol, {ceq, keq});

  feFunction *funSource   = new feFunction(fConstant, {fcode});
  feFunction *funZero     = new feFunction(fConstant, {0.});
  feFunction *funOne     = new feFunction(fConstant, {1.});
  feFunction *diffusivity = new feFunction(fConstant, {kcode});
  feFunction *scalarVelocity = new feFunction(fConstant, {ccode});
  feVectorFunction *velocity = new feVectorFunction(fVelocity, {ccode});

  int deg = 1, degreeQuadrature = 4;
  double xa = 0.;
  double xb = 1.;

  int nConv = 6;
  int nElm[8] = {10, 20, 40, 80, 160, 320, 640, 1280};
  std::vector<double> L2error(2 * nConv, 0.0);

  for(int i = 0; i < nConv; ++i){

    feInfo("Element Peclet = %f", fabs(ceq)*(xb-xa)/nElm[i] / (2.*fabs(keq)) );
    feInfo("Global  Peclet = %f", fabs(ceq)*(xb-xa) / (fabs(keq)) );

    feMesh1DP1 mesh(xa, xb, nElm[i], "BXA", "BXB", "Domaine");

    feSpace *uG, *uD, *uDomaine;
    feCheck(createFiniteElementSpace(uG, &mesh, elementType::LAGRANGE, deg, "U", "BXA", degreeQuadrature, funZero));
    feCheck(createFiniteElementSpace(uD, &mesh, elementType::LAGRANGE, deg, "U", "BXB", degreeQuadrature, funZero));
    feCheck(createFiniteElementSpace(uDomaine, &mesh, elementType::LAGRANGE, deg, "U", "Domaine", degreeQuadrature, funZero));

    std::vector<feSpace*> spaces = {uG, uD, uDomaine};
    std::vector<feSpace*> essentialSpaces = {uG, uD};

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);

    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
    sol.initializeUnknowns(&mesh);

    feBilinearForm *adv, *diff, *source, *supg;
    feCheck(createBilinearForm(    adv, {uDomaine}, new feSysElm_Advection<1>(velocity))    );
    feCheck(createBilinearForm(   diff, {uDomaine}, new feSysElm_Diffusion<1>(diffusivity)) );
    feCheck(createBilinearForm( source, {uDomaine}, new feSysElm_Source(funSource))          );
    feCheck(createBilinearForm(   supg, {uDomaine}, new feSysElm_1D_SUPGStab(-kcode, scalarVelocity, funSource)) );

    feLinearSystem *system;
    feCheck(createLinearSystem(system, PETSC, {adv, diff, source, supg}, numbering.getNbUnknowns(), argc, argv));

    feNorm normU(L2_ERROR, {uDomaine}, &sol, funSol);
    std::vector<feNorm *> norms = {&normU};

    feExporter *exporter;
    feCheck(createVisualizationExporter(exporter, VTK, &numbering, &sol, &mesh, spaces));
    int exportEveryNSteps = 1;
    std::string vtkFileRoot = "output";
    feExportData exportData = {nullptr, exportEveryNSteps, vtkFileRoot};

    TimeIntegrator *solver;
    feTolerances tol{1e-10, 1e-10, 10};
    feCheck(createTimeIntegrator(solver, STATIONARY, tol, system, &numbering, &sol, &mesh, norms, exportData));
    feCheck(solver->makeStep());

    double umax = 0.;
    for(auto val : sol.getSolutionReference())
      umax = fmax(umax, fabs(val));
    feInfo("umax = %f", umax);
    feInfo("umax = %f", umax);
    feInfo("umax = %f", umax);
    feInfo("umax = %f", umax);

    // Get the error
    L2error[2*i] = normU.compute();

    // if(i == nConv-1){
      int nInteriorPlotNodes = 40;
      feBasicViewer viewer("test", mesh.getNumInteriorElements(), nInteriorPlotNodes);

      // double xLim[2] = {xa, xb};
      // double yLim[2] = {0, 0.2};
      // viewer.setAxesLimits(xLim, yLim);

      do{
        viewer.reshapeWindowBox(0, 1., -0.1, 1.3);
        viewer.draw1DCurve(mesh, numbering, sol, uDomaine, funSol, ceq);
        viewer.windowUpdate();
      } while(!viewer.windowShouldClose());
    // }
  }

  // Compute the convergence rate
  for(int i = 1; i < nConv; ++i) {
    L2error[2 * i + 1] = -log(L2error[2 * i] / L2error[2 * (i - 1)]) / log(nElm[i] / nElm[i - 1]);
  }
  printf("%12s \t %12s \t %12s \n", "nElm", "||E_U||", "rate");
  for(int i = 0; i < nConv; ++i){
    printf("%12d \t %12.6e \t %12.6e\n", nElm[i], L2error[2 * i], L2error[2 * i + 1]);
  }

  petscFinalize();
  return 0;
}