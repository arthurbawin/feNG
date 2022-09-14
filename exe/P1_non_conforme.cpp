#include "feAPI.h"

double fSol(const double t, const std::vector<double> &x, const std::vector<double> par)
{
  double x1 = x[0];
  double y1 = x[1];
  double a = par[1];
  return pow(t, a) * (x1 * x1 + 1) * y1;
  return pow(t, a) * (x1 * x1 + 1) * y1;
}

double fSource(const double t, const std::vector<double> &x, const std::vector<double> par)
{
  double x1 = x[0];
  double y1 = x[1];
  double c1 = par[0];
  double a = par[1];
  double beta = (a == 0.) ? 0. : a * pow(t, a - 1.);
  return -beta * (x1 * x1 + 1) * y1 + c1 * 2 * pow(t, a) * y1;
}

double dfdy(const double t, const std::vector<double> &x, const std::vector<double> par)
{
  double x1 = x[0];
  double y1 = x[1];
  double c1 = par[0];
  double a = par[1];
  return pow(t, a) * (x1 * x1 + 1);
}

double dfdx(const double t, const std::vector<double> &x, const std::vector<double> par)
{
  double x1 = x[0];
  double y1 = x[1];
  double c1 = par[0];
  double a = par[1];
  return pow(t, a) * y1 * (2 * x1);
}

double fSolDot(const double t, const std::vector<double> &x, const std::vector<double> par)
{
  double c1 = par[0];
  double a = par[1];
  double beta = (a == 0.) ? 0. : a * pow(t, a - 1.);
  double x1 = x[0];
  double y1 = x[1];
  return beta * (x1 * x1 + 1) * y1;
}

double f0(const double t, const std::vector<double> &x, const std::vector<double> par) { return 0; }

void fgrad(const double t, const std::vector<double> &x, const std::vector<double> &par,
           std::vector<double> &res)
{
  double x1 = x[0];
  double y1 = x[1];
  double c1 = par[0];
  double a = par[1];
  res[0] = pow(t, a) * y1 * (2 * x1);
  res[1] = pow(t, a) * (x1 * x1 + 1);
}

int main(int argc, char **argv)
{
  // petscInitialize(argc, argv);
  double exposant = 5.;
  double kd = 0.1;
  std::vector<double> par = {kd, exposant};
  feFunction *funSol = new feFunction(fSol, par);
  feFunction *funSolDot = new feFunction(fSolDot, par);
  feFunction *funSource = new feFunction(fSource, par);
  feFunction *dfunSoldx = new feFunction(dfdx, par);
  feFunction *dfunSoldy = new feFunction(dfdy, par);
  feFunction *funZero = new feFunction(f0, par);
  feVectorFunction *fungrad = new feVectorFunction(fgrad, par);

  int nIter = 2;
  std::vector<double> normL2_BDF2(2 * nIter, 0.0);
  std::vector<double> normL2_DC3(2 * nIter, 0.0);
  std::vector<double> normBDF2_bord(2 * nIter, 0.0);
  std::vector<double> normDC3_bord(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);
  std::vector<int> TT;
  TT.resize(nIter);

  for(int iter = 0; iter < nIter; ++iter) {
    // Set the default parameters.
    // const char *meshFile ="../../data/VonKarman/vonKarman1.msh";
    std::string meshFile = "../../data/Random/random" + std::to_string(iter * 2 + 1) + ".msh";
    const char *meshName = meshFile.c_str();
    // Maillage
    int verbosity = 2;
    int order = 1;
    int degreeQuadrature = 5;

    // Create an option parser and parse the command line arguments.
    // If a command line argument is provided, it will overwrite the default parameter.
    // Each "addOption" adds an optional field, although it can be made required by
    // setting the 5th argument of addOption to "true".
    feOptionsParser options(argc, argv);
    options.addOption(&meshName, "-m", "--mesh", "Mesh file");
    options.addOption(&order, "-o", "--order", "Finite element space order");
    options.addOption(&degreeQuadrature, "-dquad", "--degreeQuadrature",
                      "Degree of the quadrature");
    options.addOption(&verbosity, "-v", "--verbosity", "Verbosity level");
    feCheck(options.parse());
    // Set the global verbosity level :
    // - 0 : No information messages, only print warnings and errors
    // - 1 : Moderate information messages (default)
    // - 2 : All information messages
    setVerbose(verbosity);
    // Create a mesh structure from a Gmsh mesh file (version 2.2, 4.1+)
    feMesh2DP1 mesh(meshName);
    nElm[iter] = mesh.getNbInteriorElems();

    int dim;
    feSpace *U_bord, *U_surface, *L_bord, *V_bord;
    feCheck(createFiniteElementSpace(U_bord, &mesh, dim = 1, LINE_CR, order, "U", "Bord",
                                     degreeQuadrature, funSol));
    feCheck(createFiniteElementSpace(U_surface, &mesh, dim = 2, TRI_CR, order, "U", "Surface",
                                     degreeQuadrature, funSol));
    feCheck(createFiniteElementSpace(L_bord, &mesh, dim = 1, LINE_CR, order, "L", "Bord",
                                     degreeQuadrature, funSol));
    feCheck(createFiniteElementSpace(V_bord, &mesh, dim = 1, LINE_CR, order, "V", "Bord",
                                     degreeQuadrature, funSol));

    // Define the set of all finite elements spaces and the set of feSpaces
    // forming the essential boundary conditions. The second set must always be
    // a subset of the first.
    std::vector<feSpace *> spaces = {U_bord, U_surface, L_bord, V_bord};
    std::vector<feSpace *> essentialSpaces = {};

    // Create a "meta" structure containing the numbering of each unknown field.
    feMetaNumber metaNumber(&mesh, spaces, essentialSpaces);

    // Create the solution structure (contains essentially the solution at all DOFs).
    feSolution sol(&mesh, spaces, essentialSpaces, &metaNumber);

    // Define (bi-)linear forms : here we define a bilinear symmetric diffusion form which
    // is the integral of k*(grad phi)_i*(grad phi)_j and a linear form which is the integral
    // of f*phi.
    std::vector<feSpace *> spaceDiffusion2D_U = {U_surface};
    std::vector<feSpace *> spaceSource2D_U = {U_surface};
    std::vector<feSpace *> spaceMasse2D_U = {U_surface};
    std::vector<feSpace *> spaceWeak_bord = {U_bord, L_bord, V_bord};

    feBilinearForm diff_U_M2D(spaceDiffusion2D_U, &mesh, degreeQuadrature,
                              new feSysElm_2D_Diffusion(kd, nullptr));
    feBilinearForm source_U_M2D(spaceSource2D_U, &mesh, degreeQuadrature,
                                new feSysElm_2D_Source(1.0, funSource));
    feBilinearForm masse_U_M2D(spaceMasse2D_U, &mesh, degreeQuadrature,
                               new feSysElm_2D_Masse(1.0, nullptr));
    feBilinearForm weakBC_bord(spaceWeak_bord, &mesh, degreeQuadrature,
                               new feSysElm_1D_weakBC_edo1(1.0, funSolDot));

    // Initialize the linear system. Assembly of the elementary matrices and RHS is
    // performed in the solve step. Two linear solvers are available :
    // MKL Pardiso (direct solver) and PETSc (collection of iterative solvers).
    feLinearSystem *system;
    feCheck(createLinearSystem(system, MKLPARDISO,
                               {&diff_U_M2D, &source_U_M2D, &masse_U_M2D, &weakBC_bord},
                               &metaNumber, &mesh, argc, argv));

    // Norme de la solution
    feComputer *normBDF2 = new feComputer(U_surface, &mesh, &metaNumber, "L2Norm_1Field", funSol);
    feComputer *normDC3 = new feComputer(U_surface, &mesh, &metaNumber, "L2Norm_1Field", funSol);
    feComputer *int_Lambda_BDF2 =
      new feComputer(L_bord, &mesh, &metaNumber, "IntFlux", funSol, fungrad);
    feComputer *int_Lambda_DC3 =
      new feComputer(L_bord, &mesh, &metaNumber, "IntFlux", funSol, fungrad);
    std::vector<feComputer *> comput = {normBDF2, normDC3, int_Lambda_BDF2, int_Lambda_DC3};

    // Create an exporter structure to write the solution for visualization. Currently the only
    // supported format for visualization is VTK.
    feExporter *exporter;
    feCheck(createVisualizationExporter(exporter, VTK, &metaNumber, &sol, &mesh, spaces));
    int exportEveryNSteps = 1;
    std::string vtkFileRoot = "../../data/ResultatVK/VonKarman";
    feExportData exportData = {exporter, exportEveryNSteps, vtkFileRoot};
    // Solve the discrete problem. Initialize a TimeIntegrator object and tolerances on the
    // Newton-Raphson nonlinear solver (tolerance on the solution correction dx, tolerance on the
    // residual, max number of iterations). Here the PDE is linear : the nonlinear solver should
    // converge in 2 iterations. The TimeIntegrator can be STATIONARY, BDF1, BDF2 or a
    // deferred-correction method (DC2F, DC3, DC3F). The solution will be exported for visualization
    // according to the exportData structure. The linear system is assembled and solved in the
    // "makeSteps()" call.
    TimeIntegrator *solver;
    double t0 = 0.;
    double t1 = 1.;
    int nTimeSteps = 5 * pow(2, iter);
    TT[iter] = nTimeSteps;
    feTolerances tol{1e-5, 1e-5, 5};
    feCheck(createTimeIntegrator(solver, DC3, tol, system, &metaNumber, &sol, &mesh, comput,
                                 exportData, t0, t1, nTimeSteps));
    feCheck(solver->makeSteps(nTimeSteps));
    // std::vector<double> &normLu = solver->getNorm(0);
    std::vector<double> &normL2BDF2 = normBDF2->getResult();
    std::vector<double> &normL2DC3 = normDC3->getResult();
    std::vector<double> &norm_BDF2_Lambda_bord = int_Lambda_BDF2->getResult();
    std::vector<double> &norm_DC3_Lambda_bord = int_Lambda_DC3->getResult();
    normL2_BDF2[2 * iter] = *std::max_element(normL2BDF2.begin(), normL2BDF2.end());
    normL2_DC3[2 * iter] = *std::max_element(normL2DC3.begin(), normL2DC3.end());
    normBDF2_bord[2 * iter] =
      *std::max_element(norm_BDF2_Lambda_bord.begin(), norm_BDF2_Lambda_bord.end());
    normDC3_bord[2 * iter] =
      *std::max_element(norm_DC3_Lambda_bord.begin(), norm_DC3_Lambda_bord.end());
    // Free the used memory

    delete solver;
    delete exporter;
    // petscFinalize();
  }

  for(int i = 1; i < nIter; ++i) {
    normL2_BDF2[2 * i + 1] = log(normL2_BDF2[2 * (i - 1)] / normL2_BDF2[2 * i]) / log(2.);
    normL2_DC3[2 * i + 1] = log(normL2_DC3[2 * (i - 1)] / normL2_DC3[2 * i]) / log(2.);
    normBDF2_bord[2 * i + 1] = log(normBDF2_bord[2 * (i - 1)] / normBDF2_bord[2 * i]) / log(2.);
    normDC3_bord[2 * i + 1] = log(normDC3_bord[2 * (i - 1)] / normDC3_bord[2 * i]) / log(2.);
  }
  printf("%12s \t %12s \t %12s \t %12s\t %12s \t %12s\n", "nSteps", "nElm", "||E-BDF||",
         "Taux BDF2", "||E-DC3||", "Taux DC3");
  for(int i = 0; i < nIter; ++i)

    printf("%12d \t %12d \t %12.6e \t %12.6e\t %12.6e \t %12.6e\n", TT[i], nElm[i],
           normL2_BDF2[2 * i], normL2_BDF2[2 * i + 1], normL2_DC3[2 * i], normL2_DC3[2 * i + 1]);
  printf("\n");
  printf("\n");
  printf("Lambda based on BDF2");
  printf("\n");
  printf("%12s \t %12s \t %16s \t %16s \n", "nTimeSteps", "nElm", "||E_lambda||", "Taux lambda");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %16.6e \t %16.6e \n", TT[i], nElm[i], normBDF2_bord[2 * i],
           normBDF2_bord[2 * i + 1]);
  printf("\n");
  printf("\n");
  printf("Lambda based on DC3");
  printf("\n");
  printf("%12s \t %12s \t %16s \t %16s \n", "nTimeSteps", "nElm", "||E_lambda||", "Taux lambda");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %16.6e \t %16.6e \n", TT[i], nElm[i], normDC3_bord[2 * i],
           normDC3_bord[2 * i + 1]);
  return 0;
}