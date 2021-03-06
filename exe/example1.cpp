
#include "feAPI.h"

double fSol(const double t, const std::vector<double> x, const std::vector<double> par)
{
  return pow(x[0], 6);
}

double fSource(const double t, const std::vector<double> x, const std::vector<double> par)
{
  double k = par[0];
  return k * 30. * pow(x[0], 4);
}

int main(int argc, char **argv)
{
  // petscInitialize(argc, argv);

  // Set the default parameters.
  const char *meshFile = "square.msh";
  int verbosity = 2;
  int order = 1;
  int degreeQuadrature = 10;

  // Create an option parser and parse the command line arguments.
  // If a command line argument is provided, it will overwrite the default parameter.
  // Each "addOption" adds an optional field, although it can be made required by
  // setting the 5th argument of addOption to "true".
  feOptionsParser options(argc, argv);
  options.addOption(&meshFile, "-m", "--mesh", "Mesh file");
  options.addOption(&order, "-o", "--order", "Finite element space order");
  options.addOption(&degreeQuadrature, "-dquad", "--degreeQuadrature", "Degree of the quadrature");
  options.addOption(&verbosity, "-v", "--verbosity", "Verbosity level");
  feCheck(options.parse());

  // Set the global verbosity level :
  // - 0 : No information messages, only print warnings and errors
  // - 1 : Moderate information messages (default)
  // - 2 : All information messages
  setVerbose(verbosity);

  // Create a mesh structure from a Gmsh mesh file (version 2.2, 4.1+)
  feMesh2DP1 mesh(meshFile);

  // Create function structures for the analytic solution and the source term
  // Here the model PDE is the stationary heat equation div(k grad u) + f = 0,
  // with k the thermal conductivity.
  double k = 1.0;
  feFunction *funSol = new feFunction(fSol, {});
  feFunction *funSource = new feFunction(fSource, {k});

  // Define a finite element space on each subdomain of the computational domain.
  // Here the mesh contains only a 2D subdomain ("Domaine") and a 1D boundary ("Bord"),
  // and the solution field is named "U". Quadrature rules are defined on the finite
  // element space, with parameter "degreeQuadrature" governing the number of quad nodes.
  // The function provided is used to initialize the degrees of freedom on the feSpace.
  int dim;
  feSpace *uBord, *uDomaine;
  feCheck(createFiniteElementSpace(uBord, &mesh, dim = 1, LINE, order, "U", "Bord",
                                   degreeQuadrature, funSol));
  feCheck(createFiniteElementSpace(uDomaine, &mesh, dim = 2, TRI, order, "U", "Domaine",
                                   degreeQuadrature, funSol));

  // Define the set of all finite elements spaces and the set of feSpaces
  // forming the essential boundary conditions. The second set must always be
  // a subset of the first.
  std::vector<feSpace *> spaces = {uBord, uDomaine};
  std::vector<feSpace *> essentialSpaces = {uBord};

  // Create a "meta" structure containing the numbering of each unknown field.
  feMetaNumber metaNumber(&mesh, spaces, essentialSpaces);

  // Create the solution structure (contains essentially the solution at all DOFs).
  feSolution sol(&mesh, spaces, essentialSpaces, &metaNumber);

  // Define (bi-)linear forms : here we define a bilinear symmetric diffusion form which
  // is the integral of k*(grad phi)_i*(grad phi)_j and a linear form which is the integral
  // of f*phi.
  feBilinearForm diffU({uDomaine}, &mesh, degreeQuadrature, new feSysElm_2D_Diffusion(k, nullptr));
  feBilinearForm sourceU({uDomaine}, &mesh, degreeQuadrature,
                         new feSysElm_2D_Source(1.0, funSource));

  // Initialize the linear system. Assembly of the elementary matrices and RHS is
  // performed in the solve step. Two linear solvers are available :
  // MKL Pardiso (direct solver) and PETSc (collection of iterative solvers).
  feLinearSystem *system;
  feCheck(
    createLinearSystem(system, MKLPARDISO, {&diffU, &sourceU}, &metaNumber, &mesh, argc, argv));

  // Define post-processing tools to compute norms and whatnot (norms will be replaced by
  // fePostProc)
  feNorm normU(uDomaine, &mesh, degreeQuadrature, funSol);
  std::vector<feNorm *> norms = {&normU};

  // Create an exporter structure to write the solution for visualization. Currently the only
  // supported format for visualization is VTK.
  feExporter *exporter;
  feCheck(createVisualizationExporter(exporter, VTK, &metaNumber, &sol, &mesh, spaces));
  int exportEveryNSteps = 1;
  std::string vtkFileRoot = "root";
  feExportData exportData = {exporter, exportEveryNSteps, vtkFileRoot};

  // Solve the discrete problem. Initialize a TimeIntegrator object and tolerances on the
  // Newton-Raphson nonlinear solver (tolerance on the solution correction dx, tolerance on the
  // residual, max number of iterations). Here the PDE is linear : the nonlinear solver should
  // converge in 2 iterations. The TimeIntegrator can be STATIONARY, BDF1, BDF2 or a
  // deferred-correction method (DC2F, DC3, DC3F). The solution will be exported for visualization
  // according to the exportData structure. The linear system is assembled and solved in the
  // "makeSteps()" call.
  TimeIntegrator *solver;
  feTolerances tol{1e-9, 1e-8, 3};
  feCheck(createTimeIntegrator(solver, STATIONARY, tol, system, &metaNumber, &sol, &mesh, norms,
                               exportData));
  feCheck(solver->makeSteps(0));

  // Post-process the solution
  fePostProc post(uDomaine, &mesh, &metaNumber, funSol);
  feInfo("Measure = %f", post.computeMeasure());
  feInfo("Integral of sol = %f", post.computeSolutionIntegral(&sol));
  feInfo("Integral of fun = %f", post.computeFunctionIntegral(funSol, 0.0));
  feInfo("L2 Error = %10.10f", post.computeL2ErrorNorm(&sol));

  // Free the used memory
  delete solver;
  delete exporter;
  delete uBord;
  delete uDomaine;
  delete funSol;
  delete funSource;

  // petscFinalize();

  return 0;
}