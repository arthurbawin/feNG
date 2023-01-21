/*

  Example 1: Solve the steady-state heat equation with constant conductivity k (Poisson equation):

    div(k grad u) + f = 0  => - k Delta u = f

  on a 2D domain subject to Dirichlet boundary conditions.
  The domain is discretized in space using a simplcial (triangle) mesh.
  
*/

#include "feAPI.h"

double fSol(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  double x = pos[0];
  // return pow(x, 6);
  double y = pos[1];
  double a = 10.0;
  double b = 2.;
  return 0.5*(1.0 + tanh(a*(x - sin(b*M_PI*y)/4.0)));
}

double fSource(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  double x = pos[0];
  double y = pos[1];
  double a = 10.0;
  double b = 2.;
  double k = par[0];
  // return k * 30. * pos[0] * pos[0] * pos[0] * pos[0];

  double T = pow(tanh(a*(x - sin(M_PI*b*y)/4)),2);
  double uxx = a*a*tanh(a*(x - sin(M_PI*b*y)/4))*(T - 1);
  double uyy = (a*a*b*b*M_PI*M_PI*tanh(a*(x - sin(M_PI*b*y)/4))*pow(cos(M_PI*b*y),2)*(T - 1))/16 - (a*b*b*M_PI*M_PI*sin(M_PI*b*y)*(T - 1))/8;
  return k * (uxx + uyy);
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  return 0.0;
}

double fConstant(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  return par[0];
}

int main(int argc, char **argv)
{
  // The executable should start with petscInitialize and end with petscFinalize.
  // If the library was not compiled with PETSc, these functions are dummy.
  petscInitialize(argc, argv);

  // Set the default parameters.
  const char *meshFile = "../data/square1.msh";
  int verbosity = 2;
  int order = 1;
  int degreeQuadrature = 10;

  // Create an option parser and parse the command line arguments.
  // Returns an error if a command line argument is ill-formed.
  // If a command line argument is provided, it will overwrite the default parameter.
  // Each parameter added with "addOption" is optional by default, to make it required,
  // set the 5th argument of addOption to "true".
  feOptionsParser options(argc, argv);
  options.addOption(&meshFile, "-m", "--mesh", "Mesh file (gmsh format .msh only)");
  options.addOption(&order, "-o", "--order", "Finite element space order (polynomial degree)");
  options.addOption(&degreeQuadrature, "-dquad", "--degreeQuadrature", "Degree of the quadrature");
  options.addOption(&verbosity, "-v", "--verbosity", "Verbosity level");
  feCheck(options.parse());

  // Set the global verbosity level :
  // - 0 : No information messages, only print warnings and errors
  // - 1 : Moderate information messages (default)
  // - 2 : All information messages
  setVerbose(verbosity);

  // Create a mesh structure from a Gmsh mesh file (version 2.2 or 4.1+)
  feMesh2DP1 mesh(meshFile);

  // Create function objects for the analytic solution and the source term.
  // Here the model PDE is the stationary heat equation div(k grad u) + f = 0,
  // with k the thermal conductivity.
  // These functions are used to initialize the degrees of freedom (unknown and
  // boundary conditions).
  double k = 1.0;
  feFunction *funSol    = new feFunction(fSol,    {});
  feFunction *funSource = new feFunction(fSource, {k});
  feFunction *funZero   = new feFunction(fZero,   {});
  feFunction *kDiffusivity = new feFunction(fConstant, {1.0});

  // Define a finite element space on each subdomain of the computational domain.
  // The subdomains are defined by the Physical Entities in gmsh, and the name given
  // for the FE space (i.e. "Bord") must match the name of the Physical Entity in the mesh file.
  // Here the mesh consists of the 2D interior ("Domaine") and the 1D boundary ("Bord"),
  // and we give the name "U" to the solution. Quadrature rules are defined on the finite
  // element space, with parameter "degreeQuadrature" governing the number of quadrature nodes.
  // The feFunction provided is used to initialize the degrees of freedom on the feSpace.
  int dim;
  feSpace *uBord, *uDomaine;
  feCheck(createFiniteElementSpace(uBord, &mesh, dim = 1, LINE, LAGRANGE, order, "U", "Bord", degreeQuadrature, funSol));
  feCheck(createFiniteElementSpace(uDomaine, &mesh, dim = 2, TRI, LAGRANGE, order, "U", "Domaine", degreeQuadrature, funZero));

  // Define the set of all finite elements spaces and the set of FE spaces
  // forming the essential (Dirichlet) boundary conditions. The second set must always be
  // a subset of the first.
  std::vector<feSpace*> spaces = {uBord, uDomaine};
  std::vector<feSpace*> essentialSpaces = {uBord};

  // Create the DOFs numbering for each field. Here there is only one field U.
  feMetaNumber numbering(&mesh, spaces, essentialSpaces);

  // Create the solution structure (it is essentially a wrapper containing the solution at each DOF).
  feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);

  // Define the (bi-)linear forms: here we define on the interior a bilinear symmetric "diffusion" form:
  // 
  //                         /
  // feSysElm_2D_Diffusion = | k*(grad phi)_i*(grad phi)_j dx
  //                         /
  //
  // and a linear "source" form:
  //
  //                      /
  // feSysElm_2D_Source = | f*phi_i dx
  //                      /
  //
  // There is no form to define on the boundary.
  feBilinearForm *diff, *source;
  feCheck(createBilinearForm(  diff, {uDomaine}, new feSysElm_2D_Diffusion(k, nullptr)  ));
  feCheck(createBilinearForm(source, {uDomaine}, new feSysElm_2D_Source(1.0, funSource) ));

  // Create the linear system. Assembly of the elementary matrices and RHS is
  // performed in the solve step below ("makeSteps"). Two linear solvers are available:
  // MKL Pardiso (direct solver) and PETSc (collection of iterative solvers).
  feLinearSystem *system;
  feCheck(createLinearSystem(system, PETSC, {diff, source}, numbering.getNbUnknowns(), argc, argv));

  // Post-processing tools to compute norms and whatnot
  feNorm normU(L2_ERROR, {uDomaine}, &sol, funSol);
  std::vector<feNorm *> norms = {&normU};

  // Create an exporter structure to write the solution for visualization. Currently the only
  // supported format for visualization is VTK.
  feExporter *exporter;
  feCheck(createVisualizationExporter(exporter, VTK, &numbering, &sol, &mesh, spaces));
  int exportEveryNSteps = 1;
  std::string vtkFileRoot = "output";
  feExportData exportData = {exporter, exportEveryNSteps, vtkFileRoot};

  // Solve the discrete problem. Initialize a TimeIntegrator object and tolerances on the
  // Newton-Raphson nonlinear solver (tolerance on the solution correction dx, tolerance on the
  // residual, max number of iterations). Here the PDE is linear in u: the nonlinear solver should
  // converge in 2 iterations. The TimeIntegrator can be STATIONARY, BDF1, BDF2 or a
  // deferred-correction method (DC2F, DC3, DC3F, still experimental). The solution will be 
  // exported for visualization according to the exportData structure. The linear system is 
  // assembled and solved in the "makeSteps()" call.
  TimeIntegrator *solver;
  feTolerances tol{1e-9, 1e-8, 10};
  feCheck(createTimeIntegrator(solver, STATIONARY, tol, system, &numbering, &sol, &mesh, norms, exportData));
  feCheck(solver->makeSteps(0));

  // Compute L2 norm of the error u-uh
  feNorm *norm;
  feCheck(createNorm(norm, L2_ERROR, {uDomaine}, &sol, funSol));

  fePostProc post(uDomaine, &mesh, &numbering, funSol);
  feInfo("L2 error = %10.10f", norm->compute(L2_ERROR));
  feInfo("L2 error = %10.10f", post.computeL2ErrorNorm(&sol));

  feInfo("L1   error = %10.10f", norm->compute(L1_ERROR));
  feInfo("Linf error = %10.10f", norm->compute(LINF_ERROR));

  feInfo("Area = %10.10f", norm->compute(AREA));
  feInfo("Area = %10.10f", post.computeMeasure());

  feInfo("Integral = %10.10f", norm->compute(INTEGRAL));
  feInfo("Integral = %10.10f", post.computeSolutionIntegral(&sol));

  // Free the used memory
  delete diff;
  delete source;
  delete solver;
  delete exporter;
  delete uDomaine;
  delete uBord;
  delete funZero;
  delete funSol;
  delete funSource;
  delete norm;

  petscFinalize();
  return 0;
}