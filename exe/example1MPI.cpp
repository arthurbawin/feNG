/*

  Example 1: Solve the steady-state heat equation with constant conductivity k (Poisson equation):

    div(k grad u) + f = 0  => - k Delta u = f

  on a 2D domain subject to Dirichlet boundary conditions.
  The domain is discretized in space using a simplcial (triangle) mesh.

*/

// #include "mpi.h"
#include "feAPI.h"
#if defined(HAVE_PETSC)
  #include "petscsys.h"
#endif

double fSol(const double /* t */, const std::vector<double> &pos, const std::vector<double>& /* par */)
{
  double x = pos[0];
  return pow(x, 6);
}

double fSource(const double /* t */, const std::vector<double> &pos, const std::vector<double>& par)
{
  double x = pos[0];
  double k = par[0];
  return -k * 30. * x * x * x * x;
}

int main(int argc, char **argv)
{
  // The executable should start with petscInitialize and end with petscFinalize.
  // If the library was not compiled with PETSc, these functions are dummy.
  initialize(argc, argv);

  // Set the default parameters.

  // int dim = 2;
  const char *meshFile = "../data/square1.msh";

  // int dim = 3;
  // const char *meshFile = "../data/cube4.msh";
  // const char *meshFile = "../data/singleTet.msh";

  int verbosity = 2;
  int order = 1;
  int degreeQuadrature = 10;

  const char *solverType = "PETSc_MUMPS";
  // Decide if the Pardiso linear system is stored on a single process (-1)
  // or distributed on the processes (>= 0)
  int ownershipSplit = -1;

#if defined(HAVE_PETSC)
  char meshFileTmp[64];
  PetscBool permute;
  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "Poisson example options", "");
  {
    PetscCall(PetscOptionsString("-m", "Mesh file", "", meshFileTmp, meshFileTmp, sizeof(meshFileTmp), NULL));
    PetscCall(PetscOptionsInt("-o", "Finite element space order (polynomial degree)", "", order, &order, NULL));
    permute  = PETSC_FALSE;
    PetscCall(PetscOptionsBool("-permute", "Reorder the matrix?", "", permute, &permute, NULL));
  }
  PetscOptionsEnd();
  meshFile = &meshFileTmp[0];
#else
  // Create an option parser and parse the command line arguments.
  // The parser is not compatible with the PETSc options for now.
  // Returns an error if a command line argument is ill-formed.
  // If a command line argument is provided, it will overwrite the default parameter.
  // Each parameter added with "addOption" is optional by default, to make it required,
  // set the 5th argument of addOption to "true".
  feOptionsParser options(argc, argv);
  options.addOption(&meshFile, "-m", "--mesh", "Mesh file (gmsh format .msh only)");
  options.addOption(&order, "-o", "--order", "Finite element space order (polynomial degree)");
  options.addOption(&degreeQuadrature, "-dquad", "--degreeQuadrature", "Degree of the quadrature");
  options.addOption(&verbosity, "-v", "--verbosity", "Verbosity level");
  options.addOption(&solverType, "-linsol", "--linear_solver_type", "Choice of linear solver (PETSc or MKL Pardiso)");
  options.addOption(&ownershipSplit, "-own_split", "--ownership_split", "Store linear system on single (-1) or multiple (>= 0) MPI processes");
  feCheck(options.parse());
#endif

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
  feFunction funSol(fSol);
  feFunction funSource(fSource, {k});
  feConstantFunction funZero(0.);
  feConstantFunction funOne(1.);
  feConstantFunction diffusivity(k);

  // Define a finite element space on each subdomain of the computational domain.
  // The subdomains are defined by the Physical Entities in gmsh, and the name given
  // for the FE space (i.e. "Bord") must match the name of the Physical Entity in the mesh file.
  // Here the mesh consists of the 2D interior ("Domaine") and the 1D boundary ("Bord"),
  // and we give the name "U" to the solution. Quadrature rules are defined on the finite
  // element space, with parameter "degreeQuadrature" governing the number of quadrature nodes.
  // The feFunction provided is used to initialize the degrees of freedom on the feSpace.
  feSpace *uBord, *uDomaine;
  feCheck(createFiniteElementSpace(uBord, &mesh, elementType::LAGRANGE, order, "U", "Bord",
                                   degreeQuadrature, &funSol));
  feCheck(createFiniteElementSpace(uDomaine, &mesh, elementType::LAGRANGE, order, "U", "Domaine",
                                   degreeQuadrature, &funZero));

  // feNorm nn(mesh.getCncGeoByName("Domaine"), &funSol);
  // feInfo("Volume = %10.10f", nn.compute(AREA));
  // feNorm mm(mesh.getCncGeoByName("Bord"), &funSol);
  // feInfo("Area of boundary = %10.10f", mm.compute(AREA));

  // Define the set of all finite elements spaces and the set of FE spaces
  // forming the essential (Dirichlet) boundary conditions. The second set must always be
  // a subset of the first.
  std::vector<feSpace *> spaces = {uBord, uDomaine};
  std::vector<feSpace *> essentialSpaces = {uBord};

  // Create the DOFs numbering for each field. Here there is only one field U.
  feMetaNumber numbering(&mesh, spaces, essentialSpaces);

  // Create the solution structure (it is essentially a wrapper containing the solution at each
  // DOF).
  feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);

  //////////////////////
  // sol.initializeUnknowns(&mesh);
  // sol.initializeEssentialBC(&mesh);

  // feExporter *exporter;
  // feCheck(createVisualizationExporter(exporter, VTK, &numbering, &sol, &mesh, spaces));
  // feCheck(exporter->writeStep("tet.vtk"));

  // finalize();
  // return 0;
  //////////////////////

  // Define the (bi-)linear forms: here we define on the interior a bilinear symmetric "diffusion"
  // form:
  //
  //                         /
  //    feSysElm_Diffusion = | k*(grad phi)_i*(grad phi)_j dx
  //                         /
  //
  // and a linear "source" form:
  //
  //                      /
  //    feSysElm_Source = | f*phi_i dx
  //                      /
  //
  // There is no form to define on the boundary.
  feBilinearForm *diff, *source;
  feCheck(createBilinearForm(diff, {uDomaine}, new feSysElm_Diffusion<2>(&diffusivity)));
  feCheck(createBilinearForm(source, {uDomaine}, new feSysElm_Source(&funSource)));

  // Create the linear system. Assembly of the elementary matrices and RHS is
  // performed in the solve step below ("makeSteps"). Two linear solvers are available:
  // MKL Pardiso (direct solver) and PETSc (collection of iterative solvers).
  feLinearSystem *system;
  if(std::string(solverType).compare("PETSc") == 0) {
    feCheck(createLinearSystem(system, PETSC, {diff, source}, numbering.getNbUnknowns()));
  } else if(std::string(solverType).compare("PETSc_MUMPS") == 0) {
    feCheck(createLinearSystem(system, PETSC_MUMPS, {diff, source}, numbering.getNbUnknowns()));
  } else if(std::string(solverType).compare("Pardiso") == 0) {
    feCheck(createLinearSystem(system, MKLPARDISO, {diff, source}, numbering.getNbUnknowns(), ownershipSplit));
  } else {
    feErrorMsg(FE_STATUS_ERROR, "Unrecognized linear solver type: %s", solverType);
    return 1;
  }

  // system->setDisplayRHSInConsole(true);

#if defined(HAVE_PETSC)
  system->setRelativeTol(1e-10);
  system->setReorderingStatus(permute);
#endif

  // Post-processing tools to compute norms and whatnot
  feNorm normU(L2_ERROR, {uDomaine}, &sol, &funOne);
  std::vector<feNorm *> norms = {&normU};

  // Create an exporter structure to write the solution for visualization. Currently the only
  // supported format for visualization is VTK.
  feExporter *exporter;
  feCheck(createVisualizationExporter(exporter, VTK, &numbering, &sol, &mesh, spaces));
  std::string vtkFileRoot = "output";
  feExportData exportData = {exporter, 1, vtkFileRoot};
  // feExportData exportData = {nullptr, 1, ""};

  // Solve the discrete problem. Initialize a TimeIntegrator object and tolerances on the
  // Newton-Raphson nonlinear solver (tolerance on the solution correction dx, tolerance on the
  // residual, max number of iterations). Here the PDE is linear in u: the nonlinear solver should
  // converge in 2 iterations. The TimeIntegrator can be STATIONARY, BDF1, BDF2 or a
  // deferred-correction method (DC2F, DC3, DC3F, still experimental). The solution will be
  // exported for visualization according to the exportData structure. The linear system is
  // assembled and solved in the "makeSteps()" call.
  TimeIntegrator *solver;
  feTolerances tol{1e-9, 1e-8, 1e4, 4};
  feCheck(createTimeIntegrator(solver, STATIONARY, tol, system, &numbering, &sol, &mesh, norms,
                               exportData));
  feCheck(solver->makeStep());

  // Compute L2 norm of the error u-uh
  feNorm *norm;
  feCheck(createNorm(norm, L2_ERROR, {uDomaine}, &sol, &funSol));
  feInfo("L2 error = %10.10f", norm->compute(L2_ERROR));

  // Free the used memory
  delete diff;
  delete source;
  delete solver;
  delete exporter;
  delete uDomaine;
  delete uBord;
  // delete norm;
  delete system;

  finalize();
  return 0;
}