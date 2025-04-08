/*

  Example 1: Solve the steady-state heat equation with constant conductivity k (Poisson equation):

    - div(k grad u) + f = 0  => - k Delta u + f = 0.

    Here we simply take f = -1.

  on a 2D or 3D domain subject to Dirichlet boundary conditions.
  The domain is discretized in space using a simplcial (triangle or tetrahedra) mesh.

*/

// Include all the functions and classes of the library API
#include "feAPI.h"

// This callback is used to evaluate the exact solution u_exact(t,x,y,z),
// with - t             : the time
//      - pos = {x,y,z} : the position vector
//      - par           : the vector of parameters given when creating the associated feFunction.
// 
// Here, the exact solution is the infinite sum :
//
//                 1                  16 
// u_exact(x,y) = ---   sum   ----------------- sin(i pi x) sin(j pi y),
//                 k  i,j odd  pi^4 ij(i^2+j^2)
//
// which is evaluated up to nMax terms.
double uSol_f(const feFunctionArguments &args, const std::vector<double> &par)
{
  double x = args.pos[0];
  double y = args.pos[1];
  double k = par[0];
  // Max number of terms in the infinite sum
  int nMax = 40;
  double Si, Cij, res = 0., pi4 = M_PI*M_PI*M_PI*M_PI;
  for(int i = 1; i < nMax; i += 2) {
    Si = sin(i*M_PI*x);
    for(int j = 1; j < nMax; j += 2) {
      Cij =  16. / (pi4 * i*j * (i*i+j*j));
      res += Cij * Si * sin(j*M_PI*y);
    }
  }
  return 1./k * res;
}

// This callback is used to evaluate the source term f = u_source(t,x,y,z)
// Here it is simply -1.;
double uSource_f(const feFunctionArguments &/*args*/, const std::vector<double> &/* par */)
{
  return -1.;
}

int main(int argc, char **argv)
{
  // Each executable should start with initialize and end with finalize.
  // If compiled with PETSc, they call PetscInitialize() and PetscFinalize()
  // to initialize MPI and PETSc environment.
  // If compiled with MPI but not PETSc, they call MPI_Init() and MPI_Finalize() instead.
  // If compiled with neither PETSc nor MPI, these are dummy functions.
  initialize(argc, argv);

  //
  // Set the default parameters.
  //

  // Space dimension and mesh file
  int dim = 2;
  const char *meshFile = (dim == 2) ? "../data/square1.msh" : "../data/cube2.msh";

  // Verbosity level (amount of information printed in console)
  int verbosity = 1;
  // Degree of the finite element space
  int order = 1;
  // Degree of the polynomial to integrate exactly
  int degreeQuadrature = 10;
  // Choice of linear solver (PETSc, PETSc_MUMPS or Pardiso)
  const char *solverType = "PETSc_MUMPS";
  // For Pardiso, specify if linear system is stored on a single MPI process (-1)
  // or distributed on the processes (>= 0)
  int ownershipSplit = -1;
  // Constant scalar diffusivity parameter
  double k = 1.;

  // Create an option parser with the supported arguments, then parse the command line arguments.
  // Returns an error if a command line argument is ill-formed or if its value is invalid.
  // If a command line argument is provided, it will overwrite the default parameter.
  // Each parameter added with "addOption" is optional by default, to make it required,
  // set the 5th argument of addOption to "true".
  // PETSc options (e.g. -ksp_view) are supported and do not need to be specified manually.
  feOptionsParser options(argc, argv);
  options.addOption(&dim, "-d", "--dim", "Space dimension (2D or 3D)");
  options.addOption(&meshFile, "-m", "--mesh", "Mesh file (gmsh format .msh only)");
  options.addOption(&order, "-o", "--order", "Finite element space order (polynomial degree)");
  options.addOption(&degreeQuadrature, "-dquad", "--degreeQuadrature", "Degree of the quadrature");
  options.addOption(&verbosity, "-v", "--verbosity", "Verbosity level");
  options.addOption(&solverType, "-linsol", "--linear_solver_type", "Choice of linear solver (PETSc or MKL Pardiso)");
  options.addOption(&ownershipSplit, "-own_split", "--ownership_split", "Store linear system on single (-1) or multiple (>= 0) MPI processes");
  options.addOption(&k, "-k", "--diffusivity", "Constant scalar diffusivity");
  feCheck(options.parse());

  if(dim != 2 && dim != 3) {
    feErrorMsg(FE_STATUS_ERROR, "Unsupported space dimension : %d", dim);
    return 1;
  }

  // Set the global verbosity level :
  // - 0 : No information messages, only print warnings and errors
  // - 1 : Moderate information messages (default)
  // - 2 : All information messages
  setVerbose(verbosity);

  // Create a mesh structure from a Gmsh mesh file (version 4.1+)
  feMesh2DP1 mesh(meshFile);

  // Create function objects for the analytic solution, source term and constant parameters.
  // These functions are used to initialize the degrees of freedom, evaluate weak forms
  // and compute e.g. error norms.
  feFunction uSol(uSol_f, {k});
  feFunction uSource(uSource_f, {k});
  // Alternatively, here the source term is constant and equal to -1,
  // so we could use the constant function defined in feAPI.h:
  // feFunction uSource = scalarConstant::minusOne;
  feConstantFunction diffusivity(k);

  // Define a finite element space on each subdomain of the computational domain.
  // The subdomains are defined by the Physical Entities in Gmsh, and the name given
  // for the FE space (i.e. "Bord") must match the name of the Physical Entity in the mesh file.
  // Here the mesh consists of the 2D interior ("Domaine") and the 1D boundary ("Bord"),
  // and we give the name "U" to the solution. Quadrature rules are defined on the finite
  // element space, with parameter "degreeQuadrature" governing the number of quadrature nodes.
  // The feFunction provided is used to initialize the degrees of freedom on the feSpace.
  feSpace *uBord = nullptr, *uDomaine = nullptr;
  feCheck(createFiniteElementSpace(   uBord, &mesh, elementType::LAGRANGE, order, "U",    "Bord", degreeQuadrature, &scalarConstant::zero));
  feCheck(createFiniteElementSpace(uDomaine, &mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &scalarConstant::zero));

  // Define the set of all finite elements spaces and the set of FE spaces
  // forming the essential (Dirichlet) boundary conditions. The second set must always be
  // a subset of the first.
  std::vector<feSpace *> spaces = {uBord, uDomaine};
  std::vector<feSpace *> essentialSpaces = {uBord};

  // Create the DOFs numbering for each field. Here there is only one field U.
  feMetaNumber numbering(&mesh, spaces, essentialSpaces);

  // Create the solution structure (it is essentially a wrapper containing the solution at each DOF).
  feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);

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
  // There is no form to define on the boundary as all DOFs are essential.
  feBilinearForm *diff = nullptr, *source = nullptr;
  if(dim == 2)
    feCheck(createBilinearForm(diff, {uDomaine}, new feSysElm_Diffusion<2>(&diffusivity)));
  if(dim == 3)
    feCheck(createBilinearForm(diff, {uDomaine}, new feSysElm_Diffusion<3>(&diffusivity)));
  feCheck(createBilinearForm(source, {uDomaine}, new feSysElm_Source(&scalarConstant::minusOne)));
  std::vector<feBilinearForm*> forms = {diff, source};

  // Create the linear system. Assembly of the elementary matrices and RHS is
  // performed in the solve step below ("makeSteps").
  feLinearSystem *system;
  if(std::string(solverType).compare("PETSc") == 0) {
    feCheck(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
  } else if(std::string(solverType).compare("PETSc_MUMPS") == 0) {
    feCheck(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
  } else if(std::string(solverType).compare("Pardiso") == 0) {
    feCheck(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns(), ownershipSplit));
  } else {
    feErrorMsg(FE_STATUS_ERROR, "Unrecognized linear solver type: %s", solverType);
    return 1;
  }

  // The global matrix and RHS can be displayed at each Newton-Raphson iteration
  // by uncommenting the following :
  // system->setDisplayMatrixInConsole(true);
  // system->setDisplayRHSInConsole(true);
  // system->setDisplayMatrixInWindow(true);

  // Post-processing tools to compute norms and whatnot
  // These are given to the solver and are performed after each time step
  // They can also be called here before or after resolution
  feNorm *errorU;
  feCheck(createNorm(errorU, L2_ERROR, {uDomaine}, &sol, &uSol));
  std::vector<feNorm*> norms = {errorU};

  // Create an exporter structure to write the solution for visualization.
  // Currently the only supported format for visualization is the legacy VTK format.
  feExporter *exporter;
  feCheck(createVisualizationExporter(exporter, VTK, &numbering, &sol, &mesh, spaces));
  std::string vtkFileRoot = "example1";
  feExportData exportData = {exporter, 1, vtkFileRoot};

  // Set tolerances on the nonlinear solver
  feNLSolverOptions NLoptions{
    1e-10, // Tolerance on the nonlinear residual NL(u)
    1e-10, // Tolerance on the correction du
    1e4,   // Divergence tolerance : solver stops if ||NL(u)|| > tol at any point
    50,    // Maximum number of iteration
    // For unsteady problems only :
    4,     // Recompute jacobian matrix after this many steps
    1e-1   // Recompute jacobian matrix if ||NL(u)||_i > val * ||NL(u)||_i-1
  };

  // Solve the discrete problem.
  // Create a TimeIntegrator object, which will solve the nonlinear problem NL(u) = 0 at each time step.
  // Here, the problem is steady and the PDE is linear in u, hence
  // we compute a single "timestep" and the solution should converge in a single Newton-Raphson iteration.
  // TimeIntegrator can be STATIONARY, BDF1 or BDF2 (experimental deferred correction methods are
  // begin reworked and currently unavailable (DC2F, DC3, DC3F)). 
  // The solution will be exported for visualization according to the exportData structure.
  // The linear system is assembled and solved during the "makeSteps(1)" call.
  TimeIntegrator *solver;
  feCheck(createTimeIntegrator(solver, timeIntegratorScheme::STATIONARY, NLoptions, system, &sol, &mesh, norms, exportData));
  feCheck(solver->makeSteps(1));

  // Evaluate the error in L2 norm
  feInfo("L2 error = %10.10f", errorU->compute());

  // Free the used memory (all pointers)
  delete errorU;
  delete solver;
  delete exporter;
  delete system;
  for(feBilinearForm* f : forms)
    delete f;
  for(feSpace *s : spaces)
    delete s;

  finalize();
  return 0;
}