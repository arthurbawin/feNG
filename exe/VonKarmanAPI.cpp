#include "feAPI.h"

static double h = 16.0;
static double cx = h;
static double cy = h / 2.;
static double r = 0.5;
static double d = 2 * r;
static double U = 1;
static double rho = 1.;
static double Re = 200;

double fSolU_inlet(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  double y = pos[1];
  double coeff = 1.;
  // if (t<2.) coeff = pow(t,2)/4.;
  // else coeff=1.;
  return ((0.125 * y - 0.5) + std::min(t / 10, 1.) * (-0.125 * y + 0.5) + 1) * coeff;
}

double fSolU_inlet_dot(const double t, const std::vector<double> &pos,
                       const std::vector<double> &par)
{
  double y = pos[1];
  if(t < 10.)
    (-0.125 * y + 0.5) / 10.;
  else
    0.;
}

double fSolP(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  double x = pos[0];
  double y = pos[1];
  return 0;
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  return 0.01;
}

int main(int argc, char **argv)
{
  // petscInitialize(argc, argv);
  int nIter = 1;
  std::vector<double> intL_u(2 * nIter, 0.0);
  std::vector<double> intL_v(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);
  std::vector<int> TT;
  TT.resize(nIter);
  for(int iter = 0; iter < nIter; ++iter) {
    // Set the default parameters.
    const char *meshFile ="../../data/VonKarman/vonKarman1_V3.msh";
    // const char *meshFile = "../../data/VonKarman/test.msh";
    int verbosity = 2;
    int order = 2;
    int degreeQuadrature = 5;

    // Create an option parser and parse the command line arguments.
    // If a command line argument is provided, it will overwrite the default parameter.
    // Each "addOption" adds an optional field, although it can be made required by
    // setting the 5th argument of addOption to "true".
    feOptionsParser options(argc, argv);
    options.addOption(&meshFile, "-m", "--mesh", "Mesh file");
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
    feMesh2DP1 mesh(meshFile);
    nElm[iter] = mesh.getNbInteriorElems();

    // Create function structures for the analytic solution and the source term
    // Here the model PDE is the stationary heat equation div(k grad u) + f = 0,
    // with k the thermal conductivity.
    double mu = rho * U * d / Re;
    std::vector<double> stokesParam = {rho, mu};
    feFunction *funSolU_inlet = new feFunction(fSolU_inlet, stokesParam);
    feFunction *funSolU_inlet_dot = new feFunction(fSolU_inlet_dot, stokesParam);
    feFunction *funSolP = new feFunction(fSolP, stokesParam);
    feFunction *funZero = new feFunction(fZero, {});

    // Define a finite element space on each subdomain of the computational domain. Quadrature rules
    // are defined on the finite element space, with parameter "degreeQuadrature" governing the
    // number of quad nodes. The function provided is used to initialize the degrees of freedom on
    // the feSpace.
    int dim;
    feSpace *U_walls, *U_inlet, *U_cylindre, *U_surface, *V_walls, *V_inlet, *V_cylindre,
      *V_surface, *P_surface, *Lu_cylinder, *Lv_cylinder, *Bu_cylinder, *Bv_cylinder;
    feCheck(createFiniteElementSpace(U_walls, &mesh, dim = 1, LINE, order, "U", "walls",
                                     degreeQuadrature, funZero));
    feCheck(createFiniteElementSpace(U_inlet, &mesh, dim = 1, LINE, order, "U", "inlet",
                                     degreeQuadrature, funSolU_inlet));
    feCheck(createFiniteElementSpace(U_cylindre, &mesh, dim = 1, LINE, order, "U", "cylinder",
                                     degreeQuadrature, funZero));
    feCheck(createFiniteElementSpace(U_surface, &mesh, dim = 2, TRI, order, "U", "surface",
                                     degreeQuadrature, funZero));
    feCheck(createFiniteElementSpace(V_walls, &mesh, dim = 1, LINE, order, "V", "walls",
                                     degreeQuadrature, funZero));
    feCheck(createFiniteElementSpace(V_inlet, &mesh, dim = 1, LINE, order, "v", "inlet",
                                     degreeQuadrature, funZero));
    feCheck(createFiniteElementSpace(V_cylindre, &mesh, dim = 1, LINE, order, "v", "cylinder",
                                     degreeQuadrature, funZero));
    feCheck(createFiniteElementSpace(V_surface, &mesh, dim = 2, TRI, order, "v", "surface",
                                     degreeQuadrature, funZero));
    feCheck(createFiniteElementSpace(P_surface, &mesh, dim = 2, TRI, order - 1, "P", "surface",
                                     degreeQuadrature, funSolP));
    feCheck(createFiniteElementSpace(Lu_cylinder, &mesh, dim = 1, LINE, order, "Lu", "cylinder",
                                     degreeQuadrature, funZero));
    feCheck(createFiniteElementSpace(Lv_cylinder, &mesh, dim = 1, LINE, order, "Lv", "cylinder",
                                     degreeQuadrature, funZero));
    feCheck(createFiniteElementSpace(Bu_cylinder, &mesh, dim = 1, LINE, order, "Bu", "cylinder",
                                     degreeQuadrature, funZero));
    feCheck(createFiniteElementSpace(Bv_cylinder, &mesh, dim = 1, LINE, order, "Bv", "cylinder",
                                     degreeQuadrature, funZero));
    // Define the set of all finite elements spaces and the set of feSpaces
    // forming the essential boundary conditions. The second set must always be
    // a subset of the first.
    std::vector<feSpace *> spaces = {U_walls,     U_inlet,     U_cylindre, U_surface, V_walls,
                                     V_inlet,     V_cylindre,  V_surface,  P_surface, Lu_cylinder,
                                     Lv_cylinder, Bu_cylinder, Bv_cylinder};
    std::vector<feSpace *> essentialSpaces = {U_inlet, V_walls, V_inlet};

    // Create a "meta" structure containing the numbering of each unknown field.
    feMetaNumber metaNumber(&mesh, spaces, essentialSpaces);

    // Create the solution structure (contains essentially the solution at all DOFs).
    feSolution sol(&mesh, spaces, essentialSpaces, &metaNumber);

    // Define (bi-)linear forms : here we define a bilinear symmetric diffusion form which
    // is the integral of k*(grad phi)_i*(grad phi)_j and a linear form which is the integral
    // of f*phi.
    std::vector<feSpace *> spacesNS2D = {U_surface, V_surface, P_surface};
    std::vector<feSpace *> spaceWeak_U = {U_cylindre, Lu_cylinder, Bu_cylinder};
    std::vector<feSpace *> spaceWeak_V = {V_cylindre, Lv_cylinder, Bv_cylinder};

    feBilinearForm NS2D(spacesNS2D, &mesh, degreeQuadrature,
                        new feSysElm_2D_NavierStokes(stokesParam, nullptr));
    feBilinearForm weakBC_UCylinder(spaceWeak_U, &mesh, degreeQuadrature,
                                    new feSysElm_1D_weakBC_edo1(1.0, funZero));
    feBilinearForm weakBC_VCylinder(spaceWeak_V, &mesh, degreeQuadrature,
                                    new feSysElm_1D_weakBC_edo1(1.0, funZero));

    // Initialize the linear system. Assembly of the elementary matrices and RHS is
    // performed in the solve step. Two linear solvers are available :
    // MKL Pardiso (direct solver) and PETSc (collection of iterative solvers).
    feLinearSystem *system;
    feCheck(createLinearSystem(system, MKLPARDISO, {&NS2D, &weakBC_UCylinder, &weakBC_VCylinder},
                               &metaNumber, &mesh, argc, argv));
    // createLinearSystem(system, PETSC, {&NS2D, &weakBC_UCylinder, &weakBC_VCylinder}, &metaNumber,
    // &mesh, argc, argv));

    feComputer *Compute_intLu = new feComputer(Lu_cylinder, &mesh, &metaNumber, "L2Norm_1Field");
    feComputer *Compute_intLv = new feComputer(Lv_cylinder, &mesh, &metaNumber, "L2Norm_1Field");
    std::vector<feComputer *> comput = {Compute_intLu,  Compute_intLv};

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
    double t1 = 60.;
    int nTimeSteps = 60 * pow(2, iter);
    TT[iter] = nTimeSteps;
    feTolerances tol{1e-5, 1e-5, 10};
    tic();
    std::string CodeIni = "BDF1/DCF";
    feCheck(createTimeIntegrator(solver, BDF2, tol, system, &metaNumber, &sol, &mesh, comput,
                                 exportData, t0, t1, nTimeSteps, CodeIni));
    feCheck(solver->makeSteps(nTimeSteps));
    toc();
    // std::vector<double> &normLu = solver->getNorm(0);
    std::vector<double> &normLu = Compute_intLu->getResult();
    std::vector<double> &normLv = Compute_intLv->getResult();
    // std::cout << "norm Lv" << normLv[2] << std::endl;
    // std::cout << "norm Lv" << normLv[2] << std::endl;

    std::cout << "Int LU" << std::endl;
    for(int i = 0; i < normLu.size(); i++) { std::cout << normLu[i] << std::endl; }

    std::cout << "Int LV" << std::endl;
    for(int i = 0; i < normLu.size(); i++) { std::cout << normLv[i] << std::endl; }

    // intL_u[2 * iter] = *std::max_element(normLu.begin(), normLu.end());
    // intL_v[2 * iter] = *std::max_element(normLv.begin(), normLv.end());
    // Post-process the solution

    // Free the used memory
    delete solver;
    delete exporter;
    delete U_walls;
    delete U_inlet;
    delete U_cylindre;
    delete U_surface;
    delete V_walls;
    delete V_inlet;
    delete V_cylindre;
    delete V_surface;
    delete P_surface;
    delete Lu_cylinder;
    delete Lv_cylinder;
    delete Bu_cylinder;
    delete Bv_cylinder;
    delete funZero;
    delete funSolU_inlet;
    delete funSolU_inlet_dot;
    delete funSolP;

    printf("\n");
    printf("\n");
    printf("%12s \t %12s \t %12s \t %12s \n", "nSteps", "nElm", "||L||_BDF2", "||L||_DC3");
    for(int i = 0; i < nIter; ++i)
      printf("%12d \t %12d \t %12.6e \t %12.6e \n", TT[i], nElm[i], intL_u[2 * i], intL_v[2 * i]);
    printf("\n");
    printf("\n");

    // petscFinalize();
  }
  return 0;
}