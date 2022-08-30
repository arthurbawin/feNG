/*

  Solving the steady-state heat equation on a triangle mesh.

*/

#include "feAPI.h"

double fSol(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  // return pow(x[0], 6);
  double x = pos[0];
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

double fNeumannBC(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  double x = pos[0];
  double u0 = 1.0;
  double L = 1.0;
  double H = 1.0;
  return u0 / H * x / L * (1 - x / L);
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  return 0.0;
}

int main(int argc, char **argv)
{
  // tic();

  petscInitialize(argc, argv);

  // Set the default parameters.
  const char *meshFile = "../data/mesh/finer_27k.msh";
  int verbosity = 2;
  int order = 2;
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

  // Create a mesh structure from a Gmsh mesh file (version 2.2 or 4.1+)
  feMesh2DP1 mesh(meshFile, true, true);

  exit(-1);
  // feMesh2DP1 mesh(meshFile, true, true);
  feInfo("There are %d elms in the mesh\n", mesh.getNbInteriorElems());

  // Create function structures for the analytic solution and the source term
  // Here the model PDE is the stationary heat equation div(k grad u) + f = 0,
  // with k the thermal conductivity.
  double k = 1.0;
  feFunction *funSol = new feFunction(fSol, {});
  feFunction *funSource = new feFunction(fSource, {k});
  feFunction *funZero = new feFunction(fZero, {});

  // Define a finite element space on each subdomain of the computational domain.
  // Here the mesh contains only a 2D subdomain ("Domaine") and a 1D boundary ("Bord"),
  // and the solution field is named "U". Quadrature rules are defined on the finite
  // element space, with parameter "degreeQuadrature" governing the number of quad nodes.
  // The function provided is used to initialize the degrees of freedom on the feSpace.
  int dim;
  // feSpace *uDirichlet, *uNeumannHaut, *uNeumannBas, *uDomaine;
  // feCheck(createFiniteElementSpace(uDirichlet, &mesh, dim = 1, LINE, order, "U", "Dirichlet",
  //                                  degreeQuadrature, funZero));
  // feCheck(createFiniteElementSpace(uNeumannHaut, &mesh, dim = 1, LINE, order, "U", "NeumannHaut",
  //                                  degreeQuadrature, funZero));
  // feCheck(createFiniteElementSpace(uNeumannBas, &mesh, dim = 1, LINE, order, "U", "NeumannBas",
  //                                  degreeQuadrature, funZero));
  // feCheck(createFiniteElementSpace(uDomaine, &mesh, dim = 2, TRI, order, "U", "Domaine",
  //                                  degreeQuadrature, funZero));

  feSpace *uBord, *uDomaine;
  feCheck(createFiniteElementSpace(uBord, &mesh, dim = 1, LINE, order, "U", "Bord",
                                   degreeQuadrature, funSol));
  feCheck(createFiniteElementSpace(uDomaine, &mesh, dim = 2, TRI, order, "U", "Domaine",
                                   degreeQuadrature, funZero));

  // Define the set of all finite elements spaces and the set of feSpaces
  // forming the essential boundary conditions. The second set must always be
  // a subset of the first.
  // std::vector<feSpace *> spaces = {uDirichlet, uNeumannHaut, uNeumannBas, uDomaine};
  std::vector<feSpace *> spaces = {uBord, uDomaine};
  // std::vector<feSpace *> essentialSpaces = {uDirichlet};

  std::vector<feSpace *> essentialSpaces = {uBord};

  // Create a "meta" structure containing the numbering of each unknown field.
  feMetaNumber metaNumber(&mesh, spaces, essentialSpaces);

  // Create the solution structure (contains essentially the solution at all DOFs).
  feSolution sol(&mesh, spaces, essentialSpaces, &metaNumber);

  // Define (bi-)linear forms : here we define a bilinear symmetric diffusion form which
  // is the integral of k*(grad phi)_i*(grad phi)_j and a linear form which is the integral
  // of f*phi_i.
  feBilinearForm diffU({uDomaine}, &mesh, degreeQuadrature, new feSysElm_2D_Diffusion(k, nullptr));
  feBilinearForm sourceU({uDomaine}, &mesh, degreeQuadrature,
                         new feSysElm_2D_Source(1.0, funSource));
  // feBilinearForm neumannHaut({uNeumannHaut}, &mesh, degreeQuadrature,
  //                            new feSysElm_1D_NeumannBC(1.0, funNeumann));
  // feBilinearForm neumannBas({uNeumannBas}, &mesh, degreeQuadrature,
  //                           new feSysElm_1D_NeumannBC(1.0, funZero));


  // Initialize the linear system. Assembly of the elementary matrices and RHS is
  // performed in the solve step. Two linear solvers are available :
  // MKL Pardiso (direct solver) and PETSc (collection of iterative solvers).
  feLinearSystem *system;
  // feCheck(createLinearSystem(system, MKLPARDISO, {&diffU, &neumannHaut, &neumannBas},
  // &metaNumber, &mesh, argc, argv));
  feCheck(createLinearSystem(system, MKLPARDISO, spaces, {&diffU, &sourceU}, &metaNumber,
                             &mesh, argc, argv));
  // feCheck(createLinearSystem(system, MKLPARDISO, {&diffU, &sourceU}, &metaNumber, &mesh, argc,
  // argv));

  // Define post-processing tools to compute norms and whatnot (norms will be replaced by
  // fePostProc)
  feComputer normU(uDomaine, &mesh, &metaNumber, "L2Norm_1Field", funSol);
  std::vector<feComputer *> norms = {&normU};
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
  feTolerances tol{1e-9, 1e-8, 3, 0, 0};
  feCheck(createTimeIntegrator(solver, STATIONARY, tol, system, &metaNumber, &sol, &mesh, norms,
                               exportData));
  feCheck(solver->makeSteps(0));

  // Print the solution to check if it is the interpolant
  std::vector<Vertex> &vertices = mesh.getVertices();
  feNumber *n = metaNumber.getNumbering(0);
  for(int iVertex = 0; iVertex < vertices.size(); ++iVertex){
    int iDOF = n->getDOFNumberAtVertex(iVertex);

    Vertex *v = mesh.getVertex(iVertex);
    std::vector<double> pos{v->x(),v->y(),0.};
    std::vector<double> par(1,0.);

    if(iDOF >= 0) {
      // There is a degree of freedom at this mesh vertex
      // output << sol[iDOF] << std::endl;
      feInfo("DOF FOUND at (%f - %f) : solEF = %+1.16e - sol = %+-1.16e - e = %+-1.16e", v->x(), v->y(), 
        sol.getSolAtDOF(iDOF), fSol(0.,pos,par), fabs(sol.getSolAtDOF(iDOF) - fSol(0.,pos,par)));
    } else {
      feInfo("ERREUR at (%f - %f)", v->x(), v->y());
    }
    // fprintf(f, "%+-16.16e\n", sol[iDOF]);
  }


  // Post-process the solution
  fePostProc post(uDomaine, &mesh, &metaNumber, funSol);
  feInfo("Measure = %f", post.computeMeasure());
  feInfo("Integral of sol = %10.10f", post.computeSolutionIntegral(&sol));
  feInfo("Integral of fun = %10.10f", post.computeFunctionIntegral(funSol, 0.0));
  feInfo("L2 Error = %10.10f", post.computeL2ErrorNorm(&sol));

  // Free the used memory
  delete solver;
  delete exporter;
  // delete uDirichlet;
  // delete uNeumannHaut;
  // delete uNeumannBas;
  delete uDomaine;
  delete uBord;
  delete funSol;
  delete funSource;

  petscFinalize();
  // toc();
  return 0;
}