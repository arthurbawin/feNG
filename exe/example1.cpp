/*

  Solving the steady-state heat equation on a triangle mesh.

*/

#include "feAPI.h"   //bibliotheque avec tous les references des fonctions

//dans cet exemple, on connait déjà la forme de la solution (u(x,y)=x⁶), on calcule alors le terme source avec l'équation de la diffusion(f=k*u''=30*k*x⁴)

double fSol(const double t, const std::vector<double> x, const std::vector<double> par)  //solution exacte  
{
  return pow(x[0], 6);
}

double fSource(const double t, const std::vector<double> pos, const std::vector<double> par)  //terme source  //par: parametre
{
  double k = par[0];
  return k * 30. * (pow(pos[0], 4) + pow(pos[1], 4));
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  return 0.0;
}

int main(int argc, char **argv)
{
  petscInitialize(argc, argv);

  // Set the default parameters.
  const char *meshFile = "square.msh";
  int verbosity = 2;
  int order = 1;                   //par defaut 1 sinon lors de la compilation -o --order XX  ? 
  int degreeQuadrature = 10;       //par defaut 10 sinon lors de la compilation -dquad --degreeQuadrature XX  ?

  // Create an option parser and parse the command line arguments.
  // If a command line argument is provided, it will overwrite the default parameter.
  // Each "addOption" adds an optional field, although it can be made required by
  // setting the 5th argument of addOption to "true".
  feOptionsParser options(argc, argv);                                                        //options pour la compilation ou l'execution? 
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
  feMesh2DP1 mesh(meshFile);                                                               //P1 car order=1 ? //on peut inclure d'autres pararametre (exemple bool curved), si pas inclut prend valeurs par defaut definies dans .h

  // std::vector<feCncGeo *> cncVec = mesh.getCncGeo();

  // feInfo("Mesh has %d connectivities", cncVec.size());

  // for(feCncGeo *cnc : cncVec){
  //   feInfo("Nom = %s", cnc->getID().c_str());
  //   feInfo("Interpolants geo = %s", cnc->getForme().c_str());
  //   feInfo("Nb noeuds = %d", cnc->getNbNodes());
  //   feInfo("Nb noeuds par element = %d", cnc->getNbNodePerElem());
  // }

  // feCncGeo *cnc = cncVec[0];
  // fePatch myPatch(cnc, &mesh);

  // std::set<int> elementConnectes = myPatch.getPatch(0);

  // for(auto e : elementConnectes)
  //   feInfo("%d", e);


  // return 0;

  // Create function structures for the analytic solution and the source term
  // Here the model PDE is the stationary heat equation div(k grad u) + f = 0,
  // with k the thermal conductivity.
  double k = 1.0;
  feFunction *funSol = new feFunction(fSol, {});                                           //solution exacte ? Comment est defini funsol alors que l'espace non def
  feFunction *funSource = new feFunction(fSource, {k});
  feFunction *funZero = new feFunction(fZero, {});

  // Define a finite element space on each subdomain of the computational domain.
  // Here the mesh contains only a 2D subdomain ("Domaine") and a 1D boundary ("Bord"),
  // and the solution field is named "U". Quadrature rules are defined on the finite
  // element space, with parameter "degreeQuadrature" governing the number of quad nodes.
  // The function provided is used to initialize the degrees of freedom on the feSpace.
  int dim;
  feSpace *uBord, *uDomaine, *vBord, *vDomaine;
  feCheck(createFiniteElementSpace(uBord, &mesh, dim = 1, LINE, order, "U", "Bord",
                                   degreeQuadrature, funSol));                                  //degreeQuadrature=10 => 10 noeuds de calcul ? 
  feCheck(createFiniteElementSpace(uDomaine, &mesh, dim = 2, TRI, order, "U", "Domaine",
                                   degreeQuadrature, funSol));                                  //où est defini le mesh ? car ligne 45 pas evident à comprendre, pourquoi pas defini comme mesh=feMesh2DP1(meshfile) 
                                                                                                //qu'est ce que creatFiniteElementSpace retourne ? la fonction semble modifier uBord et uDomain (space) ? 
  
  // Define the set of all finite elements spaces and the set of feSpaces
  // forming the essential boundary conditions. The second set must always be
  // a subset of the first.
  std::vector<feSpace *> spaces = {uBord, uDomaine};                                             //creation d'un vecteur avec toutes le noeuds (bord+domaine ??)
  std::vector<feSpace *> essentialSpaces = {uBord};                                              //correspond au condition limite (comment est definit essentialSpace avec c.l Neumann, ou dirichlet non imposé sur tous les bords)  

  // Create a "meta" structure containing the numbering of each unknown field.
  feMetaNumber metaNumber(&mesh, spaces, essentialSpaces);                                       //numérotation DDL ?? comment sont gerer les C.L ? 

  // Create the solution structure (contains essentially the solution at all DOFs).
  feSolution sol(&mesh, spaces, essentialSpaces, &metaNumber);

  // Define (bi-)linear forms : here we define a bilinear symmetric diffusion form which
  // is the integral of k*(grad phi)_i*(grad phi)_j and a linear form which is the integral
  // of f*phi.
  feBilinearForm diffU({uDomaine}, &mesh, degreeQuadrature, new feSysElm_2D_Diffusion(k, nullptr));  //diffU est une forme billinaire mais comment marche la fonction ? Qu'est-ce qu'elle retourne  exactement? 
  feBilinearForm sourceU({uDomaine}, &mesh, degreeQuadrature, new feSysElm_2D_Source(1.0, funSource));                

  // Initialize the linear system. Assembly of the elementary matrices and RHS is
  // performed in the solve step. Two linear solvers are available :
  // MKL Pardiso (direct solver) and PETSc (collection of iterative solvers).
  feLinearSystem *system;
  //feCheck(
  //  createLinearSystem(system, MKLPARDISO, {&diffU, &sourceU}, &metaNumber, &mesh, argc, argv));
  feCheck(createLinearSystem(system, PETSC, {&diffU, &sourceU}, &metaNumber, &mesh, argc, argv));

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
                               exportData));                                                           //BFD1, BFD2, etc ???
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

  petscFinalize();

  return 0;
}