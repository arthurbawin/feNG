#include "feAPI.h"

double fSol(const double t, const std::vector<double> x, const std::vector<double> par)  //solution exacte  
{
  return pow(x[0], 6);
}

double fSource(const double t, const std::vector<double> pos, const std::vector<double> par)  //terme source  //par: parametre
{
  double k = par[0];
  return k * 30. * (pow(pos[0], 4));
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  return 0.0;
}

int main(int argc, char **argv) {
  petscInitialize(argc, argv);

  // Set the default parameters.
  const char *meshFile = "square2.msh";
  int verbosity = 2;
  int order = 2;                   //par defaut 1 sinon lors de la compilation -o --order XX  ? 
  int degreeQuadrature = 10;       //par defaut 10 sinon lors de la compilation -dquad --degreeQuadrature XX  ?

  feOptionsParser options(argc, argv);
  options.addOption(&meshFile, "-m", "--mesh", "Mesh file");
  options.addOption(&order, "-o", "--order", "Finite element space order");
  options.addOption(&degreeQuadrature, "-dquad", "--degreeQuadrature", "Degree of the quadrature");
  options.addOption(&verbosity, "-v", "--verbosity", "Verbosity level");
  feCheck(options.parse());

  setVerbose(verbosity);
  feMesh2DP1 mesh(meshFile);

  // std::vector <int> list=mesh.getList();  
  // for (int i=0;i<list.size();++i){
  //   feInfo("%d - %d",i,list[i]);
  // }
  // feInfo("%d",mesh.getNbInteriorElems());

  // std::vector <int> nbElmPerColor=mesh.getNbElmPerColor();  
  // std::vector<int> index=mesh.getIndexStartColorInList();

  // for (int i=0;i<nbElmPerColor.size();++i){
  //   feInfo("%d = %d - %d",i,nbElmPerColor[i],index[i]);
  // }
  // feInfo("%d",mesh.getNbColor());

  // int nbColor=mesh.getNbColor();
  // std::vector <int> colorElm=mesh.getColorElm();
  // for (int i=0;i<colorElm.size();++i){
  //   feInfo("%d = %d",i,colorElm[i]);
  // }
  // feInfo("%d",mesh.getNbInteriorElems());
  
  // int nbColor=mesh.getNbColor();
  // feInfo("%d",nbColor);

  //  #pragma omp parallel
  
  // {
  //   int ID = omp_get_thread_num();
  //   printf("hello(%d/%d)\n", ID, omp_get_num_threads());
  //   printf("world(%d/%d)\n", ID, omp_get_num_threads());

    
  // }

  // std::vector<std::vector<int>> myVec = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9,10,11,12},{13,14,15,16}};

  // // int cnt = 0;
  // #pragma omp parallel
  // {
  //   int cnt = 0;
  //   #pragma omp for schedule(static)
  //   for(int i = 0; i < 100; ++i){
  //   // printf("Printing %3d from thread %d\n", i, omp_get_thread_num());
  //   printf("Outer loop : Thread %d has printed %d times\n", omp_get_thread_num(), cnt++);

  //   for(int j = 0; j < 4; ++j){
  //     printf("Inner loop : myVec[%d][%d] = %d\n", omp_get_thread_num(), j, myVec[omp_get_thread_num()][j]);
  //   }

  // }
  // }
  

  
  double k = 1.0;
  feFunction *funSol = new feFunction(fSol, {});
  feFunction *funSource = new feFunction(fSource, {k});
  feFunction *funZero = new feFunction(fZero, {});

  int dim;
  feSpace *uBord, *uDomaine, *vBord, *vDomaine;
  feCheck(createFiniteElementSpace(uBord, &mesh, dim = 1, LINE, order, "U", "Bord",
                                   degreeQuadrature, funSol));
  feCheck(createFiniteElementSpace(uDomaine, &mesh, dim = 2, TRI, order, "U", "Domaine",
                                   degreeQuadrature, funZero));

  std::vector<feSpace *> spaces = {uBord, uDomaine};
  std::vector<feSpace *> essentialSpaces = {uBord};

  feMetaNumber metaNumber(&mesh, spaces, essentialSpaces);
  feSolution sol(&mesh, spaces, essentialSpaces, &metaNumber);

  feBilinearForm sourceU({uDomaine}, &mesh, degreeQuadrature, new feSysElm_2D_Source(1.0, funSource));        
  feBilinearForm diffU({uDomaine}, &mesh, degreeQuadrature, new feSysElm_2D_Diffusion(k, nullptr));

  feLinearSystem *system;
  
  
  //feCheck(
  //  createLinearSystem(system, MKLPARDISO, {&diffU, &sourceU}, &metaNumber, &mesh, argc, argv));
  feCheck(createLinearSystem(system, PETSC, spaces, {&diffU, &sourceU}, &metaNumber, &mesh, argc, argv));
  

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

  // feInfo("test");
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

  //sol.printSol("solution_2_thread");

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

  
  //std::vector<feBilinearForm *> formResiduals=system->getformResiduals();

  // feInfo("Size formMatrices :%d" ,formResiduals.size());

  // #pragma omp parallel
  // {
  //   feBilinearForm diffU({uDomaine}, &mesh, degreeQuadrature, new feSysElm_2D_Diffusion(k, nullptr));
  //   #pragma omp for
  //   for(int i = 0; i < mesh.getNbInteriorElems(); ++i){
  //     diffU.initialize(&metaNumber, &mesh, &sol, i);
  //     feInfo("From thread %d/%d : Lecture elm %d : [%3d %3d %3d %3d %3d %3d]", 
  //       omp_get_thread_num(),
  //       omp_get_num_threads(), i,
  //       diffU._adr[0],
  //       diffU._adr[1],
  //       diffU._adr[2],
  //       diffU._adr[3],
  //       diffU._adr[4],
  //       diffU._adr[5]);
  //   }
  // }

 

  return 0;
}