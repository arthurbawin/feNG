/*

  Solving the steady-state heat equation on a triangle mesh.

*/
#include <iostream>
#include <fstream>

#include "feAPI.h"

double fSol(const double t, const std::vector<double> &x, const std::vector<double> &par)
{
  return pow(x[0],2) + pow(x[1],2) + pow(t,6);
}

double fSource(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  double k = par[0];
  return -6.*pow(t,5) + k*4.;
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
  
  int verbosity = 2;
  int order = 3;
  int degreeQuadrature = 10;
  int nIter = 1;

  double k = 0.0001;
  feFunction *funSol = new feFunction(fSol, {});
  feFunction *funSource = new feFunction(fSource, {k});
  feFunction *funZero = new feFunction(fZero, {});

  feOptionsParser options(argc, argv);
  options.addOption(&nIter, "-n", "--nMeshes", "Number of meshes for the convergence study");
  options.addOption(&order, "-o", "--order", "Finite element space order");
  options.addOption(&degreeQuadrature, "-dquad", "--degreeQuadrature", "Degree of the quadrature");
  options.addOption(&verbosity, "-v", "--verbosity", "Verbosity level");
  feCheck(options.parse());
  setVerbose(verbosity);

  std::vector<double> normBDF(2 * nIter, 0.0);
  std::vector<double> normDC2(2 * nIter, 0.0);
  std::vector<double> normDC3(2 * nIter, 0.0);
  std::vector<double> normDC4(2 * nIter, 0.0);
  std::vector<double> normDC5(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);
  std::vector<int> TT(nIter, 0);


  for(int iter = 0; iter < nIter; ++iter) {

    // std::string meshName ="../data/Convergence/finer" + std::to_string(iter+1) + ".msh";
    std::string meshName = "../data/Convergence/squareDirichlet1.msh";

    feMesh2DP1 mesh(meshName);
    nElm[iter] = mesh.getNbInteriorElems();


    int dim;
    feSpace *uBord, *uDomaine;
    feCheck(createFiniteElementSpace(uBord, &mesh, dim = 1, LINE, order, "U", "Bord", degreeQuadrature, funSol));
    feCheck(createFiniteElementSpace(uDomaine, &mesh, dim = 2, TRI, order, "U", "Surface", degreeQuadrature, funSol));

    std::vector<feSpace *> spaces = {uBord, uDomaine};
    std::vector<feSpace *> essentialSpaces = {uBord};

    feMetaNumber metaNumber(&mesh, spaces, essentialSpaces);

    feSolution sol(&mesh, spaces, essentialSpaces, &metaNumber);

    feBilinearForm masseU({uDomaine}, &mesh, degreeQuadrature, new feSysElm_2D_Masse(1.0, nullptr));
    feBilinearForm diffU({uDomaine}, &mesh, degreeQuadrature, new feSysElm_2D_Diffusion(k, nullptr));
    feBilinearForm sourceU({uDomaine}, &mesh, degreeQuadrature, new feSysElm_2D_Source(1.0, funSource));


    feLinearSystem *system;
    feCheck(createLinearSystem(system, PETSC, spaces, {&diffU, &sourceU, &masseU}, &metaNumber, &mesh, argc, argv));
    // feCheck(createLinearSystem(system, MKLPARDISO, spaces, {&diffU, &sourceU, &masseU}, &metaNumber, &mesh, argc, argv));

    feNormV2  *normU = new feNormV2({uDomaine}, &mesh, &metaNumber, "L2Norm_1Field", {funSol});
    std::vector<feNormV2 *> norms = {normU};

    feExporter *exporter;
    feCheck(createVisualizationExporter(exporter, VTK, &metaNumber, &sol, &mesh, spaces));
    int exportEveryNSteps = 1;
    std::string vtkFileRoot = "outDIFF_";
    feExportData exportData = {nullptr, exportEveryNSteps, vtkFileRoot};


    // Solver
    feSolverV2 *solver;
    feTolerancesV2 tol{1e-11, 1e-11, 50};
    // double criteriaNormU; double criteriaNormRes; double ratioNormU; double ratioNormRes; double nbIterBeforeRecompute; double convergedBeforeNIter;
    feMatrixComputeCriteria criteria{1e-0, 1e-0, 5, 5, 4, 5};
    feCheck(createSolver(solver, system, &metaNumber, Newton, tol, &criteria));


    //Resolution
    feTimeIntegrator *timeIntegrator;
    double t0 = 0.;
    double t1 = 1.;
    std::string codeDt = "Alt";
    int nTimeSteps = 10* pow(2, iter);
    TT[iter] = nTimeSteps;
    // feCheck(createTimeIntegrator(solver, DC5, tol, system, &metaNumber, &sol, &mesh, comput, exportData, t0, t1, nTimeSteps, CodeIni));
    feCheck(createTimeIntegratorV2(timeIntegrator, &mesh, solver, system, &metaNumber, &sol, DC5F_V2, "SOLDOT", t0, t1, nTimeSteps, norms, exportData, codeDt));
    feCheck(timeIntegrator->makeSteps());


    std::vector<double> normIterBDF = timeIntegrator->getINormBDF(0);
    std::vector<double> normIterDC2 = timeIntegrator->getINormDC2(0);
    std::vector<double> normIterDC3 = timeIntegrator->getINormDC3(0);
    std::vector<double> normIterDC4 = timeIntegrator->getINormDC4(0);
    std::vector<double> normIterDC5 = timeIntegrator->getINormDC5(0);
    normBDF[2*iter] = *std::max_element(normIterBDF.begin(), normIterBDF.end());
    normDC2[2*iter] = *std::max_element(normIterDC2.begin(), normIterDC2.end());
    normDC3[2*iter] = *std::max_element(normIterDC3.begin(), normIterDC3.end());
    normDC4[2*iter] = *std::max_element(normIterDC4.begin(), normIterDC4.end());
    normDC5[2*iter] = *std::max_element(normIterDC5.begin(), normIterDC5.end());

    // printf("\n");
    // feInfo("normBDF");
    // for (size_t i=0; i<normIterBDF.size(); i++){
    //   feInfo("%f", normIterBDF[i]);
    // }

    // printf("\n"); 
    // feInfo("normDC2");
    // for (size_t i=0; i<normIterDC2.size(); i++){
    //   feInfo("%f", normIterDC2[i]);
    // }

    // printf("\n"); 
    // feInfo("normDC3");
    // for (size_t i=0; i<normIterDC3.size(); i++){
    //   feInfo("%f", normIterDC3[i]);
    // }

    // printf("\n"); 
    // feInfo("normDC4");
    // for (size_t i=0; i<normIterDC4.size(); i++){
    //   feInfo("%f", normIterDC4[i]);
    // }

    // printf("\n"); 
    // feInfo("normDC5");
    // for (size_t i=0; i<normIterDC5.size(); i++){
    //   feInfo("%f", normIterDC5[i]);
    // }

    // printf("\n"); 
    // feInfo("normDC6");
    // for (size_t i=0; i<normIterDC6.size(); i++){
    //   feInfo("%f", normIterDC6[i]);
    // }



    delete solver;
    delete exporter;
    delete uDomaine;
    delete uBord;
  }


  // Calcul du taux de convergence
  for(int i = 1; i < nIter; ++i) {
    normBDF[2 * i + 1] = -log(normBDF[2 * i] / normBDF[2 * (i - 1)]) / log(2.);
    normDC2[2 * i + 1] = -log(normDC2[2 * i] / normDC2[2 * (i - 1)]) / log(2.);
    normDC3[2 * i + 1] = -log(normDC3[2 * i] / normDC3[2 * (i - 1)]) / log(2.);
    normDC4[2 * i + 1] = -log(normDC4[2 * i] / normDC4[2 * (i - 1)]) / log(2.);
    normDC5[2 * i + 1] = -log(normDC5[2 * i] / normDC5[2 * (i - 1)]) / log(2.);
    
  }

  printf("BDF");
  printf("%12s \t %12s \t %12s\n", "nElm", "||u-uh||", "tauxU");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12.6e \t %12.6e \n", nElm[i], normBDF[2 * i], normBDF[2 * i + 1]);


  printf("");
  printf("DC2");
  printf("%12s \t %12s \t %12s\n", "nElm", "||u-uh||", "tauxU");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12.6e \t %12.6e \n", nElm[i], normDC2[2 * i], normDC2[2 * i + 1]);


  printf("");
  printf("DC3");
  printf("%12s \t %12s \t %12s\n", "nElm", "||u-uh||", "tauxU");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12.6e \t %12.6e \n", nElm[i], normDC3[2 * i], normDC3[2 * i + 1]);


  printf("");
  printf("DC4");
  printf("%12s \t %12s \t %12s\n", "nElm", "||u-uh||", "tauxU");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12.6e \t %12.6e \n", nElm[i], normDC4[2 * i], normDC4[2 * i + 1]);


  printf("");
  printf("DC5");
  printf("%12s \t %12s \t %12s\n", "nElm", "||u-uh||", "tauxU");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12.6e \t %12.6e \n", nElm[i], normDC5[2 * i], normDC5[2 * i + 1]);
 
  delete funSol;
  delete funSource;
  delete funZero;

  petscFinalize();

  return 0;
}