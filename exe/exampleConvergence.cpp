/*

  Convergence study for the heat equation.

*/

#include "feAPI.h"

double fSol(const double t, const std::vector<double> pos, const std::vector<double> par)
{
  double x = pos[0];
  // double y = pos[1];
  return pow(x, 6);
}

double fSource(const double t, const std::vector<double> x, const std::vector<double> par)
{
  double k = par[0];
  return k * 30. * pow(x[0], 4);
}


double fZero(const double t, const std::vector<double> x, const std::vector<double> par)
{
  return 0.0;
}


int main(int argc, char **argv)
{
  petscInitialize(argc, argv);
  int order = 1;
  int nMesh = 1;

  feOptionsParser options(argc, argv);
  options.addOption(&order, "-o", "--order", "Finite element space order");
  options.addOption(&nMesh, "-n", "--nMeshes", "Number of meshes for the convergence study");
  feCheck(options.parse());

  setVerbose(2);

  double k = 1.0;
  feFunction *funSol = new feFunction(fSol, {});
  feFunction *funSource = new feFunction(fSource, {k});
  feFunction *funZero = new feFunction(fZero, {});

  std::vector<int> nElm(nMesh, 0);
  std::vector<double> L2errorU(2 * nMesh, 0.0);

  for(int i = 0; i < nMesh; ++i){
    std::string meshFile = "../data/Convergence/m" + std::to_string(i+1) + ".msh";

    feMesh2DP1 mesh(meshFile);
    nElm[i] = mesh.getNbInteriorElems();

    // feSpace *uBord, *uDomaine;
    int degreeQuadrature = 30;
    int dim;
    feSpace *uBord, *uDomaine;
    feCheck(createFiniteElementSpace(uBord, &mesh, 1, LINE, order, "U", "Bord", degreeQuadrature, funSol));
    feCheck(createFiniteElementSpace(uDomaine, &mesh, 2, TRI, order, "U", "Domaine", degreeQuadrature, funZero));

    std::vector<feSpace *> spaces = {uBord, uDomaine};
    std::vector<feSpace *> essentialSpaces = {uBord};

    feMetaNumber metaNumber(&mesh, spaces, essentialSpaces);
    feSolution sol(&mesh, spaces, essentialSpaces, &metaNumber);

    feBilinearForm diffU({uDomaine}, &mesh, degreeQuadrature, new feSysElm_2D_Diffusion(k, nullptr));
    feBilinearForm sourceU({uDomaine}, &mesh, degreeQuadrature, new feSysElm_2D_Source(1.0, funSource));

    feLinearSystem *system;
    feCheck(createLinearSystem(system, MKLPARDISO, spaces, {&diffU, &sourceU}, &metaNumber, &mesh, argc, argv));

    // Norms will be replaced by postproc
    feNorm normU(uDomaine, &mesh, degreeQuadrature, funSol);
    std::vector<feNorm *> norms = {&normU};

    feExporter *exporter;
    feCheck(createVisualizationExporter(exporter, VTK, &metaNumber, &sol, &mesh, spaces));
    int exportEveryNSteps = 1;
    std::string vtkFileRoot = "out" + std::to_string(i);
    feExportData exportData = {exporter, exportEveryNSteps, vtkFileRoot};

    TimeIntegrator *solver;
    feTolerances tol{1e-9, 1e-8, 3};
    feCheck(createTimeIntegrator(solver, STATIONARY, tol, system, &metaNumber, &sol, &mesh, norms,
                                 exportData));
    feCheck(solver->makeSteps(0));

    fePostProc post(uDomaine, &mesh, &metaNumber, funSol);
    L2errorU[2*i] = post.computeL2ErrorNorm(&sol);

    delete solver;
    delete exporter;
    delete uBord;
    delete uDomaine;
  }

  // Compute the convergence rate
  for(int i = 1; i < nMesh; ++i) {
    L2errorU[2 * i + 1] = -log(L2errorU[2 * i] / L2errorU[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
  }
  printf("%12s \t %12s \t %12s \n", "nElm", "||E_U||", "rate");
  for(int i = 0; i < nMesh; ++i){
    printf("%12d \t %12.6e \t %12.6e\n", nElm[i], L2errorU[2 * i], L2errorU[2 * i + 1]);
  }

  delete funSol;
  delete funSource;

  petscFinalize();

  return 0;
}