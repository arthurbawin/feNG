/*

  Convergence study for the heat equation
  - div(k grad u) + f = 0

  Manufactured solution u(x,y) = sin(pi*x)*sin(pi*y)

*/

#include "feAPI.h"


double uSol_f(const double /* t */, const std::vector<double> &pos, const std::vector<double> &/*par*/)
{
  double x = pos[0];
  double y = pos[1];
  return sin(M_PI*x)*sin(M_PI*y);
}

double uSource_f(const double /* t */, const std::vector<double> &pos, const std::vector<double> &/*par*/)
{
  double x = pos[0];
  double y = pos[1];
  return -2. * M_PI*M_PI * sin(M_PI*x)*sin(M_PI*y);
}

int main(int argc, char **argv)
{
  initialize(argc, argv);
  setVerbose(1);

  int dim = 2;
  int degreeQuadrature = 10;
  const char *solverType = "PETSc_MUMPS";
  int ownershipSplit = -1;
  double k = 1.;

  // Number of meshes for the convergence study
  int nMesh = 4;
  int order = 1;

  feOptionsParser options(argc, argv);
  options.addOption(&dim, "-d", "--dim", "Space dimension (2D or 3D)");
  options.addOption(&nMesh, "-n", "--nMesh", "Number of meshes for the convergence study");
  options.addOption(&order, "-o", "--order", "Finite element space polynomial order");
  feCheck(options.parse());

  if(dim != 2 && dim != 3) {
    feErrorMsg(FE_STATUS_ERROR, "Unsupported space dimension : %d", dim);
    return 1;
  }

  feFunction uSol(uSol_f, {k});
  feFunction uSource(uSource_f, {k});
  feConstantFunction zero(0);
  feConstantFunction one(1.);
  feConstantFunction diffusivity(k);

  // Number of topmost dimension elements in mesh and error
  std::vector<int> nElm(nMesh, 0);
  std::vector<double> L2errorU(nMesh, 0.0);
  std::vector<double> convRate(nMesh, 0.0);

  // Compute the solution on each mesh
  // By default, the meshes are square*.msh in 2D and cube*.msh in 3D
  // The mesh cube1.msh has no vertex inside the domain (only faces),
  // so start convergence study with cube2.msh
  for(int i = 0; i < nMesh; ++i)
  {
    std::string meshFile;
    if(dim == 2) meshFile = "../data/square" + std::to_string(i + 1) + ".msh";
    if(dim == 3) meshFile = "../data/cube"   + std::to_string(i + 2) + ".msh";

    //
    // Solve the heat equation on this mesh
    //
    feMesh2DP1 mesh(meshFile);
    // Get the number of topmost dimension (non-boundary) elements
    nElm[i] = mesh.getNumInteriorElements();

    feSpace *uBord, *uDomaine;
    feCheck(createFiniteElementSpace(   uBord, &mesh, elementType::LAGRANGE, order, "U",    "Bord", degreeQuadrature, &uSol));
    feCheck(createFiniteElementSpace(uDomaine, &mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &zero));

    std::vector<feSpace *> spaces = {uDomaine, uBord};
    std::vector<feSpace *> essentialSpaces = {uBord};

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);

    feBilinearForm *diff = nullptr, *source = nullptr;
    if(dim == 2)
      feCheck(createBilinearForm(diff, {uDomaine}, new feSysElm_Diffusion<2>(&diffusivity)));
    if(dim == 3)
      feCheck(createBilinearForm(diff, {uDomaine}, new feSysElm_Diffusion<3>(&diffusivity)));
    feCheck(createBilinearForm(source, {uDomaine}, new feSysElm_Source(&uSource)));
    std::vector<feBilinearForm*> forms = {diff, source};

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

    feNorm *errorU;
    feCheck(createNorm(errorU, L2_ERROR, {uDomaine}, &sol, &uSol));
    std::vector<feNorm*> norms = {errorU};

    feExporter *exporter;
    feCheck(createVisualizationExporter(exporter, VTK, &numbering, &sol, &mesh, spaces));
    std::string vtkFileRoot = "example2_mesh" + std::to_string(i + 1);
    feExportData exportData = {exporter, 1, vtkFileRoot};

    feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 50, 4, 1e-1};

    TimeIntegrator *solver;
    feCheck(createTimeIntegrator(solver, timeIntegratorScheme::STATIONARY, NLoptions, system, &sol, &mesh, norms, exportData));
    feCheck(solver->makeSteps(1));

    std::vector<std::vector<double>> postProc = solver->getPostProcessingData();
    L2errorU[i] = postProc[0][0];

    delete errorU;
    delete solver;
    delete exporter;
    delete system;
    for(feBilinearForm* f : forms)
      delete f;
    for(feSpace *s : spaces)
      delete s;
  }

  // Compute the convergence rate assuming isotropic elements
  for(int i = 1; i < nMesh; ++i) {
    convRate[i] = -log(L2errorU[i] / L2errorU[i-1]) / log( pow( nElm[i] / nElm[i-1], 1./(double) dim ) );
  }
  feInfo("%12s \t %12s \t %12s", "nElm", "||E_U||", "rate");
  for(int i = 0; i < nMesh; ++i) {
    feInfo("%12d \t %12.6e \t %12.6e", nElm[i], L2errorU[i], convRate[i]);
  }

  finalize();
  return 0;
}