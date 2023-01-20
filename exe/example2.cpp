/*

  Convergence study for the heat equation.

*/

#include "feAPI.h"

double fSol(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  double x = pos[0];
  // double y = pos[1];
  return pow(x, 6);

  // // Solution for the Neumann problem
  // double u0 = 1.0;
  // double L = 2.0;
  // double H = 2.0;
  // double u = 0.0;
  // for(int n = 1; n < 100; n += 2){
  //   u += 1./pow(n*M_PI,4) * sin(n*M_PI*x/L) * cosh(n*M_PI*y/L) / sinh(n*M_PI*H/L); 
  // }
  // u*= 8.*u0*L/H;
  // return u;
}

double fSource(const double t, const std::vector<double> &x, const std::vector<double> &par)
{
  double k = par[0];
  return k * 30. * pow(x[0], 4);
}

// double fNeumannBC(const double t, const std::vector<double> pos, const std::vector<double> par)
// {
//   double x = pos[0];
//   double u0 = 1.0;
//   double L = 2.0;
//   double H = 2.0;
//   return u0/H*(1+x/L)*(1-x/L);
// }

// double fNeumannYBC(const double t, const std::vector<double> pos, const std::vector<double> par)
// {
//   double y = pos[1];
//   double u0 = 1.0;
//   double L = 2.0;
//   double H = 2.0;
//   return u0/H*(1+y/L)*(1-y/L);
// }

double fZero(const double t, const std::vector<double> &x, const std::vector<double> &par)
{
  return 0.0;
}

double fConstant(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  return par[0];
}

// double fOne(const double t, const std::vector<double> x, const std::vector<double> par)
// {
//   double c = par[0];
//   return c;
// }

int main(int argc, char **argv)
{
  // ProfilerStart();
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
  feFunction *kDiffusivity = new feFunction(fConstant, {1.0});
  // feFunction *funNeumann = new feFunction(fNeumannBC, {});
  // feFunction *funNeumannY = new feFunction(fNeumannYBC, {});
  // feFunction *funOne = new feFunction(fOne, {1.0});
  // feFunction *funTwo = new feFunction(fOne, {2.0});
  // feFunction *funThree = new feFunction(fOne, {3.0});

  std::vector<int> nElm(nMesh, 0);
  std::vector<double> L2errorU(2 * nMesh, 0.0);

  bool hasAnalyticSolution = true;

  // Compute the solution on each mesh
  for(int i = 0; i < nMesh; ++i){
    // std::string meshFile = "../data/Convergence/m" + std::to_string(i+3) + ".msh";
    std::string meshFile = "/home/bawina/Code/feNG/data/square" + std::to_string(i+1) + ".msh";
    // std::string meshFile = "wideCrack.msh";
    // std::string meshFile = "transfinite" + std::to_string(i+1) + "_P2.msh";

    feMesh2DP1 mesh(meshFile);
    nElm[i] = mesh.getNbInteriorElems();

    feSpace *uBord, *uDomaine;
    int degreeQuadrature = 12;
    int dim;
    feCheck(createFiniteElementSpace(uBord, &mesh, dim = 1, LINE, LAGRANGE, order, "U", "Bord", degreeQuadrature, funSol));
    feCheck(createFiniteElementSpace(uDomaine, &mesh, dim = 2, TRI, LAGRANGE, order, "U", "Domaine", degreeQuadrature, funZero));
    // feSpace *uSaut, *uHaut, *uBas, *uDroite, *uGauche, *uDomaine;
    // feCheck(createFiniteElementSpace(uHaut,    &mesh, dim = 1, LINE, order, "U", "Haut",    degreeQuadrature, funZero));
    // feCheck(createFiniteElementSpace(uBas,     &mesh, dim = 1, LINE, order, "U", "Bas",     degreeQuadrature, funZero));
    // feCheck(createFiniteElementSpace(uDroite,  &mesh, dim = 1, LINE, order, "U", "Droite",  degreeQuadrature, funZero));
    // feCheck(createFiniteElementSpace(uGauche,  &mesh, dim = 1, LINE, order, "U", "Gauche",  degreeQuadrature, funZero));
    // feCheck(createFiniteElementSpace(uSaut,    &mesh, dim = 1, LINE, order, "U", "Saut",    degreeQuadrature, funZero));
    // feCheck(createFiniteElementSpace(uDomaine, &mesh, dim = 2, TRI,  order, "U", "Domaine", degreeQuadrature, funZero));

    // std::vector<feSpace *> spaces = {uSaut, uHaut, uBas, uDroite, uGauche, uDomaine};
    // std::vector<feSpace *> spaces = {uHaut, uBas, uDroite, uGauche, uDomaine};
    std::vector<feSpace *> spaces = {uBord, uDomaine};
    // std::vector<feSpace *> essentialSpaces = {uDroite, uGauche};
    // std::vector<feSpace *> essentialSpaces = {uHaut, uBas, uGauche, uDroite};
    std::vector<feSpace *> essentialSpaces = {uBord};

    feMetaNumber metaNumber(&mesh, spaces, essentialSpaces);
    feSolution sol(metaNumber.getNbDOFs(), spaces, essentialSpaces);

    feBilinearForm diffU({uDomaine}, new feSysElm_2D_Diffusion(k, nullptr));
    feBilinearForm sourceU({uDomaine}, new feSysElm_2D_Source(1.0, funSource));

    feLinearSystem *system;
    feCheck(createLinearSystem(system, PETSC, {&diffU, &sourceU}, metaNumber.getNbUnknowns(), argc, argv));
    // feCheck(createLinearSystem(system, MKLPARDISO, spaces, {&diffU, &sourceU}, &metaNumber, &mesh, argc, argv));
    // feCheck(createLinearSystem(system, MKLPARDISO, spaces, {&diffU, &n1, &n2}, &metaNumber, &mesh, argc, argv));
    // feCheck(createLinearSystem(system, MKLPARDISO, spaces, {&diffU, &sautU}, &metaNumber, &mesh, argc, argv));
    // feCheck(createLinearSystem(system, MKLPARDISO, spaces, {&diffU, &n5}, &metaNumber, &mesh, argc, argv));

    // Norms will be replaced by postproc
    feNorm normU({uDomaine}, &mesh, degreeQuadrature, funSol);
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

    if(hasAnalyticSolution){
      fePostProc post(uDomaine, &mesh, &metaNumber, funSol);
      L2errorU[2*i] = post.computeL2ErrorNorm(&sol);
    } else{
      // Export solution for error estimation
      std::string solutionFile = "solution_" + std::to_string(i+1) + ".txt";
      sol.printSol(solutionFile);
    }

    delete solver;
    delete exporter;
    delete uBord;
    delete uDomaine;
  }

  // Error estimation if the analytic solution is not available
  if(!hasAnalyticSolution){
    for(int i = 0; i < nMesh-1; ++i){
      // Load current and next solution files
      std::string solFile1 = "solution_" + std::to_string(i+1) + ".txt";
      std::string solFile2 = "solution_" + std::to_string(i+2) + ".txt";
      feSolution *sol1 = new feSolution(solFile1);
      feSolution *sol2 = new feSolution(solFile2);

      // Create meshes, spaces and numberings for interpolation
      std::string mshFile1 = "m" + std::to_string(i+1) + ".msh";
      std::string mshFile2 = "m" + std::to_string(i+2) + ".msh";
      feMesh2DP1 mesh1(mshFile1);
      feMesh2DP1 mesh2(mshFile2);

      int dim, degreeQuadrature = 30;
      feSpace *uSaut1, *uHaut1, *uBas1, *uDroite1, *uGauche1, *uDomaine1;
      feSpace *uSaut2, *uHaut2, *uBas2, *uDroite2, *uGauche2, *uDomaine2;
      feCheck(createFiniteElementSpace(uHaut1,    &mesh1, dim = 1, LINE, LAGRANGE, order, "U", "Haut",    degreeQuadrature, funZero));
      feCheck(createFiniteElementSpace(uBas1,     &mesh1, dim = 1, LINE, LAGRANGE, order, "U", "Bas",     degreeQuadrature, funZero));
      feCheck(createFiniteElementSpace(uDroite1,  &mesh1, dim = 1, LINE, LAGRANGE, order, "U", "Droite",  degreeQuadrature, funZero));
      feCheck(createFiniteElementSpace(uGauche1,  &mesh1, dim = 1, LINE, LAGRANGE, order, "U", "Gauche",  degreeQuadrature, funZero));
      feCheck(createFiniteElementSpace(uDomaine1, &mesh1, dim = 2, TRI,  LAGRANGE, order, "U", "Domaine", degreeQuadrature, funZero));
      feCheck(createFiniteElementSpace(uSaut1,    &mesh1, dim = 1, LINE, LAGRANGE, order, "U", "Saut",    degreeQuadrature, funZero));
      feCheck(createFiniteElementSpace(uHaut2,    &mesh2, dim = 1, LINE, LAGRANGE, order, "U", "Haut",    degreeQuadrature, funZero));
      feCheck(createFiniteElementSpace(uBas2,     &mesh2, dim = 1, LINE, LAGRANGE, order, "U", "Bas",     degreeQuadrature, funZero));
      feCheck(createFiniteElementSpace(uDroite2,  &mesh2, dim = 1, LINE, LAGRANGE, order, "U", "Droite",  degreeQuadrature, funZero));
      feCheck(createFiniteElementSpace(uGauche2,  &mesh2, dim = 1, LINE, LAGRANGE, order, "U", "Gauche",  degreeQuadrature, funZero));
      feCheck(createFiniteElementSpace(uDomaine2, &mesh2, dim = 2, TRI,  LAGRANGE, order, "U", "Domaine", degreeQuadrature, funZero));
      feCheck(createFiniteElementSpace(uSaut2,    &mesh2, dim = 1, LINE, LAGRANGE, order, "U", "Saut",    degreeQuadrature, funZero));

      std::vector<feSpace *> spaces1 = {uSaut1, uHaut1, uBas1, uDroite1, uGauche1, uDomaine1};
      std::vector<feSpace *> spaces2 = {uSaut2, uHaut2, uBas2, uDroite2, uGauche2, uDomaine2};
      std::vector<feSpace *> essentialSpaces1 = {uHaut1, uBas1, uGauche1, uDroite1};
      std::vector<feSpace *> essentialSpaces2 = {uHaut2, uBas2, uGauche2, uDroite2};

      feMetaNumber metaNumber1(&mesh1, spaces1, essentialSpaces1);
      feMetaNumber metaNumber2(&mesh2, spaces2, essentialSpaces2);

      // fePostProc post(uDomaine1, &mesh1, &metaNumber1, funSol);
      feNorm norm({uDomaine1}, &mesh1, degreeQuadrature);
      norm.computeErrorNormFromExternalSolution(&metaNumber1, sol1, &mesh1, 
        &metaNumber2, sol2, &mesh2, spaces2);
      L2errorU[2*i] = norm.getNorm();

      delete sol1;
      delete sol2;
    }
  }

  for(auto val : L2errorU){
    feInfo("val = %f", val);
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

  // ProfilerStop();
  return 0;
}