/*

  Convergence study for the heat equation.

*/

#include "feAPI.h"

double fSol(const double t, const std::vector<double> pos, const std::vector<double> par)
{
  double x = pos[0];
  double y = pos[1];
  return pow(x, 6);

  // Solution for the Neumann problem
  double u0 = 1.0;
  double L = 2.0;
  double H = 2.0;
  double u = 0.0;
  for(int n = 1; n < 100; n += 2){
    u += 1./pow(n*M_PI,4) * sin(n*M_PI*x/L) * cosh(n*M_PI*y/L) / sinh(n*M_PI*H/L); 
  }
  u*= 8.*u0*L/H;
  return u;
}

double fSource(const double t, const std::vector<double> x, const std::vector<double> par)
{
  double k = par[0];
  return k * 30. * pow(x[0], 4);
}

double fNeumannBC(const double t, const std::vector<double> pos, const std::vector<double> par)
{
  double x = pos[0];
  double u0 = 1.0;
  double L = 2.0;
  double H = 2.0;
  return u0/H*(1+x/L)*(1-x/L);
}

double fNeumannYBC(const double t, const std::vector<double> pos, const std::vector<double> par)
{
  double y = pos[1];
  double u0 = 1.0;
  double L = 2.0;
  double H = 2.0;
  return u0/H*(1+y/L)*(1-y/L);
}

double fZero(const double t, const std::vector<double> x, const std::vector<double> par)
{
  return 0.0;
}

double fOne(const double t, const std::vector<double> x, const std::vector<double> par)
{
  double c = par[0];
  return c;
}

int main(int argc, char **argv)
{
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
  feFunction *funNeumann = new feFunction(fNeumannBC, {});
  feFunction *funNeumannY = new feFunction(fNeumannYBC, {});
  feFunction *funOne = new feFunction(fOne, {1.0});
  feFunction *funTwo = new feFunction(fOne, {2.0});
  feFunction *funThree = new feFunction(fOne, {3.0});

  std::vector<int> nElm(nMesh, 0);
  std::vector<double> L2errorU(2 * nMesh, 0.0);

  for(int i = 0; i < nMesh; ++i){
    std::string meshFile = "m" + std::to_string(i+1) + ".msh";
    // std::string meshFile = "wideCrack.msh";
    // std::string meshFile = "transfinite" + std::to_string(i+1) + "_P2.msh";

    feMesh2DP1 mesh(meshFile);
    nElm[i] = mesh.getNbInteriorElems();

    // feSpace *uBord, *uDomaine;
    int degreeQuadrature = 30;
    int dim;
    // feCheck(createFiniteElementSpace(uBord, &mesh, 1, LINE, order, "U", "Bord", degreeQuadrature, funSol));
    // feCheck(createFiniteElementSpace(uDomaine, &mesh, 2, TRI, order, "U", "Domaine", degreeQuadrature, funZero));
    feSpace *uSaut, *uHaut, *uBas, *uDroite, *uGauche, *uDomaine;
    feCheck(createFiniteElementSpace(uHaut,    &mesh, dim = 1, LINE, order, "U", "Haut",    degreeQuadrature, funZero));
    feCheck(createFiniteElementSpace(uBas,     &mesh, dim = 1, LINE, order, "U", "Bas",     degreeQuadrature, funZero));
    feCheck(createFiniteElementSpace(uDroite,  &mesh, dim = 1, LINE, order, "U", "Droite",  degreeQuadrature, funZero));
    feCheck(createFiniteElementSpace(uGauche,  &mesh, dim = 1, LINE, order, "U", "Gauche",  degreeQuadrature, funZero));
    feCheck(createFiniteElementSpace(uSaut,    &mesh, dim = 1, LINE, order, "U", "Saut",    degreeQuadrature, funZero));
    feCheck(createFiniteElementSpace(uDomaine, &mesh, dim = 2, TRI,  order, "U", "Domaine", degreeQuadrature, funZero));

    std::vector<feSpace *> spaces = {uSaut, uHaut, uBas, uDroite, uGauche, uDomaine};
    // std::vector<feSpace *> essentialSpaces = {uDroite, uGauche};
    std::vector<feSpace *> essentialSpaces = {uHaut, uBas, uGauche, uDroite};

    feMetaNumber metaNumber(&mesh, spaces, essentialSpaces);
    feSolution sol(&mesh, spaces, essentialSpaces, &metaNumber);

    // feBilinearForm diffU({uDomaine}, &mesh, degreeQuadrature, new feSysElm_2D_Diffusion(k, nullptr));
    // feBilinearForm sourceU({uDomaine}, &mesh, degreeQuadrature, new feSysElm_2D_Source(1.0, funSource));

    feBilinearForm diffU({uDomaine}, &mesh, degreeQuadrature, new feSysElm_2D_Diffusion(k, nullptr));
    feBilinearForm sourceU({uDomaine}, &mesh, degreeQuadrature, new feSysElm_2D_Source(1.0, funSource));
    // feBilinearForm n1({uDroite}, &mesh, degreeQuadrature, new feSysElm_1D_NeumannBC(1.0, funNeumannY));
    // feBilinearForm n2({uGauche}, &mesh, degreeQuadrature, new feSysElm_1D_NeumannBC(1.0, funZero));
    // feBilinearForm n3({uHaut},   &mesh, degreeQuadrature, new feSysElm_1D_NeumannBC(1.0, funNeumann));
    // feBilinearForm n4({uBas},    &mesh, degreeQuadrature, new feSysElm_1D_NeumannBC(1.0, funZero));
    feBilinearForm n5({uSaut},   &mesh, degreeQuadrature, new feSysElm_1D_NeumannBC(1.0, funOne));
    // feBilinearForm sautU({uSaut}, &mesh, degreeQuadrature, new feSysElm_1D_NeumannBC(1.0, funOne));

    feLinearSystem *system;
    // feCheck(createLinearSystem(system, MKLPARDISO, {&diffU, &sourceU}, &metaNumber, &mesh, argc, argv));
    // feCheck(createLinearSystem(system, MKLPARDISO, {&diffU, &n1, &n2}, &metaNumber, &mesh, argc, argv));
    // feCheck(createLinearSystem(system, MKLPARDISO, {&diffU, &sautU}, &metaNumber, &mesh, argc, argv));
    feCheck(createLinearSystem(system, MKLPARDISO, {&diffU, &n5}, &metaNumber, &mesh, argc, argv));

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
    // delete uBord;
    // delete uDirichlet;
    // delete uNeumannHaut;
    // delete uNeumannBas;
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

  return 0;
}