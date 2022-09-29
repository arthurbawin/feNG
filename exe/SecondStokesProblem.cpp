#include <iostream>
#include <fstream>

#include "feAPI.h"

double fSolU(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  double rho = par[0];
  double nu = par[1];
  double U = par[2];
  double omega = par[3];
  double rac = sqrt(omega/(2*nu));

  return U*exp(-rac*y)*cos(omega*t -rac*y);
}

double fSolV(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  double rho = par[0];
  double nu = par[1];

  return 0.0;
}

double fSolP(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  double rho = par[0];
  double nu = par[1];

  return 0.0;
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  return 0.0;
}

double fOne(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  return 1.0;
}

void fVecVeloc(const double t, const std::vector<double> &pos, const std::vector<double> &par, std::vector<double> &res) {
  double x = pos[0];
  double y = pos[1];
  double rho = par[0];
  double nu = par[1];

}

void fSource(const double t, const std::vector<double> &pos, const std::vector<double> &par, std::vector<double> &res) {
  double x = pos[0];
  double y = pos[1];
  double rho = par[0];
  double nu = par[1];
  double U = par[2];
  double omega = par[3];
  double rac = sqrt(omega/(2*nu));

  res[0] = +nu*omega*U*exp(-rac*y)*sin(omega*t-rac*y);
  res[1] = 0.0;
}

int main(int argc, char **argv) {

  petscInitialize(argc, argv);
  int orderUS=2;
  int orderUB=2;
  int orderVS=2;
  int orderVB=2;
  int orderPP=0;
  int orderPS=1;

  int nIter=1;

  feOptionsParser options(argc, argv);
  options.addOption(&nIter, "-n", "--nMeshes", "Number of meshes for the convergence study");
  feCheck(options.parse());
  setVerbose(2);

  double rho = 1000;
  double nu = 0.001;
  double U = 1;
  double omega = 1.0;
  std::vector<double> stokesParam = {rho, nu, U, omega};
  feFunction *funSolU = new feFunction(fSolU, stokesParam);
  feFunction *funSolV = new feFunction(fSolV, stokesParam);
  feFunction *funSolP = new feFunction(fSolP, stokesParam);
  feFunction *funZero = new feFunction(fZero, {});
  // feFunction *funOne = new feFunction(fOne, {});
  feVectorFunction *funSource = new feVectorFunction(fSource, stokesParam);
  feVectorFunction *funVecVeloc = new feVectorFunction(fVecVeloc, stokesParam);

  std::vector<double> normL2_U(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);
  std::vector<int> TT(nIter, 0);

  for(int iter = 0; iter < nIter; ++iter) {    
    std::string meshName = "../data/Cavite" + std::to_string(iter+1) + ".msh";
    // std::string meshName = "KelvinCatEye" + std::to_string(iter+1) + ".msh";

    feMesh2DP1 mesh(meshName);
    nElm[iter] = mesh.getNbInteriorElems();

    int degreeQuadrature = 10;

    feSpace *P_point, *P_surface;
    feSpace *U_haut, *U_bas, *U_droite, *U_gauche, *U_surface;
    feSpace *V_haut, *V_bas, *V_droite, *V_gauche,*V_surface;

    feCheck(createFiniteElementSpace(U_surface, &mesh, 2, TRI, orderUS, "U", "Surface", degreeQuadrature, funSolU));
    feCheck(createFiniteElementSpace(U_haut, &mesh, 1, LINE, orderUB, "U", "Haut", degreeQuadrature, funSolU));     
    feCheck(createFiniteElementSpace(U_bas, &mesh, 1, LINE, orderUB, "U", "Bas", degreeQuadrature, funSolU));
    feCheck(createFiniteElementSpace(U_droite, &mesh, 1, LINE, orderUB, "U", "Droite", degreeQuadrature, funSolU));
    feCheck(createFiniteElementSpace(U_gauche, &mesh, 1, LINE, orderUB, "U", "Gauche", degreeQuadrature, funSolU));     

    feCheck(createFiniteElementSpace(V_surface, &mesh, 2, TRI, orderVS, "V", "Surface", degreeQuadrature, funSolV)); 
    feCheck(createFiniteElementSpace(V_haut, &mesh, 1, LINE, orderVB, "V", "Haut", degreeQuadrature, funSolV));   
    feCheck(createFiniteElementSpace(V_bas, &mesh, 1, LINE, orderVB, "V", "Bas", degreeQuadrature, funSolV));
    feCheck(createFiniteElementSpace(V_droite, &mesh, 1, LINE, orderVB, "V", "Droite", degreeQuadrature, funSolV));
    feCheck(createFiniteElementSpace(V_gauche, &mesh, 1, LINE, orderVB, "V", "Gauche", degreeQuadrature, funSolV));      

    feCheck(createFiniteElementSpace(P_point, &mesh, 0, POINT, orderPP, "P", "PointPression", degreeQuadrature, funSolP));
    feCheck(createFiniteElementSpace(P_surface, &mesh, 2, TRI, orderPS, "P", "Surface", degreeQuadrature, funSolP));

    std::vector<feSpace *> spaces = {U_haut, U_bas, U_droite, U_gauche, U_surface, V_haut, V_bas, V_droite, V_gauche, V_surface, P_point, P_surface};
    std::vector<feSpace *> feEssBC = {U_haut, U_bas, U_droite, U_gauche, V_haut, V_bas, V_droite, V_gauche, V_surface, P_point, P_surface};

    feMetaNumber metaNumber(&mesh, spaces, feEssBC);

    feSolution sol(&mesh, spaces, feEssBC, &metaNumber);

    // Formes (bi)lineaires
    std::vector<feSpace *> spacesNS2D = {U_surface, V_surface, P_surface};
    feBilinearForm NS2D(spacesNS2D, &mesh, degreeQuadrature, new feSysElm_2D_NavierStokes(stokesParam, funSource));
  
    //System
    feLinearSystem *system;
    feCheck(createLinearSystem(system, PETSC, spaces, {&NS2D}, &metaNumber, &mesh, argc, argv));
    // feCheck(createLinearSystem(system, MKLPARDISO, spaces, {&NS2D}, &metaNumber, &mesh, argc, argv));

    //Norm
    std::vector<feSpace *> SpaceNorm = {U_surface, V_surface};
    feComputer *normU = new feComputer(U_surface, &mesh, &metaNumber, "L2Norm_1Field", funSolU);
    std::vector<feComputer *> norms = {normU};

    //Exporter
    feExporter *exporter;
    feCheck(createVisualizationExporter(exporter, VTK, &metaNumber, &sol, &mesh, spaces));
    int exportEveryNSteps = 1;
    std::string vtkFileRoot = "outNS_" + std::to_string(iter) + "_";
    feExportData exportData = {exporter, exportEveryNSteps, vtkFileRoot};

    //Resolution
    TimeIntegrator *solver;
    feTolerances tol{1e-9, 1e-9, 100};
    double t0 = 0;
    double t1 = 1;
    int nTimeSteps = 10*pow(2, iter);
    TT[iter]=nTimeSteps;
    // std::string CodeIni = "BDF1/DCF";
    std::string CodeIni = " ";
    feCheck(createTimeIntegrator(solver, BDF1, tol, system, &metaNumber, &sol, &mesh, norms, exportData, t0, t1, nTimeSteps, CodeIni));
    feCheck(solver->makeSteps(nTimeSteps));

    std::vector<double> &normL2_normU = normU->getResult();
    normL2_U[2*iter] = *std::max_element(normL2_normU.begin(), normL2_normU.end());

    delete solver;
    delete exporter;
    delete U_haut; 
    delete U_bas; 
    delete U_droite; 
    delete U_gauche; 
    delete U_surface;
    delete V_haut; 
    delete V_bas; 
    delete V_droite; 
    delete V_gauche;
    delete V_surface; 
    delete P_surface;
    delete P_point;
  }


  // Calcul du taux de convergence
  for(int i = 1; i < nIter; ++i) {
    normL2_U[2 * i + 1] = -log(normL2_U[2 * i] / normL2_U[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
  }

  printf("%12s\t %12s \t %12s \t %12s \n", "nSteps", "nElm", "||u-uh||", "tauxU");
  
  for(int i = 0; i < nIter; ++i)
    printf("%12d\t %12d \t %12.6e \t %12.6e \n", TT[i], nElm[i], normL2_U[2 * i], normL2_U[2 * i + 1]);

 
  delete funSource;
  delete funSolU;
  delete funSolV;
  delete funSolP;


  petscFinalize();

  return 0;
}
