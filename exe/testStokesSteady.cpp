#include <iostream>
#include <fstream>

#include "feAPI.h"

double fSolU(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  // return 2. * exp(x) * (x - 1.) * (x - 1.) * x * x * (y * y - y) * (2. * y - 1.);
  return x * x * x * x * y * y * y * y;
}

double fSolV(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  // return exp(x) * (x - 1.) * x * (x * x + 3. * x - 2.) * (y * y - 1.) * y * y;
  return -4. / 5. * x * x * x * y * y * y * y * y;
}

double fSolP(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  // return -424. + 156. * exp(1) +
  //        (y * y - y) *
  //          (-456. +
  //           exp(x) * (456. + x * x * (228. - 5. * (y * y - y)) + 2. * x * (y * y - y - 228.) +
  //                     2. * x * x * x * (y * y - y - 36.) + x * x * x * x * (y * y - y + 12.)));
  return x * x * y * y;
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  return 0.0;
}

void fSource(const double t, const std::vector<double> &pos, const std::vector<double> &par,
             std::vector<double> &res) {
  // double rho = par[0];
  double mu = par[1];
  // double nu = mu/rho;
  double x = pos[0];
  double y = pos[1];
  // double dudx = 2*exp(x)*x*y*(2*y*y-3*y+1)*(x*x*x + 2*x*x - 5*x + 2);
  // double dudy = 2*exp(x)*x*x*(x-1)*(x-1)*(6*y*y-6*y+1);
  // double dvdx = exp(x)*y*y*(y-1)*(y-1)*(x*x*x*x + 6*x*x*x + x*x -8*x+2); // Signe ? Pour avoir div = 0 
  // double dvdy = 2*exp(x)*x*y*(2*y*y-3*y+1)*(x*x*x + 2*x*x - 5*x + 2); 
  // double d2udx2 =2*exp(x)*y*(x*x*x*x + 6*x*x*x + x*x -8*x + 2) * (2*y*y - 3*y + 1); 
  // double d2udxy =2*exp(x)*x*(x*x*x + 2*x*x - 5*x + 2) * (6*y*y-6*y+1); 
  // double d2udy2 = 12*exp(x)*x*x*(x-1)*(x-1)*(2*y-1); double d2vdx2 = -exp(x)*y*y*(x*x*x*x + 10*x*x*x + 19*x*x - 6*x - 6)*(y-1)*(y-1);
  // double d2vdxy = -2*exp(x)*y*(x*x*x*x + 6*x*x*x + x*x -8*x + 2) * (2*y*y -3*y + 1); 
  // double d2vdy2 = -2*exp(x)*x*(x-1)*(x*x+3*x-2)*(6*y*y-6*y+1); 
  // double ddxxdx = 2 *d2udx2; 
  // double ddyxdy = d2udy2 + d2vdxy;
  // double ddxydx = d2udxy + d2vdx2; 
  // double ddyydy = 2 *d2vdy2;

  // double dpdx = (y*y - y)*(exp(x)*(4*x*x*x*(y*y - y + 12) + 6*x*x*(y*y - y - 36) + 2*(228 -5*(y*y - y)) + 2*(y*y - y - 228))+ exp(x)*(x*x*x*x*(y*y - y + 12) + 2*x*x*x*(y*y - y - 36) + x*x*(228 - 5*(y*y - y)) + 2*x*(y*y- y - 228) + 456));
  // double dpdy = exp(x)*(y*y - y)*(x*x*x*x*(2*y - 1) + 2*x*x*x* (2*y - 1) - 5*x*x*(2*y - 1) +2*x*(2*y - 1))+ (2*y - 1)*(exp(x)*(x*x*x*x*(y*y - y + 12) + 2*x*x*x*(y*y - y - 36) + x*x*(228 - 5*(y*y - y))+ 2*x*(y*y - y - 228) + 456) - 456);

  // -div (nu*d) + grad(p)
  // res[0] = - (ddxxdx + ddyxdy) + dpdx;
  // res[1] = - (ddxydx + ddyydy) + dpdy;

  res[0] = -mu * (12 * (x * x * y * y * y * y + x * x * x * x * y * y)) + 2 * x * y * y;
  res[1] = -mu * (-8. / 5. * x * y * y * y * (10 * x * x + 3 * y * y)) + 2 * y * x * x;
}

int main(int argc, char **argv) {

  petscInitialize(argc, argv);
  int orderUA=2;
  int orderUH=2;
  int orderUG=2;
  int orderUS=2;
  int orderVA=2;
  int orderVH=2;
  int orderVG=2;
  int orderVS=2;
  int orderPH=1;
  int orderPS=1;
  int nMesh = 1;

  feOptionsParser options(argc, argv);
  options.addOption(&nMesh, "-n", "--nMeshes", "Number of meshes for the convergence study");
  feCheck(options.parse());

  setVerbose(2);

  double rho = 1.0;
  double mu = 1.0;
  std::vector<double> stokesParam = {rho, mu};
  feFunction *funSolU = new feFunction(fSolU, {});
  feFunction *funSolV = new feFunction(fSolV, {});
  feFunction *funSolP = new feFunction(fSolP, {});
  feFunction *funZero = new feFunction(fZero, {});
  feVectorFunction *funSource = new feVectorFunction(fSource, stokesParam);

  std::vector<double> normL2_U(2 * nMesh, 0.0);
  std::vector<double> normL2_V(2 * nMesh, 0.0);
  std::vector<double> normL2_P(2 * nMesh, 0.0);
  std::vector<int> nElm(nMesh, 0);

  for(int iter = 0; iter < nMesh; ++iter) {
    // std::string meshName = "../../data/squarePression" + std::to_string(iter+1) + ".msh";
    std::string meshName = "../data/Convergence/squareStokes" + std::to_string(iter + 1) + ".msh";
    // std::string meshName = "../../data/squareCoarse" + std::to_string(iter+1) + ".msh";

    feMesh2DP1 mesh(meshName);
    nElm[iter] = mesh.getNbInteriorElems();

    int degreeQuadrature = 5;
    int dim;
    feSpace *U_angle, *U_haut, *U_gauche, *U_surface;
    feSpace *V_angle, *V_haut, *V_gauche, *V_surface;
    feSpace *P_haut, *P_surface;

    feCheck(createFiniteElementSpace(U_angle, &mesh, 1, LINE, orderUA, "U", "Angle", degreeQuadrature, funSolU));
    feCheck(createFiniteElementSpace(U_haut, &mesh, 1, LINE, orderUH, "U", "Haut", degreeQuadrature, funSolU));
    feCheck(createFiniteElementSpace(U_gauche, &mesh, 1, LINE, orderUG, "U", "Gauche", degreeQuadrature, funSolU));
    feCheck(createFiniteElementSpace(U_surface, &mesh, 2, TRI, orderUS, "U", "Surface", degreeQuadrature, funZero));

    feCheck(createFiniteElementSpace(V_angle, &mesh, 1, LINE, orderVA, "V", "Angle", degreeQuadrature, funSolV));
    feCheck(createFiniteElementSpace(V_haut, &mesh, 1, LINE, orderVH, "V", "Haut", degreeQuadrature, funSolV));
    feCheck(createFiniteElementSpace(V_gauche, &mesh, 1, LINE, orderVG, "V", "Gauche", degreeQuadrature, funSolV));
    feCheck(createFiniteElementSpace(V_surface, &mesh, 2, TRI, orderVS, "V", "Surface", degreeQuadrature, funZero));

    feCheck(createFiniteElementSpace(P_haut, &mesh, 1, LINE, orderPH, "P", "Haut", degreeQuadrature, funSolP));
    feCheck(createFiniteElementSpace(P_surface, &mesh, 2, TRI, orderPS, "P", "Surface", degreeQuadrature, funZero));

    std::vector<feSpace *> spaces = {U_angle, U_haut, U_gauche, U_surface, V_angle, V_haut, V_gauche, V_surface, P_surface, P_haut};
    std::vector<feSpace *> feEssBC = {U_angle, U_haut, U_gauche, V_angle, V_haut, V_gauche, P_haut};

    feMetaNumber metaNumber(&mesh, spaces, feEssBC);

    feSolution sol(&mesh, spaces, feEssBC, &metaNumber);

    // Formes (bi)lineaires
    std::vector<feSpace *> spacesStokes2D = {U_surface, V_surface, P_surface};
    feBilinearForm stokes2D(spacesStokes2D, &mesh, degreeQuadrature, new feSysElm_2D_Stokes(stokesParam, funSource));

    feLinearSystem *system;
    feCheck(createLinearSystem(system, PETSC, spaces, {&stokes2D}, &metaNumber, &mesh, argc, argv));

    //Norms
    feComputer normU(U_surface, &mesh, &metaNumber, "L2Norm_1Field", funSolU);
    feComputer normV(V_surface, &mesh, &metaNumber, "L2Norm_1Field", funSolV);
    feComputer normP(P_surface, &mesh, &metaNumber, "L2Norm_1Field", funSolP);
    std::vector<feComputer *> norms = {&normU, &normV, &normP};

    //Exporter
    feExporter *exporter;
    feCheck(createVisualizationExporter(exporter, VTK, &metaNumber, &sol, &mesh, spaces));
    int exportEveryNSteps = 1;
    std::string vtkFileRoot = "out" + std::to_string(iter);
    feExportData exportData = {exporter, exportEveryNSteps, vtkFileRoot};

    //Resolution
    TimeIntegrator *solver;
    feTolerances tol{1e-9, 1e-8, 10};
    feCheck(createTimeIntegrator(solver, STATIONARY, tol, system, &metaNumber, &sol, &mesh, norms, exportData));
    feCheck(solver->makeSteps(0));

    //Post-porcess
    fePostProc postU(U_surface, &mesh, &metaNumber, funSolU);
    fePostProc postV(V_surface, &mesh, &metaNumber, funSolV);
    fePostProc postP(P_surface, &mesh, &metaNumber, funSolP);
    normL2_U[2 * iter] = postU.computeL2ErrorNorm(&sol);
    normL2_V[2 * iter] = postV.computeL2ErrorNorm(&sol);
    normL2_P[2 * iter] = postP.computeL2ErrorNorm(&sol);

    delete solver;
    delete exporter;
    delete U_surface;
    delete U_haut;
    delete U_angle;
    delete U_gauche;
    delete V_surface;
    delete V_haut;
    delete V_angle;
    delete V_gauche;
    delete P_surface;
    delete P_haut;
  }
  // Calcul du taux de convergence
  for(int i = 1; i < nMesh; ++i) {
    normL2_U[2 * i + 1] =
      -log(normL2_U[2 * i] / normL2_U[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
    normL2_V[2 * i + 1] =
      -log(normL2_V[2 * i] / normL2_V[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
    normL2_P[2 * i + 1] =
      -log(normL2_P[2 * i] / normL2_P[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
  }
  printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s\n", "nElm", "||u-uh||", "tauxU",
         "||v-vh||", "tauxV", "||p-ph||", "tauxP");
  for(int i = 0; i < nMesh; ++i)
    printf("%12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n", nElm[i],
           normL2_U[2 * i], normL2_U[2 * i + 1], normL2_V[2 * i], normL2_V[2 * i + 1],
           normL2_P[2 * i], normL2_P[2 * i + 1]);

  delete funSource;
  delete funSolU;
  delete funSolV;
  delete funSolP;
  delete funZero;
  
  petscFinalize();

  return 0;
}
