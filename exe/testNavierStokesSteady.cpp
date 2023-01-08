#include <iostream>
#include <fstream>

#include "feAPI.h"

double fSolU(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  double F = exp(-2.0 * t);
  return F * sin(x) * cos(y); // Taylor-Green
  // return pow(x,4) + pow(y,4);
}

double fSolV(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  double F = exp(-2.0 * t);
  return -F * cos(x) * sin(y); // Taylor-Green
  // return -4.0 * pow(x,3) * y;
}

double fSolP(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  // double rho = par[0];
  double F = exp(-2.0 * t);
  return F * F * 1.0 / 4.0 * (cos(2. * x) + cos(2. * y)); // Taylor-Green
  // return pow(x,4) * pow(y,4);
}

double fOne(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  return 1.0;
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

  // res[0] = -mu * (12*(x*x*y*y*y*y + x*x*x*x*y*y)) + 2*x*y*y;
  // res[1] = -mu * (-8./5. *x*y*y*y*(10*x*x + 3*y*y)) + 2*y*x*x;

  // res[0] = 4.0*pow(x,3)*(pow(x,4) + pow(y,4)) - 12.0*pow(x,3)*pow(y,4) - 12.0*x*x - 12.0*y*y;
  // res[1] = 4.0*pow(x,4)*pow(y,3) + 24.0*x*y + 16.0*pow(x,6)*y - 12.0*x*x*y*(pow(x,4) + pow(y,4));
  if(fabs(t) <= 1e-10) {
    res[0] = 2.0 * cos(y) * sin(x);
    res[1] = -2.0 * cos(x) * sin(y);
  } else {
    res[0] = 0.0;
    res[1] = 0.0;
  }
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

  int nIter=1;

  feOptionsParser options(argc, argv);
  options.addOption(&nIter, "-n", "--nMeshes", "Number of meshes for the convergence study");
  feCheck(options.parse());
  setVerbose(2);

  double rho = 1.0;
  double mu = 1.0;
  std::vector<double> stokesParam = {rho, mu};
  feFunction *funSolU = new feFunction(fSolU, {});
  feFunction *funSolV = new feFunction(fSolV, {});
  feFunction *funSolP = new feFunction(fSolP, stokesParam);
  feFunction *funZero = new feFunction(fZero, {});
  // feFunction *funOne = new feFunction(fOne, {});
  feVectorFunction *funSource = new feVectorFunction(fSource, stokesParam);

  std::vector<double> normL2_U(2 * nIter, 0.0);
  std::vector<double> normL2_V(2 * nIter, 0.0);
  std::vector<double> normL2_P(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);

  for(int iter = 0; iter < nIter; ++iter) {
    // std::string meshName = "../../data/squareTaylorGreen2.msh";
    std::string meshName = "../data/Convergence/squareNS" + std::to_string(iter + 1) + ".msh";
    // std::string meshName = "../data/Convergence/squarePression" + std::to_string(iter+1) + ".msh";

    feMesh2DP1 mesh(meshName);
    nElm[iter] = mesh.getNbInteriorElems();

    int degreeQuadrature = 5;
    feSpace *U_angle, *U_haut, *U_gauche, *U_surface;
    feSpace *V_angle, *V_haut, *V_gauche, *V_surface;
    feSpace *P_haut, *P_surface;
    // feSpace *P_point, *P_surface;

    feCheck(createFiniteElementSpace(U_angle, &mesh, 1, LINE, orderUA, "U", "Angle", degreeQuadrature, funSolU));
    feCheck(createFiniteElementSpace(U_haut, &mesh, 1, LINE, orderUH, "U", "Haut", degreeQuadrature, funSolU));
    feCheck(createFiniteElementSpace(U_gauche, &mesh, 1, LINE, orderUG, "U", "Gauche", degreeQuadrature, funSolU));
    feCheck(createFiniteElementSpace(U_surface, &mesh, 2, TRI, orderUS, "U", "Surface", degreeQuadrature, funZero));
    
    feCheck(createFiniteElementSpace(V_angle, &mesh, 1, LINE, orderVA, "V", "Angle", degreeQuadrature, funSolV));
    feCheck(createFiniteElementSpace(V_haut, &mesh, 1, LINE, orderVH, "V", "Haut", degreeQuadrature, funSolV));
    feCheck(createFiniteElementSpace(V_gauche, &mesh, 1, LINE, orderVG, "V", "Gauche", degreeQuadrature, funSolV));
    feCheck(createFiniteElementSpace(V_surface, &mesh, 2, TRI, orderVS, "V", "Surface", degreeQuadrature, funZero));

    feCheck(createFiniteElementSpace(P_haut, &mesh, 1, LINE, orderPH, "P", "Haut", degreeQuadrature, funSolP));
    // feCheck(createFiniteElementSpace(P_point, &mesh, 0, POINT, orderPH, "P", "PointPression", degreeQuadrature, funSolP));
    feCheck(createFiniteElementSpace(P_surface, &mesh, 2, TRI, orderPS, "P", "Surface", degreeQuadrature, funZero));

    std::vector<feSpace *> spaces = {U_angle, U_haut, U_gauche, U_surface, V_angle, V_haut, V_gauche, V_surface, P_haut, P_surface};
    // std::vector<feSpace *> spaces = {U_angle, U_haut, U_gauche, U_surface, V_angle, V_haut, V_gauche, V_surface, P_point, P_surface};
    std::vector<feSpace *> feEssBC = {U_angle, U_haut, U_gauche, V_angle, V_haut, V_gauche, P_haut};
    // std::vector<feSpace *> feEssBC = {U_angle, U_haut, U_gauche, V_angle, V_haut, V_gauche, P_point};

    feMetaNumber metaNumber(&mesh, spaces, feEssBC);

    feSolution sol(&mesh, spaces, feEssBC, &metaNumber);

    // Formes (bi)lineaires
    std::vector<feSpace *> spacesNS2D = {U_surface, V_surface, P_surface};
    feBilinearForm NS2D(spacesNS2D, &mesh, degreeQuadrature, new feSysElm_2D_NavierStokes(stokesParam, funSource));
  
    //System
    feLinearSystem *system;
    feCheck(createLinearSystem(system, PETSC, spaces, {&NS2D}, &metaNumber, &mesh, argc, argv));

    //Norm
    feNorm normU(U_surface, &mesh, degreeQuadrature, funSolU);
    feNorm normV(V_surface, &mesh, degreeQuadrature, funSolV);
    feNorm normP(P_surface, &mesh, degreeQuadrature, funSolP);
    std::vector<feNorm *> norms = {&normU, &normV, &normP};

    //Exporter
    feExporter *exporter;
    feCheck(createVisualizationExporter(exporter, VTK, &metaNumber, &sol, &mesh, spaces));
    int exportEveryNSteps = 1;
    std::string vtkFileRoot = "outNS" + std::to_string(iter);
    feExportData exportData = {exporter, exportEveryNSteps, vtkFileRoot};

    //Resolution
    TimeIntegrator *solver;
    feTolerances tol{1e-9, 1e-8, 10};
    feCheck(createTimeIntegrator(solver, STATIONARY, tol, system, &metaNumber, &sol, &mesh, norms, exportData));
    feCheck(solver->makeSteps(0));

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
    // delete P_point;
    // delete system;
  }
  // Calcul du taux de convergence
  for(int i = 1; i < nIter; ++i) {
    normL2_U[2 * i + 1] =
      -log(normL2_U[2 * i] / normL2_U[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
    normL2_V[2 * i + 1] =
      -log(normL2_V[2 * i] / normL2_V[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
    normL2_P[2 * i + 1] =
      -log(normL2_P[2 * i] / normL2_P[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
  }

  printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s\n", "nElm",
         "||u-uh||", "tauxU", "||v-vh||", "tauxV", "||p-ph||", "tauxP");
  
  for(int i = 0; i < nIter; ++i)
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
