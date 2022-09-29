#include <iostream>
#include <functional>
#include <vector>
#include <cmath>

#include "feMesh.h"
#include "feQuadrature.h"
#include "feSpace.h"
#include "feCncGeo.h"
#include "feNumber.h"
#include "feSolution.h"
#include "feSysElm.h"
#include "feBilinearForm.h"
#include "feSolver.h"
#define USING_PETSC
#ifdef HAVE_PETSC

#include "feLinearSystemPETSc.h"
#endif

// Test a source nulle
double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  return (x1 * x1 - 25.) * pow(t, 6) + pow(x1, 2) * pow(t, 6) + 1.;
}

double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  return 6.*(x1 * x1 - 25.)  * pow(t, 5.) + 6.*x1*x1*pow(t,5) -4.*pow(t,6);
}

double flambda_A(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double c1 = par[0];
  // int a = par[1];
  return 4.*x1*pow(t,6);
}
double flambda_B(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double c1 = par[0];
  return -4.*x1*pow(t,6);
}

double fSolDot(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  return 6*(x1*x1-25)*pow(t,5) + 6*x1*x1*pow(t,5);
}

int main(int argc, char **argv)
{

  petscInitialize(argc, argv);

  double xa = 0.;
  double xb = 5.;

  double kd = 1;
  std::vector<double> par = {kd, 6.};
  feFunction *funSol = new feFunction(fSol, par);
  feFunction *funSolDot = new feFunction(fSolDot, par);
  feFunction *funSource = new feFunction(fSource, par);
  feFunction *funLambda_A = new feFunction(flambda_A, par);
  feFunction *funLambda_B = new feFunction(flambda_B, par);

  int nIter = 1;
  // std::vector<double> normL2_BDF2(2 * nIter, 0.0);
  // std::vector<double> normL2_DC3(2 * nIter, 0.0);
  // std::vector<double> normBDF2_A(2 * nIter, 0.0);
  // std::vector<double> normBDF2_B(2 * nIter, 0.0);
  // std::vector<double> normDC3_A(2 * nIter, 0.0);
  // std::vector<double> normDC3_B(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);
  std::vector<int> TT;
  TT.resize(nIter);

  int degreeQuadrature = 10;
  for(int iter = 0; iter < nIter; ++iter) {
    nElm[iter] = 3 * pow(2, iter);

    feMesh1DP1 mesh(xa, xb, nElm[iter], "BXA", "BXB", "M1D");
    

    int dim;
    feSpace *U_BXA, *U_BXB, *U_M1D, *L_BXA, *L_BXB, *V_BXA, *V_BXB;
    feCheck(createFiniteElementSpace(U_BXA, &mesh, dim=0, POINT, 0, "U", "BXA", degreeQuadrature, funSol));
    feCheck(createFiniteElementSpace(U_BXB, &mesh, dim=0, POINT, 0, "U", "BXB", degreeQuadrature, funSol));
    feCheck(createFiniteElementSpace(U_M1D, &mesh, dim=1, LINE, 1, "U", "M1D", degreeQuadrature, funSol));
    feCheck(createFiniteElementSpace(L_BXA, &mesh, dim=0, POINT, 0, "L", "BXA", degreeQuadrature, funLambda_A));
    feCheck(createFiniteElementSpace(L_BXB, &mesh, dim=0, POINT, 0, "L", "BXB", degreeQuadrature, funLambda_B));
    feCheck(createFiniteElementSpace(V_BXA, &mesh, dim=0, POINT, 0, "V", "BXA", degreeQuadrature, funSol));
    feCheck(createFiniteElementSpace(V_BXB, &mesh, dim=0, POINT, 0, "V", "BXB", degreeQuadrature, funSol));
    
    std::vector<feSpace *> fespace = {U_BXA, U_BXB, U_M1D, L_BXA, L_BXB, V_BXA, V_BXB};
    std::vector<feSpace *> feEssBC = {};
    
    // Numerotations
    feMetaNumber metaNumber(&mesh, fespace, feEssBC);

    // Solution
    feSolution sol(&mesh, fespace, feEssBC, &metaNumber);
    
    // Formes (bi)lineaires
    std::vector<feSpace *> spaceWeak_A = {U_BXA, L_BXA, V_BXA};
    std::vector<feSpace *> spaceWeak_B = {U_BXB, L_BXB, V_BXB};

    feBilinearForm diff({U_M1D}, &mesh, degreeQuadrature, new feSysElm_1D_Diffusion(kd, nullptr));
    feBilinearForm source({U_M1D}, &mesh, degreeQuadrature, new feSysElm_1D_Source(1.0, funSource));
    feBilinearForm masse({U_M1D}, &mesh, degreeQuadrature, new feSysElm_1D_Masse(1.0, nullptr));
    feBilinearForm weakBC_A(spaceWeak_A, &mesh, degreeQuadrature, new feSysElm_0D_weakBC_edo1(1.0, funSolDot));
    feBilinearForm weakBC_B(spaceWeak_B, &mesh, degreeQuadrature, new feSysElm_0D_weakBC_edo1(1.0, funSolDot));

    // Systeme lineaire
    feLinearSystem *system;
    // feCheck(createLinearSystem(system, MKLPARDISO, spaces, {&NS2D, &weakBC_UCylinder, &weakBC_VCylinder},&metaNumber, &mesh, argc, argv));
    feCheck(createLinearSystem(system, PETSC, fespace, {&diff, &source, &masse, &weakBC_A, &weakBC_B}, &metaNumber, &mesh, argc, argv));

    // Norme de la solution
    // feNorm *normBDF2 = new feNorm(&U_M1D, mesh, degQuad, funSol);
    // feNorm *normDC3 = new feNorm(&U_M1D, mesh, degQuad, funSol);
    // feNorm *normBDF2Lambda_A = new feNorm(&L_BXA, mesh, degQuad, funLambda_A);
    // feNorm *normBDF2Lambda_B = new feNorm(&L_BXB, mesh, degQuad, funLambda_B);
    // feNorm *normDC3Lambda_A = new feNorm(&L_BXA, mesh, degQuad, funLambda_A);
    // feNorm *normDC3Lambda_B = new feNorm(&L_BXB, mesh, degQuad, funLambda_B);
    std::vector<feComputer *> comput = {};
    
    //exporter
    feExporter *exporter;
    feCheck(createVisualizationExporter(exporter, VTK, &metaNumber, &sol, &mesh, fespace));
    int exportEveryNSteps = 1;
    std::string vtkFileRoot = "../../data/ResultatVK/VonKarman";
    feExportData exportData = {nullptr, exportEveryNSteps, vtkFileRoot};


    TimeIntegrator *solver;
    feTolerances tol{1e-9, 1e-9, 0};
    double t0 = 0.;
    double t1 = 1.;
    std::string CodeIni = "";
    int nTimeSteps = 10 * pow(2, iter);
    TT[iter] = nTimeSteps;

    feCheck(createTimeIntegrator(solver, BDF1, tol, system, &metaNumber, &sol, &mesh, comput, exportData, t0, t1, nTimeSteps, CodeIni));
    feCheck(solver->makeSteps(nTimeSteps));
    


    // std::vector<double> &normL2BDF2 = solver.getNorm(0);
    // std::vector<double> &normL2DC3 = solver.getNorm(1);
    // std::vector<double> &normL2_BDF2_Lambda_A = solver.getNorm(2);
    // std::vector<double> &normL2_BDF2_Lambda_B = solver.getNorm(3);
    // std::vector<double> &normL2_DC3_Lambda_A = solver.getNorm(4);
    // std::vector<double> &normL2_DC3_Lambda_B = solver.getNorm(5);
    // normL2_BDF2[2 * iter] = *std::max_element(normL2BDF2.begin(), normL2BDF2.end());
    // normL2_DC3[2 * iter] = *std::max_element(normL2DC3.begin(), normL2DC3.end());
    // normBDF2_A[2 * iter] =
    //   *std::max_element(normL2_BDF2_Lambda_A.begin(), normL2_BDF2_Lambda_A.end());
    // normBDF2_B[2 * iter] =
    //   *std::max_element(normL2_BDF2_Lambda_A.begin(), normL2_BDF2_Lambda_A.end());
    // normDC3_A[2 * iter] = *std::max_element(normL2_DC3_Lambda_A.begin(), normL2_DC3_Lambda_A.end());
    // normDC3_B[2 * iter] = *std::max_element(normL2_DC3_Lambda_A.begin(), normL2_DC3_Lambda_A.end());


    delete U_BXA;
    delete U_BXB;
    delete U_M1D; 
    delete L_BXA; 
    delete L_BXB;
    delete V_BXA;
    delete V_BXB;
  }
  delete funSource;
  delete funSol;
  delete funSolDot;
  delete funLambda_A;
  delete funLambda_B;

  // Calcul du taux de convergence
//   for(int i = 1; i < nIter; ++i) {
//     normL2_BDF2[2 * i + 1] = log(normL2_BDF2[2 * (i - 1)] / normL2_BDF2[2 * i]) / log(2.);
//     normL2_DC3[2 * i + 1] = log(normL2_DC3[2 * (i - 1)] / normL2_DC3[2 * i]) / log(2.);
//     normBDF2_A[2 * i + 1] = log(normBDF2_A[2 * (i - 1)] / normBDF2_A[2 * i]) / log(2.);
//     normBDF2_B[2 * i + 1] = log(normBDF2_B[2 * (i - 1)] / normBDF2_B[2 * i]) / log(2.);
//     normDC3_A[2 * i + 1] = log(normDC3_A[2 * (i - 1)] / normDC3_A[2 * i]) / log(2.);
//     normDC3_B[2 * i + 1] = log(normDC3_B[2 * (i - 1)] / normDC3_B[2 * i]) / log(2.);
//   }
//   printf("%12s \t %12s \t %12s \t %12s\t %12s \t %12s\n", "nSteps", "nElm", "||E-BDF||",
//          "Taux BDF2", "||E-DC3||", "Taux DC3");
//   for(int i = 0; i < nIter; ++i)

//     printf("%12d \t %12d \t %12.6e \t %12.6e\t %12.6e \t %12.6e\n", TT[i], nElm[i],
//            normL2_BDF2[2 * i], normL2_BDF2[2 * i + 1], normL2_DC3[2 * i], normL2_DC3[2 * i + 1]);
//   printf("\n");
//   printf("\n");
//   printf("Lambda based on BDF2");
//   printf("\n");
//   printf("%12s \t %12s \t %16s \t %16s \t %16s \t %16s \n", "nTimeSteps", "nElm", "||E_lambda_A||",
//          "Taux lambda_A", "||E_lambda_B||", "Taux lambda_B");
//   for(int i = 0; i < nIter; ++i)
//     printf("%12d \t %12d \t %16.6e \t %16.6e \t %16.6e \t %16.6e  \n", TT[i], nElm[i],
//            normBDF2_A[2 * i], normBDF2_A[2 * i + 1], normBDF2_B[2 * i], normBDF2_B[2 * i + 1]);
//   printf("\n");
//   printf("\n");
//   printf("Lambda based on DC3");
//   printf("\n");
//   printf("%12s \t %12s \t %16s \t %16s \t %16s \t %16s \n", "nTimeSteps", "nElm", "||E_lambda_A||",
//          "Taux lambda_A", "||E_lambda_B||", "Taux lambda_B");
//   for(int i = 0; i < nIter; ++i)
//     printf("%12d \t %12d \t %16.6e \t %16.6e \t %16.6e \t %16.6e  \n", TT[i], nElm[i],
//            normDC3_A[2 * i], normDC3_A[2 * i + 1], normDC3_B[2 * i], normDC3_B[2 * i + 1]);
// #ifdef USING_PETSC
//   petscFinalize();
// #endif
  return 0;
}