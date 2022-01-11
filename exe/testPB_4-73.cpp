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

double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  return x1 * x1 * x1 * (x1 * x1 - 25.) * pow(t, 6) + pow(x1, 5) * pow(t, 6) + 1.;
}

double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double c1 = par[0];
  double x1 = x[0];
  return -x1 * x1 * x1 * (x1 * x1 - 25.) * 6. * pow(t, 5.) - 6. * pow(x1, 5) * pow(t, 5) -
         c1 * 150. * x1 * pow(t, 6.) + c1 * 40. * x1 * x1 * x1 * pow(t, 6);
}

double flambda_A(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double c1 = par[0];
  return c1 * (10. * pow(x1, 4.) * pow(t, 6) - 75. * x1 * x1 * pow(t, 6));
}
double flambda_B(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double c1 = par[0];
  return -c1 * (10. * pow(x1, 4.) * pow(t, 6) - 75. * x1 * x1 * pow(t, 6));
}

int main(int argc, char **argv) {
#ifdef USING_PETSC
  petscInitialize(argc, argv);
#endif
  double xa = 0.;
  double xb = 5.;
  double kd = 0.1;
  std::vector<double> par = {kd, 0.};
  feFunction *funSol = new feFunction(fSol, par);
  feFunction *funSource = new feFunction(fSource, par);
  feFunction *funLambda_A = new feFunction(flambda_A, par);
  feFunction *funLambda_B = new feFunction(flambda_B, par);

  int nIter = 4;
  std::vector<double> normL2_BDF2(2 * nIter, 0.0);
  std::vector<double> normL2_DC3(2 * nIter, 0.0);
  std::vector<double> normBDF2_A(2 * nIter, 0.0);
  std::vector<double> normBDF2_B(2 * nIter, 0.0);
  std::vector<double> normDC3_A(2 * nIter, 0.0);
  std::vector<double> normDC3_B(2 * nIter, 0.0);

  std::vector<int> nElm(nIter, 0);
  std::vector<int> TT;
  TT.resize(nIter);

  for(int iter = 0; iter < nIter; ++iter) {
    nElm[iter] = 40 * pow(2, iter);
    // Maillage
    feMesh1DP1 *mesh = new feMesh1DP1(xa, xb, nElm[iter], "BXA", "BXB", "M1D");
    // Espaces d'interpolation
    feSpace1DP0 U_BXA(mesh, "U", "BXA", funSol);
    feSpace1DP0 U_BXB(mesh, "U", "BXB", funSol);
    feSpace1DP4 U_M1D(mesh, "U", "M1D", funSol);
    feSpace1DP0 L_BXA(mesh, "L", "BXA", funLambda_A);
    feSpace1DP0 L_BXB(mesh, "L", "BXB", funLambda_B);

    std::vector<feSpace *> fespace = {&U_BXA, &U_BXB, &U_M1D, &L_BXA, &L_BXB};
    std::vector<feSpace *> feEssBC = {};
    // Numerotations
    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);
    // Solution
    double t0 = 0.;
    double t1 = 1.;
    int nTimeSteps = 40 * pow(2, iter);
    TT[iter] = nTimeSteps;
    feSolution *solDC3 = new feSolution(mesh, fespace, feEssBC, metaNumber);
    // solDC3->initializeTemporalSolution(t0, t1, nTimeSteps);
    // solDC3->initializeUnknowns(mesh, metaNumber);
    // solDC3->initializeEssentialBC(mesh, metaNumber);
    // Formes (bi)lineaires
    int degQuad = 11;
    std::vector<feSpace *> spaceDiffusion1D_U = {&U_M1D};
    std::vector<feSpace *> spaceSource1D_U = {&U_M1D};
    std::vector<feSpace *> spaceMasse1D_U = {&U_M1D};
    std::vector<feSpace *> spaceWeak_A = {&U_BXA, &L_BXA};
    std::vector<feSpace *> spaceWeak_B = {&U_BXB, &L_BXB};

    feBilinearForm *diff_U_M1D =
      new feBilinearForm(spaceDiffusion1D_U, mesh, degQuad, new feSysElm_1D_Diffusion(kd, nullptr));
    feBilinearForm *source_U_M1D =
      new feBilinearForm(spaceSource1D_U, mesh, degQuad, new feSysElm_1D_Source(1.0, funSource));
    feBilinearForm *masse_U_M1D =
      new feBilinearForm(spaceMasse1D_U, mesh, degQuad, new feSysElm_1D_Masse(1.0, nullptr));
    feBilinearForm *weakBC_U_A =
      new feBilinearForm(spaceWeak_A, mesh, degQuad, new feSysElm_0D_weakBC(1.0, funSol));
    feBilinearForm *weakBC_U_B =
      new feBilinearForm(spaceWeak_B, mesh, degQuad, new feSysElm_0D_weakBC(1.0, funSol));

    std::vector<feBilinearForm *> formMatrices = {diff_U_M1D, masse_U_M1D, weakBC_U_A, weakBC_U_B};
    std::vector<feBilinearForm *> formResiduals = {diff_U_M1D, masse_U_M1D, source_U_M1D,
                                                   weakBC_U_A, weakBC_U_B};
    // Norme de la solution
    feNorm *normBDF2 = new feNorm(&U_M1D, mesh, degQuad, funSol);
    feNorm *normDC3 = new feNorm(&U_M1D, mesh, degQuad, funSol);
    feNorm *normBDF2Lambda_A = new feNorm(&L_BXA, mesh, degQuad, funLambda_A);
    feNorm *normBDF2Lambda_B = new feNorm(&L_BXB, mesh, degQuad, funLambda_B);
    feNorm *normDC3Lambda_A = new feNorm(&L_BXA, mesh, degQuad, funLambda_A);
    feNorm *normDC3Lambda_B = new feNorm(&L_BXB, mesh, degQuad, funLambda_B);
    std::vector<feNorm *> norms = {normBDF2,         normDC3,         normBDF2Lambda_A,
                                   normBDF2Lambda_B, normDC3Lambda_A, normDC3Lambda_B};
    // Systeme lineaire
    feLinearSystemPETSc *linearSystem =
      new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_PETSC
    linearSystem->initialize();
    // Resolution
    feTolerances tol{1e-9, 1e-9, 20};
    std::string CodeIni = "BDF1/DC"; // Define the way of initialization |"SolEx"->for exact
                                     // solution|  |"BDF1/DCF"->using only initial conditions|
    DC3Solver solver(tol, metaNumber, linearSystem, solDC3, norms, mesh, t0, t1, nTimeSteps,
                     CodeIni);
    solver.makeSteps(nTimeSteps, fespace);
    std::vector<double> &normL2BDF2 = solver.getNorm(0);
    std::vector<double> &normL2DC3 = solver.getNorm(1);
    std::vector<double> &normL2_BDF2_Lambda_A = solver.getNorm(2);
    std::vector<double> &normL2_BDF2_Lambda_B = solver.getNorm(3);
    std::vector<double> &normL2_DC3_Lambda_A = solver.getNorm(4);
    std::vector<double> &normL2_DC3_Lambda_B = solver.getNorm(5);
    normL2_BDF2[2 * iter] = *std::max_element(normL2BDF2.begin(), normL2BDF2.end());
    normL2_DC3[2 * iter] = *std::max_element(normL2DC3.begin(), normL2DC3.end());
    normBDF2_A[2 * iter] =
      *std::max_element(normL2_BDF2_Lambda_A.begin(), normL2_BDF2_Lambda_A.end());
    normBDF2_B[2 * iter] =
      *std::max_element(normL2_BDF2_Lambda_B.begin(), normL2_BDF2_Lambda_B.end());
    normDC3_A[2 * iter] = *std::max_element(normL2_DC3_Lambda_A.begin(), normL2_DC3_Lambda_A.end());
    normDC3_B[2 * iter] = *std::max_element(normL2_DC3_Lambda_B.begin(), normL2_DC3_Lambda_B.end());
#endif

    delete normBDF2;
    delete normDC3;
    delete normBDF2Lambda_A;
    delete normBDF2Lambda_B;
    delete normDC3Lambda_A;
    delete normDC3Lambda_B;
    delete linearSystem;
    delete masse_U_M1D;
    delete source_U_M1D;
    delete diff_U_M1D;
    delete solDC3;
    delete metaNumber;
    delete mesh;
  }
  delete funSource;
  delete funSol;

  // Calcul du taux de convergence
  for(int i = 1; i < nIter; ++i) {
    normL2_BDF2[2 * i + 1] = log(normL2_BDF2[2 * (i - 1)] / normL2_BDF2[2 * i]) / log(2.);
    normL2_DC3[2 * i + 1] = log(normL2_DC3[2 * (i - 1)] / normL2_DC3[2 * i]) / log(2.);
    normBDF2_A[2 * i + 1] = log(normBDF2_A[2 * (i - 1)] / normBDF2_A[2 * i]) / log(2.);
    normBDF2_B[2 * i + 1] = log(normBDF2_A[2 * (i - 1)] / normBDF2_A[2 * i]) / log(2.);
    normDC3_A[2 * i + 1] = log(normDC3_A[2 * (i - 1)] / normDC3_A[2 * i]) / log(2.);
    normDC3_B[2 * i + 1] = log(normDC3_B[2 * (i - 1)] / normDC3_B[2 * i]) / log(2.);
  }
  printf("%12s \t %12s \t %12s \t %12s\t %12s \t %12s\n", "nSteps", "nElm", "||E-BDF||",
         "Taux BDF2", "||E-DC3||", "Taux DC3");
  for(int i = 0; i < nIter; ++i)

    printf("%12d \t %12d \t %12.6e \t %12.6e\t %12.6e \t %12.6e\n", TT[i], nElm[i],
           normL2_BDF2[2 * i], normL2_BDF2[2 * i + 1], normL2_DC3[2 * i], normL2_DC3[2 * i + 1]);
  printf("\n");
  printf("\n");
  printf("Lambda based on BDF2");
  printf("\n");
  printf("%12s \t %12s \t %16s \t %16s \t %16s \t %16s \n", "nTimeSteps", "nElm", "||E_lambda_A||",
         "Taux lambda_A", "||E_lambda_B||", "Taux lambda_B");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %16.6e \t %16.6e \t %16.6e \t %16.6e  \n", TT[i], nElm[i],
           normBDF2_A[2 * i], normBDF2_A[2 * i + 1], normBDF2_B[2 * i], normBDF2_B[2 * i + 1]);
  printf("\n");
  printf("\n");
  printf("Lambda based on DC3");
  printf("\n");
  printf("%12s \t %12s \t %16s \t %16s \t %16s \t %16s \n", "nTimeSteps", "nElm", "||E_lambda_A||",
         "Taux lambda_A", "||E_lambda_B||", "Taux lambda_B");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %16.6e \t %16.6e \t %16.6e \t %16.6e  \n", TT[i], nElm[i],
           normDC3_A[2 * i], normDC3_A[2 * i + 1], normDC3_B[2 * i], normDC3_B[2 * i + 1]);

#ifdef USING_PETSC
  petscFinalize();
#endif
  return 0;
}