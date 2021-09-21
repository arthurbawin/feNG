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
// #include "feGetSolFromFile.h"
#define USING_PETSC
#ifdef HAVE_PETSC
#include "feLinearSystemPETSc.h"
#endif

// double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double c2 = par[1];
//   double x1 = x[0];
//   return x1 * (x1 * x1 - 25.) * pow(t, 2) +1;
// }

// double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double c1 = par[0];
//   double c2 = par[1];
//   double beta = (c2 == 0.) ? 0. : c2 * pow(t, c2 - 1.);
//   double x1 = x[0];
//   return -x1 * (x1 * x1 - 25.) *  2.*pow(t, 1.)  + c1 * 6. * x1 * pow(t, 2.);
// }

// double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   return sin(2*M_PI*t)*cos(2*M_PI*x1) +1;
// }

// double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double c1 = par[0];
//   return -2*M_PI*cos(2*M_PI*t)*cos(2*M_PI*x1) - c1*2*M_PI*2*M_PI*sin(2*M_PI*t)*cos(2*M_PI*x1);
// }

double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  return pow(t, 3)*x1*(x1-5) + 1;
}

double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double c1 = par[0];
  return -3.*pow(t, 2)*x1*(x1-5) + c1*2*pow(t, 3);
}

double flambda_A(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double c1 = par[0];
  return  c1* pow(t, 3)*(2*x1 - 5.) ; 
}
double flambda_B(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double c1 = par[0];
  return - c1* pow(t, 3)*(2*x1 -5.); 
}

double fSolDot(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  return 3*pow(t, 2)*x1*(x1-5) ;
}

double fSolDotDot(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  return 6*pow(t, 1)*x1*(x1-5) ;
}

int main(int argc, char **argv) {
#ifdef USING_PETSC
  petscInitialize(argc, argv);
#endif
  double xa = 0.;
  double xb = 5.;

  double kd = 1;
  std::vector<double> par = {kd, 6.};
  feFunction *funSol = new feFunction(fSol, par);
  feFunction *funSolDot = new feFunction(fSolDot, par);
  feFunction *funSolDotDot = new feFunction(fSolDotDot, par);
  feFunction *funSource = new feFunction(fSource, par);
  feFunction *funLambda_A = new feFunction(flambda_A, par);
  feFunction *funLambda_B = new feFunction(flambda_B, par);


  int nIter = 4;
  std::vector<double> normL2_BDF1(2 * nIter, 0.0);
  std::vector<double> normL2_DC2F(2 * nIter, 0.0);
  std::vector<double> normL2_DC3F(2 * nIter, 0.0);
  std::vector<double> normBDF1_A(2 * nIter, 0.0);
  std::vector<double> normBDF1_B(2 * nIter, 0.0);
  std::vector<double> normDC2F_A(2 * nIter, 0.0);
  std::vector<double> normDC2F_B(2 * nIter, 0.0);
  std::vector<double> normDC3F_A(2 * nIter, 0.0);
  std::vector<double> normDC3F_B(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);
  std::vector<int> TT;
  TT.resize(nIter);
  for(int iter = 0; iter < nIter; ++iter) {
    nElm[iter] = 20 * pow(3, iter);
    // Maillage
    feMesh1DP1 *mesh = new feMesh1DP1(xa, xb, nElm[iter], "BXA", "BXB", "M1D");
    // Espaces d'interpolation
    // feSpace1DP0 U_BXA = feSpace1DP0(mesh, "U", "BXA", funSol);
    // feSpace1DP0 U_BXB = feSpace1DP0(mesh, "U", "BXB", funSol);
    // feSpace1DP3 U_M1D = feSpace1DP3(mesh, "U", "M1D", funSol);
    feSpace1DP0 U_BXA(mesh, "U", "BXA", funSol);
    feSpace1DP0 U_BXB(mesh, "U", "BXB", funSol);
    feSpace1DP3 U_M1D(mesh, "U", "M1D", funSol);
    feSpace1DP0 L_BXA(mesh, "L", "BXA", funLambda_A);
    feSpace1DP0 L_BXB(mesh, "L", "BXB", funLambda_B);
    feSpace1DP0 V_BXA(mesh, "V", "BXA", funSol);
    feSpace1DP0 V_BXB(mesh, "V", "BXB", funSol);
    feSpace1DP0 W_BXA(mesh, "W", "BXA", funSolDot);
    feSpace1DP0 W_BXB(mesh, "W", "BXB", funSolDot);
    std::vector<feSpace *> fespace = {&U_BXA, &U_BXB, &U_M1D, &L_BXA, &L_BXB,&V_BXA,&V_BXB, &W_BXA, &W_BXB};
    std::vector<feSpace *> feEssBC = {};
    // Numerotations
    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);
    // Solution
    double t0 = 0.;
    double t1 = 1.;
    int nTimeSteps = 20 * pow(2, iter);
    TT[iter] = nTimeSteps;
    feSolution *solDC3F = new feSolution(mesh, fespace, feEssBC, metaNumber);
    // sol->initializeTemporalSolution(t0,t1,nTimeSteps);
    // sol->initializeUnknowns(fmesh, metaNumber);
    // sol->initializeEssentialBC(mesh, metaNumber);
    // Formes (bi)lineaires
    int nQuad = 9; // TODO : change to deg
    std::vector<feSpace *> spaceDiffusion1D_U = {&U_M1D};
    std::vector<feSpace *> spaceSource1D_U = {&U_M1D};
    std::vector<feSpace *> spaceMasse1D_U = {&U_M1D};
    std::vector<feSpace *> spaceWeak_A = {&U_BXA, &L_BXA, &V_BXA,&W_BXA};
    std::vector<feSpace *> spaceWeak_B = {&U_BXB, &L_BXB, &V_BXB,&W_BXB};

    feBilinearForm *diff_U_M1D =
      new feBilinearForm(spaceDiffusion1D_U, mesh, nQuad, new feSysElm_1D_Diffusion(kd, nullptr));
    feBilinearForm *source_U_M1D =
      new feBilinearForm(spaceSource1D_U, mesh, nQuad, new feSysElm_1D_Source(1.0, funSource));
    feBilinearForm *masse_U_M1D =
      new feBilinearForm(spaceMasse1D_U, mesh, nQuad, new feSysElm_1D_Masse(1.0, nullptr));
    feBilinearForm *weakBC_U_A =
      new feBilinearForm(spaceWeak_A, mesh, nQuad, new feSysElm_0D_weakBC_edo2(1.0, funSolDotDot));
    feBilinearForm *weakBC_U_B =
      new feBilinearForm(spaceWeak_B, mesh, nQuad, new feSysElm_0D_weakBC_edo2(1.0, funSolDotDot));

    std::vector<feBilinearForm *> formMatrices = {diff_U_M1D, masse_U_M1D,weakBC_U_A,weakBC_U_B};
    std::vector<feBilinearForm *> formResiduals = {diff_U_M1D, masse_U_M1D, source_U_M1D, weakBC_U_A, weakBC_U_B};
    // std::vector<feBilinearForm*> formMatrices  = {masse_U_M1D};
    // std::vector<feBilinearForm*> formResiduals  = {masse_U_M1D, source_U_M1D};
    // Norme de la solution
    feNorm *normBDF1 = new feNorm(&U_M1D, mesh, nQuad, funSol);
    feNorm *normDC2F = new feNorm(&U_M1D, mesh, nQuad, funSol);
    feNorm *normDC3F = new feNorm(&U_M1D, mesh, nQuad, funSol);
    feNorm *normBDF1Lambda_A = new feNorm(&L_BXA, mesh, nQuad, funLambda_A);
    feNorm *normBDF1Lambda_B = new feNorm(&L_BXB, mesh, nQuad, funLambda_B);
    feNorm *normDC2FLambda_A = new feNorm(&L_BXA, mesh, nQuad, funLambda_A);
    feNorm *normDC2FLambda_B = new feNorm(&L_BXB, mesh, nQuad, funLambda_B);
    feNorm *normDC3FLambda_A = new feNorm(&L_BXA, mesh, nQuad, funLambda_A);
    feNorm *normDC3FLambda_B = new feNorm(&L_BXB, mesh, nQuad, funLambda_B);
    std::vector<feNorm *> norms = {normBDF1, normDC2F, normDC3F, normBDF1Lambda_A,normBDF1Lambda_B, normDC2FLambda_A, normDC2FLambda_B, normDC3FLambda_A, normDC3FLambda_B};
    // Systeme lineaire
    feLinearSystemPETSc *linearSystem =
      new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_PETSC
    // linearSystem->initialize();
    // Resolution
    feTolerances tol{1e-9, 1e-9, 20};

    DC3FSolver solver(tol, metaNumber, linearSystem, solDC3F, norms, mesh, t0, t1, nTimeSteps);
    solver.makeSteps(nTimeSteps, fespace);
    std::vector<double> &normL2BDF1 = solver.getNorm(0);
    std::vector<double> &normL2DC2F = solver.getNorm(1);
    std::vector<double> &normL2DC3F = solver.getNorm(2);
    std::vector<double> &normL2_BDF1_Lambda_A = solver.getNorm(3);
    std::vector<double> &normL2_BDF1_Lambda_B = solver.getNorm(4);
    std::vector<double> &normL2_DC2F_Lambda_A = solver.getNorm(5);
    std::vector<double> &normL2_DC2F_Lambda_B = solver.getNorm(6);
    std::vector<double> &normL2_DC3F_Lambda_A = solver.getNorm(7);
    std::vector<double> &normL2_DC3F_Lambda_B = solver.getNorm(8);

    normL2_BDF1[2 * iter] = *std::max_element(normL2BDF1.begin(), normL2BDF1.end());
    normL2_DC2F[2 * iter] = *std::max_element(normL2DC2F.begin(), normL2DC2F.end());
    normL2_DC3F[2 * iter] = *std::max_element(normL2DC3F.begin(), normL2DC3F.end());
    normBDF1_A[2 * iter]  = *std::max_element(normL2_BDF1_Lambda_A.begin(), normL2_BDF1_Lambda_A.end());
    normBDF1_B[2 * iter]  = *std::max_element(normL2_BDF1_Lambda_B.begin(), normL2_BDF1_Lambda_B.end());
    normDC2F_A[2 * iter]  = *std::max_element(normL2_DC2F_Lambda_A.begin(), normL2_DC2F_Lambda_A.end());
    normDC2F_B[2 * iter]  = *std::max_element(normL2_DC2F_Lambda_B.begin(), normL2_DC2F_Lambda_B.end());
    normDC3F_A[2 * iter]  = *std::max_element(normL2_DC3F_Lambda_A.begin(), normL2_DC3F_Lambda_A.end());
    normDC3F_B[2 * iter]  = *std::max_element(normL2_DC3F_Lambda_B.begin(), normL2_DC3F_Lambda_B.end());
#endif

    delete normBDF1;
    delete normDC2F;
    delete normDC3F;
    delete normBDF1Lambda_A;
    delete normBDF1Lambda_B;
    delete normDC2FLambda_A;
    delete normDC2FLambda_B;
    delete normDC3FLambda_A;
    delete normDC3FLambda_B;
    delete linearSystem;
    delete masse_U_M1D;
    delete source_U_M1D;
    delete diff_U_M1D;
    delete solDC3F;
    delete metaNumber;
    delete mesh;
  }
  delete funSource;
  delete funSol;

  // Calcul du taux de convergence
  for(int i = 1; i < nIter; ++i) {
    normL2_BDF1[2 * i + 1] = log(normL2_BDF1[2 * (i - 1)] / normL2_BDF1[2 * i]) / log(2.);
    normL2_DC2F[2 * i + 1] = log(normL2_DC2F[2 * (i - 1)] / normL2_DC2F[2 * i]) / log(2.);
    normL2_DC3F[2 * i + 1] = log(normL2_DC3F[2 * (i - 1)] / normL2_DC3F[2 * i]) / log(2.);
    normBDF1_A[2 * i + 1] = log(normBDF1_A[2 * (i - 1)] / normBDF1_A[2 * i]) / log(2.);
    normBDF1_B[2 * i + 1] = log(normBDF1_B[2 * (i - 1)] / normBDF1_B[2 * i]) / log(2.);
    normDC2F_A[2 * i + 1] = log(normDC2F_A[2 * (i - 1)] / normDC2F_A[2 * i]) / log(2.);
    normDC2F_B[2 * i + 1] = log(normDC2F_B[2 * (i - 1)] / normDC2F_B[2 * i]) / log(2.);
    normDC3F_A[2 * i + 1] = log(normDC3F_A[2 * (i - 1)] / normDC3F_A[2 * i]) / log(2.);
    normDC3F_B[2 * i + 1] = log(normDC3F_B[2 * (i - 1)] / normDC3F_B[2 * i]) / log(2.);
  }
  printf("\n");
  printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \n", "nTimeSteps", "nElm",
         "||E_BDF1||", "Taux BDF1", "||E_DC2F||", "Taux DC2F", "||E_DC3F||", "Taux DC3F");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \n", TT[i],
           nElm[i], normL2_BDF1[2 * i], normL2_BDF1[2 * i + 1], normL2_DC2F[2 * i],
           normL2_DC2F[2 * i + 1], normL2_DC3F[2 * i], normL2_DC3F[2 * i + 1]);
  printf("\n");
  printf("\n");
  printf("Lambda based on BDF1");
  printf("\n");
  printf("%12s \t %12s \t %16s \t %16s \t %16s \t %16s \n", "nTimeSteps", "nElm",
         "||E_lambda_A||", "Taux lambda_A", "||E_lambda_B||", "Taux lambda_B");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %16.6e \t %16.6e \t %16.6e \t %16.6e  \n", TT[i],
           nElm[i], normBDF1_A[2 * i], normBDF1_A[2 * i + 1], normBDF1_B[2 * i],
           normBDF1_B[2 * i + 1]);
  printf("\n");
  printf("\n");
  printf("Lambda based on DC2F");
  printf("\n");
  printf("%12s \t %12s \t %16s \t %16s \t %16s \t %16s \n", "nTimeSteps", "nElm",
         "||E_lambda_A||", "Taux lambda_A", "||E_lambda_B||", "Taux lambda_B");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %16.6e \t %16.6e \t %16.6e \t %16.6e  \n", TT[i],
           nElm[i], normDC2F_A[2 * i], normDC2F_A[2 * i + 1], normDC2F_B[2 * i],
           normDC2F_B[2 * i + 1]);
  printf("\n");
  printf("\n");
  printf("Lambda based on DC3F");
  printf("\n");
  printf("%12s \t %12s \t %16s \t %16s \t %16s \t %16s \n", "nTimeSteps", "nElm",
         "||E_lambda_A||", "Taux lambda_A", "||E_lambda_B||", "Taux lambda_B");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %16.6e \t %16.6e \t %16.6e \t %16.6e  \n", TT[i],
           nElm[i], normDC3F_A[2 * i], normDC3F_A[2 * i + 1], normDC3F_B[2 * i],
           normDC3F_B[2 * i + 1]);
#ifdef USING_PETSC
  petscFinalize();
#endif
  return 0;
}
