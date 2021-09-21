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

#ifdef HAVE_PETSC
#include "feLinearSystemPETSc.h"
#endif

// double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   return sin(2 * M_PI * t) * cos(2 * M_PI * x1) + 1;
// }

// double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double c1 = par[0];
//   return -2 * M_PI * cos(2 * M_PI * t) * cos(2 * M_PI * x1) -
//          c1 * 2 * M_PI * 2 * M_PI * sin(2 * M_PI * t) * cos(2 * M_PI * x1);
// }

double f0(const double t, const std::vector<double> &x, const std::vector<double> par) {
  return 0.;
}

// double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   // double c2 = par[1];
//   double x1 = x[0];
//   return x1 * (x1 - 25.) * pow(t, 2) ;
// }

// double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double c1 = par[0];
//   // double c2 = par[1];
//   // double beta = (c2 == 0.) ? 0. : c2 * pow(t, c2 - 1.);
//   double x1 = x[0];
//   return -x1 * (x1 - 25.)*2.* pow(t, 1) + c1 * 2. * pow(t, 2.);
// }

double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  return pow(t, 3)*x1*x1 + 1;
}

double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double c1 = par[0];
  return -3.*pow(t, 2)*x1*x1 + c1*2*pow(t, 3);
}

double flambda_A(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double c1 = par[0];
  return  c1* 2*pow(t, 3)*x1 ; 
}
double flambda_B(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double c1 = par[0];
  return - c1* 2*pow(t, 3)*x1 ; 
}

int main(int argc, char **argv) {
#ifdef HAVE_PETSC
  petscInitialize(argc, argv);
#endif

  double xa = 0.;
  double xb = 5.;

  double kd = 1;
  std::vector<double> par = {kd, 3.};
  feFunction *funSol = new feFunction(fSol, par);
  feFunction *funSource = new feFunction(fSource, par);
  feFunction *fun0 = new feFunction(f0, {});
  feFunction *funLambda_A = new feFunction(flambda_A, par);
  feFunction *funLambda_B = new feFunction(flambda_B, par);

  int nIter = 6;
  std::vector<double> normL2(2 * nIter, 0.0);
  std::vector<double> normA(2 * nIter, 0.0);
  std::vector<double> normB(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);
  std::vector<int> TT;
  TT.resize(nIter);

  for(int iter = 0; iter < nIter; ++iter) {
    nElm[iter] = 5 * pow(2, iter);
    // Maillage
    feMesh1DP1 *mesh = new feMesh1DP1(xa, xb, nElm[iter], "BXA", "BXB", "M1D");
    // Espaces d'interpolation
    // feSpace1DP0 U_BXA = feSpace1DP0(mesh, "U", "BXA", funSol);
    // feSpace1DP0 U_BXB = feSpace1DP0(mesh, "U", "BXB", funSol);
    // feSpace1DP3 U_M1D = feSpace1DP3(mesh, "U", "M1D", fun0);
    feSpace1DP0 U_BXA(mesh, "U", "BXA", funSol);
    feSpace1DP0 U_BXB(mesh, "U", "BXB", funSol);
    feSpace1DP1 U_M1D(mesh, "U", "M1D", funSol);
    feSpace1DP0 L_BXA(mesh, "L", "BXA", funLambda_A);
    feSpace1DP0 L_BXB(mesh, "L", "BXB", funLambda_B);

    std::vector<feSpace *> fespace = {&U_BXA, &U_BXB, &U_M1D, &L_BXA, &L_BXB};
    std::vector<feSpace *> feEssBC = {};
    // Numerotations
    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);
    // Solution
    double t0 = 0.;
    double t1 = 1.;
    int nTimeSteps = 5 * pow(2, iter);
    TT[iter] = nTimeSteps;
    feSolution *sol = new feSolution(mesh, fespace, feEssBC, metaNumber);
    // sol->initializeTemporalSolution(t0,t1,nTimeSteps);
    // sol->initializeUnknowns(mesh, metaNumber);
    // sol->initializeEssentialBC(mesh, metaNumber);
    // Formes (bi)lineaires

    int nQuad = 6; // TODO : change to deg
    std::vector<feSpace *> spaceDiffusion1D_U = {&U_M1D};
    std::vector<feSpace *> spaceSource1D_U = {&U_M1D};
    std::vector<feSpace *> spaceMasse1D_U = {&U_M1D};
    std::vector<feSpace *> spaceWeak_A = {&U_BXA, &L_BXA};
    std::vector<feSpace *> spaceWeak_B = {&U_BXB, &L_BXB};

    feBilinearForm *diff_U_M1D =
      new feBilinearForm(spaceDiffusion1D_U, mesh, nQuad, new feSysElm_1D_Diffusion(kd, nullptr));
    feBilinearForm *source_U_M1D =
      new feBilinearForm(spaceSource1D_U, mesh, nQuad, new feSysElm_1D_Source(1.0, funSource));
    feBilinearForm *masse_U_M1D =
      new feBilinearForm(spaceMasse1D_U, mesh, nQuad, new feSysElm_1D_Masse(1.0, nullptr));
    feBilinearForm *weakBC_U_A =
      new feBilinearForm(spaceWeak_A, mesh, nQuad, new feSysElm_0D_weakBC(1.0, funSol));
    feBilinearForm *weakBC_U_B =
      new feBilinearForm(spaceWeak_B, mesh, nQuad, new feSysElm_0D_weakBC(1.0, funSol));

    std::vector<feBilinearForm *> formMatrices = {diff_U_M1D, masse_U_M1D,weakBC_U_A,weakBC_U_B};
    std::vector<feBilinearForm *> formResiduals = {diff_U_M1D, masse_U_M1D, source_U_M1D, weakBC_U_A, weakBC_U_B};

    // Norme de la solution
    feNorm *norm = new feNorm(&U_M1D, mesh, nQuad, funSol);
    feNorm *normLambda_A = new feNorm(&L_BXA, mesh, nQuad, funLambda_A);
    feNorm *normLambda_B = new feNorm(&L_BXB, mesh, nQuad, funLambda_B);
    std::vector<feNorm *> norms = {norm,normLambda_A,normLambda_B};
    // Systeme lineaire
    feLinearSystemPETSc *linearSystem =
      new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_PETSC
    // linearSystem->initialize();
    // Resolution
    feTolerances tol{1e-9, 1e-9, 20};
    // std::vector<double> normL2BDF(nTimeSteps,0.0);
    std::string CodeIni = "BDF1/DC"; // Define the way of initialization |"SolEx"->for exact
                                      // solution|  |"BDF1/DCF"->using only initial conditions|
    BDF2Solver solver(tol, metaNumber, linearSystem, sol, norms, mesh, t0, t1, nTimeSteps, CodeIni);
    solver.makeSteps(nTimeSteps, fespace);
    std::vector<double> &normL2BDF = solver.getNorm(0);
    std::vector<double> &normL2Lambda_A = solver.getNorm(1);
    std::vector<double> &normL2Lambda_B = solver.getNorm(2);

    // solveBDF2(normL2BDF, tol, metaNumber, linearSystem, formMatrices, formResiduals, sol, norms,
    // mesh, fespace);
    normL2[2 * iter] = *std::max_element(normL2BDF.begin(), normL2BDF.end());
    normA[2 * iter] = *std::max_element(normL2Lambda_A.begin(), normL2Lambda_A.end());
    normB[2 * iter] = *std::max_element(normL2Lambda_B.begin(), normL2Lambda_B.end());

#endif

    delete norm;
    delete linearSystem;
    delete masse_U_M1D;
    delete source_U_M1D;
    delete diff_U_M1D;
    delete sol;
    delete metaNumber;
    delete mesh;
  }
  delete funSource;
  delete funSol;

  // Calcul du taux de convergence
  for(int i = 1; i < nIter; ++i) {
    normL2[2 * i + 1] = log(normL2[2 * (i - 1)] / normL2[2 * i]) / log(2.);
    normA[2 * i + 1] = log(normA[2 * (i - 1)] / normA[2 * i]) / log(2.);
    normB[2 * i + 1] = log(normB[2 * (i - 1)] / normB[2 * i]) / log(2.);
  }
  printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s\n", "nTimeSteps", "nElm", "||E||", "Taux BDF2","||E_lambda_A||", "Taux lambda_A", "||E_lambda_B||", "Taux lambda_B" );
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n", TT[i], nElm[i], normL2[2 * i], normL2[2 * i + 1],normA[2 * i],normA[2 * i+1],normB[2 * i],normB[2 * i+1] );

#ifdef HAVE_PETSC
  petscFinalize();
#endif

  return 0;
}
