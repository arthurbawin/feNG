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

// double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double c2 = par[1];
//   double x1 = x[0];
//   return x1 * (x1 * x1 - 25.) * pow(t, 6) + x1 * x1 * x1 * pow(t, c2) +1.;
// }

// double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double c1 = par[0];
//   double c2 = par[1];
//   double beta = (c2 == 0.) ? 0. : c2 * pow(t, c2 - 1.);
//   double x1 = x[0];
//   return -x1 * (x1 * x1 - 25.) * 6. * pow(t, 5.) + c1 * 6. * x1 * pow(t, 6.) +
//          c1 * 6. * x1 * pow(t, c2) - x1 * x1 * x1 * beta;
// }

double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  return x1 * x1 * x1 * pow(t, 6.);
}

double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double c1 = par[0];
  return -6. * pow(t, 5.) * x1 * x1 * x1 + 6 * x1 * c1 * pow(t, 6.);
}

// double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double a = par[1];
//   return (x1*x1*+1)*pow(t, a);
// }

// double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double c1 = par[0];
//   double a = par[1];
//   double beta = (a == 0.) ? 0. : a * pow(t, a - 1.);
//   return -beta*(x1*x1*+1) + c1*2*pow(t, a);
// }

// double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double a = par[1];
//   return pow(t, a);
// }

// double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double c1 = par[0];
//   double a = par[1];
//   double beta = (a == 0.) ? 0. : a * pow(t, a - 1.);
//   return -beta;
// }

int main(int argc, char **argv) {
#ifdef USING_PETSC
  petscInitialize(argc, argv);
#endif
  double xa = 0.;
  double xb = 5.;
  double exposant = 6.;
  double kd = 0.1;
  std::vector<double> par = {kd, exposant};
  feFunction *funSol = new feFunction(fSol, par);
  feFunction *funSource = new feFunction(fSource, par);

  int nIter = 4;
  std::vector<double> normL2_BDF2(2 * nIter, 0.0);
  std::vector<double> normL2_DC3(2 * nIter, 0.0);

  std::vector<int> nElm(nIter, 0);
  std::vector<int> TT;
  TT.resize(nIter);

  for(int iter = 0; iter < nIter; ++iter) {
    nElm[iter] = 40 * pow(2, iter);
    // nElm[iter] = 2;
    // Maillage
    feMesh1DP1 *mesh = new feMesh1DP1(xa, xb, nElm[iter], "BXA", "BXB", "M1D");
    // Espaces d'interpolation
    feSpace1DP0 U_BXA = feSpace1DP0(mesh, "U", "BXA", funSol);
    feSpace1DP0 U_BXB = feSpace1DP0(mesh, "U", "BXB", funSol);
    feSpace1DP2 U_M1D = feSpace1DP2(mesh, "U", "M1D", funSol);
    std::vector<feSpace *> fespace = {&U_BXA, &U_BXB, &U_M1D};
    std::vector<feSpace *> feEssBC = {&U_BXA, &U_BXB};
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

    feBilinearForm *diff_U_M1D =
      new feBilinearForm(spaceDiffusion1D_U, mesh, degQuad, new feSysElm_1D_Diffusion(kd, nullptr));
    feBilinearForm *source_U_M1D =
      new feBilinearForm(spaceSource1D_U, mesh, degQuad, new feSysElm_1D_Source(1.0, funSource));
    feBilinearForm *masse_U_M1D =
      new feBilinearForm(spaceMasse1D_U, mesh, degQuad, new feSysElm_1D_Masse(1.0, nullptr));

    std::vector<feBilinearForm *> formMatrices = {diff_U_M1D, masse_U_M1D};
    std::vector<feBilinearForm *> formResiduals = {diff_U_M1D, masse_U_M1D, source_U_M1D};
    // Norme de la solution
    feNorm *normBDF2 = new feNorm(&U_M1D, mesh, degQuad, funSol);
    feNorm *normDC3 = new feNorm(&U_M1D, mesh, degQuad, funSol);
    std::vector<feNorm *> norms = {normBDF2, normDC3};
    // Systeme lineaire
    feLinearSystemPETSc *linearSystem =
      new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_PETSC
    linearSystem->initialize();
    // Resolution
    feTolerances tol{1e-9, 1e-9, 20};
    std::string CodeIni = "BDF1/DCF"; // Define the way of initialization |"SolEx"->for exact
                                      // solution|  |"BDF1/DCF"->using only initial conditions|
    DC3Solver solver(tol, metaNumber, linearSystem, solDC3, norms, mesh, t0, t1, nTimeSteps,
                     CodeIni);
    solver.makeSteps(nTimeSteps, fespace);
    std::vector<double> &normL2BDF2 = solver.getNorm(0);
    std::vector<double> &normL2DC3 = solver.getNorm(1);
    normL2_BDF2[2 * iter] = *std::max_element(normL2BDF2.begin(), normL2BDF2.end());
    normL2_DC3[2 * iter] = *std::max_element(normL2DC3.begin(), normL2DC3.end());
    for(int k = 0; k < 5; k++) std::cout << normL2DC3[k] << std::endl;
#endif

    delete normBDF2;
    delete normDC3;
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
  }
  printf("%12s \t %12s \t %12s \t %12s\t %12s \t %12s\n", "nSteps", "nElm", "||E-BDF||",
         "Taux BDF2", "||E-DC3||", "Taux DC3");
  for(int i = 0; i < nIter; ++i)

    printf("%12d \t %12d \t %12.6e \t %12.6e\t %12.6e \t %12.6e\n", TT[i], nElm[i],
           normL2_BDF2[2 * i], normL2_BDF2[2 * i + 1], normL2_DC3[2 * i], normL2_DC3[2 * i + 1]);
#ifdef USING_PETSC
  petscFinalize();
#endif
  return 0;
}
