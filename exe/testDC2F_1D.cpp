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
//   // double x1 = x[0];
//   return sin(2*M_PI*t) +1;
// }

// double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   // double x1 = x[0];
//   // double c1 = par[0];
//   return -2*M_PI*cos(2*M_PI*t) ;
// }

double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  return sin(2 * M_PI * t) * cos(2 * M_PI * x1) + 1;
}

double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double c1 = par[0];
  return -2 * M_PI * cos(2 * M_PI * t) * cos(2 * M_PI * x1) -
         c1 * 2 * M_PI * 2 * M_PI * sin(2 * M_PI * t) * cos(2 * M_PI * x1);
}

int main(int argc, char **argv) {
#ifdef USING_PETSC
  petscInitialize(argc, argv);
#endif
  double xa = 0.;
  double xb = 5.;

  double kd = 0.1;
  std::vector<double> par = {kd, 6.};
  feFunction *funSol = new feFunction(fSol, par);
  feFunction *funSource = new feFunction(fSource, par);

  int nIter = 3;
  std::vector<double> normL2_BDF1(2 * nIter, 0.0);
  std::vector<double> normL2_DC2F(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);
  std::vector<int> TT;
  TT.resize(nIter);
  for(int iter = 0; iter < nIter; ++iter) {
    nElm[iter] = 40 * pow(3, iter);
    // Maillage
    feMesh1DP1 *mesh = new feMesh1DP1(xa, xb, nElm[iter], "BXA", "BXB", "M1D");
    // Espaces d'interpolation
    // feSpace1DP0 U_BXA = feSpace1DP0(mesh, "U", "BXA", funSol);
    // feSpace1DP0 U_BXB = feSpace1DP0(mesh, "U", "BXB", funSol);
    // feSpace1DP1 U_M1D = feSpace1DP1(mesh, "U", "M1D", funSol);
    feSpace1DP0 U_BXA(mesh, "U", "BXA", funSol);
    feSpace1DP0 U_BXB(mesh, "U", "BXB", funSol);
    feSpace1DP2 U_M1D(mesh, "U", "M1D", funSol);
    std::vector<feSpace *> fespace = {&U_BXA, &U_BXB, &U_M1D};
    std::vector<feSpace *> feEssBC = {&U_BXA, &U_BXB};
    // Numerotations
    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);
    // Solution
    double t0 = 0.;
    double t1 = 1.;
    int nTimeSteps = 40 * pow(2, iter);
    TT[iter] = nTimeSteps;
    feSolution *solDC2 = new feSolution(mesh, fespace, feEssBC, metaNumber);
    // sol->initializeTemporalSolution(t0,t1,nTimeSteps);
    // sol->initializeUnknowns(fmesh, metaNumber);
    // sol->initializeEssentialBC(mesh, metaNumber);
    // Formes (bi)lineaires
    int nQuad = 6; // TODO : change to deg
    std::vector<feSpace *> spaceDiffusion1D_U = {&U_M1D};
    std::vector<feSpace *> spaceSource1D_U = {&U_M1D};
    std::vector<feSpace *> spaceMasse1D_U = {&U_M1D};

    feBilinearForm *diff_U_M1D =
      new feBilinearForm(spaceDiffusion1D_U, mesh, nQuad, new feSysElm_1D_Diffusion(kd, nullptr));
    feBilinearForm *source_U_M1D =
      new feBilinearForm(spaceSource1D_U, mesh, nQuad, new feSysElm_1D_Source(1.0, funSource));
    feBilinearForm *masse_U_M1D =
      new feBilinearForm(spaceMasse1D_U, mesh, nQuad, new feSysElm_1D_Masse(1.0, nullptr));

    std::vector<feBilinearForm *> formMatrices = {diff_U_M1D, masse_U_M1D};
    std::vector<feBilinearForm *> formResiduals = {diff_U_M1D, masse_U_M1D, source_U_M1D};
    // std::vector<feBilinearForm*> formMatrices  = {masse_U_M1D};
    // std::vector<feBilinearForm*> formResiduals  = {masse_U_M1D, source_U_M1D};
    // Norme de la solution
    feNorm *normBDF1 = new feNorm(&U_M1D, mesh, nQuad, funSol);
    feNorm *normDC2F = new feNorm(&U_M1D, mesh, nQuad, funSol);
    std::vector<feNorm *> norms = {normBDF1, normDC2F};
    // Systeme lineaire
    feLinearSystemPETSc *linearSystem =
      new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_PETSC
    linearSystem->initialize();
    // Resolution
    feTolerances tol{1e-10, 1e-10, 10};

    DC2FSolver solver(tol, metaNumber, linearSystem, solDC2, norms, mesh, t0, t1, nTimeSteps);
    solver.makeSteps(nTimeSteps, fespace);
    std::vector<double> &normL2BDF1 = solver.getNorm(0);
    std::vector<double> &normL2DC2F = solver.getNorm(1);

    normL2_BDF1[2 * iter] = *std::max_element(normL2BDF1.begin(), normL2BDF1.end());
    normL2_DC2F[2 * iter] = *std::max_element(normL2DC2F.begin(), normL2DC2F.end());
#endif

    delete normBDF1;
    delete normDC2F;
    delete linearSystem;
    delete masse_U_M1D;
    delete source_U_M1D;
    delete diff_U_M1D;
    delete solDC2;
    delete metaNumber;
    delete mesh;
  }
  delete funSource;
  delete funSol;

  // Calcul du taux de convergence
  for(int i = 1; i < nIter; ++i) {
    normL2_BDF1[2 * i + 1] = log(normL2_BDF1[2 * (i - 1)] / normL2_BDF1[2 * i]) / log(2.);
    normL2_DC2F[2 * i + 1] = log(normL2_DC2F[2 * (i - 1)] / normL2_DC2F[2 * i]) / log(2.);
  }
  printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \n", "nTimeSteps", "nElm", "||E_BDF1||",
         "Taux BDF1", "||E_DC2F||", "Taux DC2F");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e \n", TT[i], nElm[i],
           normL2_BDF1[2 * i], normL2_BDF1[2 * i + 1], normL2_DC2F[2 * i], normL2_DC2F[2 * i + 1]);
#ifdef USING_PETSC
  petscFinalize();
#endif
  return 0;
}
