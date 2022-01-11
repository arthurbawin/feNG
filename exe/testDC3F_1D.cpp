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
//   return sin(4*t)*exp(2*x1) +1;
// }

// double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double c1 = par[0];
//   return -4*cos(4*t)*exp(2*x1) + c1*4*sin(4*t)*exp(2*x1);
// }

double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double a = par[1];
  return x1 * x1 * pow(t, a);
}

double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double c1 = par[0];
  double a = par[1];
  double beta = (a == 0.) ? 0. : a * pow(t, a - 1.);
  return -beta * x1 * x1 + c1 * 2 * pow(t, a);
  // return -beta*x1;
}

int main(int argc, char **argv) {
#ifdef USING_PETSC
  petscInitialize(argc, argv);
#endif
  double xa = 0.;
  double xb = 5.;
  double exposant = 2.;
  double kd = 0.1;
  std::vector<double> par = {kd, exposant};
  feFunction *funSol = new feFunction(fSol, par);
  feFunction *funSource = new feFunction(fSource, par);

  int nIter = 4;
  std::vector<double> normL2_BDF1(2 * nIter, 0.0);
  std::vector<double> normL2_DC2F(2 * nIter, 0.0);
  std::vector<double> normL2_DC3F(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);
  std::vector<int> TT;
  TT.resize(nIter);
  for(int iter = 0; iter < nIter; ++iter) {
    // nElm[iter] = 40 * pow(3, iter);
    nElm[iter] = 2;
    // Maillage
    feMesh1DP1 *mesh = new feMesh1DP1(xa, xb, nElm[iter], "BXA", "BXB", "M1D");
    // Espaces d'interpolation
    feSpace1DP0 U_BXA = feSpace1DP0(mesh, "U", "BXA", funSol);
    feSpace1DP0 U_BXB = feSpace1DP0(mesh, "U", "BXB", funSol);
    feSpace1DP2 U_M1D = feSpace1DP2(mesh, "U", "M1D", funSol);
    // feSpace1DP0 U_BXA(mesh, "U", "BXA", funSol);
    // feSpace1DP0 U_BXB(mesh, "U", "BXB", funSol);
    // feSpace1DP3 U_M1D(mesh, "U", "M1D", funSol);
    std::vector<feSpace *> fespace = {&U_BXA, &U_BXB, &U_M1D};
    std::vector<feSpace *> feEssBC = {&U_BXA, &U_BXB};
    // Numerotations
    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);
    // Solution
    double t0 = 0.;
    double t1 = 1.;
    int nTimeSteps = 5 * pow(2, iter);
    TT[iter] = nTimeSteps;
    feSolution *solDC3F = new feSolution(mesh, fespace, feEssBC, metaNumber);
    // sol->initializeTemporalSolution(t0,t1,nTimeSteps);
    // sol->initializeUnknowns(fmesh, metaNumber);
    // sol->initializeEssentialBC(mesh, metaNumber);
    // Formes (bi)lineaires
    int nQuad = 3; // TODO : change to deg
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
    feNorm *normDC3F = new feNorm(&U_M1D, mesh, nQuad, funSol);
    std::vector<feNorm *> norms = {normBDF1, normDC2F, normDC3F};
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
    std::cout << "size of normL2BDF1 " << normL2BDF1.size() << "size of normL2DC2F "
              << normL2DC2F.size() << "size of normL2DC3F " << normL2DC3F.size() << std::endl;
    for(int k = 0; k < normL2DC3F.size(); k++) std::cout << normL2DC3F[k] << std::endl;

    normL2_BDF1[2 * iter] = *std::max_element(normL2BDF1.begin(), normL2BDF1.end());
    normL2_DC2F[2 * iter] = *std::max_element(normL2DC2F.begin(), normL2DC2F.end());
    normL2_DC3F[2 * iter] = *std::max_element(normL2DC3F.begin(), normL2DC3F.end());
#endif

    delete normBDF1;
    delete normDC2F;
    delete normDC3F;
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
  }
  printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \n", "nTimeSteps", "nElm",
         "||E_BDF1||", "Taux BDF1", "||E_DC2F||", "Taux DC2F", "||E_DC3F||", "Taux DC3F");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \n", TT[i],
           nElm[i], normL2_BDF1[2 * i], normL2_BDF1[2 * i + 1], normL2_DC2F[2 * i],
           normL2_DC2F[2 * i + 1], normL2_DC3F[2 * i], normL2_DC3F[2 * i + 1]);
#ifdef USING_PETSC
  petscFinalize();
#endif
  return 0;
}
