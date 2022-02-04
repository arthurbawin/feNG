#include <iostream>
#include <functional>
#include <vector>
#include <cmath>

#include "feNG.h"
#include "feFunction.h"
#include "feMesh.h"
#include "feQuadrature.h"
#include "feSpace.h"
#include "feCncGeo.h"
#include "feNumber.h"
#include "feSolution.h"
#include "feSysElm.h"
#include "feBilinearForm.h"
#include "feSolver.h"
#include "feLinearSystemPETSc.h"

double fSol(const double t, const std::vector<double> x, const std::vector<double> par)
{
  return exp(x[0]) + 1;
}

double f0(const double t, const std::vector<double> x, const std::vector<double> par)
{
  return 0.0;
}

double fSource(const double t, const std::vector<double> x, const std::vector<double> par)
{
  double kd = par[0];
  return exp(x[0]);
  // return kd * 2. * pow(x[0], 0)*0.;
}

int main(int argc, char **argv)
{
#ifdef HAVE_PETSC
  petscInitialize(argc, argv);
#endif

  double xa = 0.;
  double xb = 5.;

  double kd = 1.0;

  int nIter = 5;
  std::vector<double> normL2(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);

  feFunction *funSol = new feFunction(fSol, {kd});
  feFunction *fun0 = new feFunction(f0, {kd});
  feFunction *funSource = new feFunction(fSource, {kd});

  for(int iter = 0; iter < nIter; ++iter) {
    nElm[iter] = 20 * pow(2, iter);
    // nElm[iter] = 10;
    // Maillage
    feMesh1DP1 *mesh = new feMesh1DP1(xa, xb, nElm[iter], "BXA", "BXB", "M1D");
    // Espaces d'interpolation
    feSpace1DP0 U_BXA(mesh, "U", "BXA", funSol);
    feSpace1DP0 U_BXB(mesh, "U", "BXB", funSol);
    feSpace1DP3 U_M1D(mesh, "U", "M1D", fun0);
    std::vector<feSpace *> fespace = {&U_BXA, &U_BXB, &U_M1D};
    std::vector<feSpace *> feEssBC = {&U_BXA, &U_BXB};
    // Numerotations
    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);
    // metaNumber->printNumberings();
    // Solution
    feSolution *sol = new feSolution(mesh, fespace, feEssBC, metaNumber);
    // sol->initializeUnknowns(mesh, metaNumber);
    // sol->initializeEssentialBC(mesh, metaNumber);
    // Formes (bi)lineaires
    int degQuad = 9;
    std::vector<feSpace *> spaceDiffusion1D_U = {&U_M1D};
    feBilinearForm *diff_U_M1D =
      new feBilinearForm(spaceDiffusion1D_U, mesh, degQuad, new feSysElm_1D_Diffusion(kd, nullptr));
    std::vector<feSpace *> spaceSource1D_U = {&U_M1D};
    feBilinearForm *source_U_M1D =
      new feBilinearForm(spaceSource1D_U, mesh, degQuad, new feSysElm_1D_Source(1.0, funSource));

    std::vector<feBilinearForm *> formMatrices = {diff_U_M1D};
    std::vector<feBilinearForm *> formResiduals = {diff_U_M1D, source_U_M1D};

    feNorm *norm = new feNorm(&U_M1D, mesh, degQuad, funSol);
    std::vector<feNorm *> norms = {norm};

    feLinearSystemPETSc *linearSystem =
      new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
    linearSystem->initialize();
    // feTolerance {tolDx,tolResidual, maxIter}
    feTolerances tol{1e-5, 1e-5, 20};
    StationarySolver solver(tol, metaNumber, linearSystem, sol, norms, mesh);
    solver.makeSteps(0, fespace);
    normL2[2 * iter] = solver.getNorm(0);

    delete norm;
    delete linearSystem;
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
  }
  printf("%12s \t %12s \t %12s\n", "nElm", "||E||", "p");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12.6e \t %12.6e\n", nElm[i], normL2[2 * i], normL2[2 * i + 1]);
}
