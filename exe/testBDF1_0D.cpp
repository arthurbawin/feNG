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

double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
  return pow(t, 2);
}

double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
  return -2 * pow(t, 1);
}

int main(int argc, char **argv) {
#ifdef USING_PETSC
  petscInitialize(argc, argv);
#endif
  double xa = 0.;

  feFunction *funSol = new feFunction(fSol, {});
  feFunction *funSource = new feFunction(fSource, {});

  int nIter = 5;
  std::vector<double> normL2BDF1(2 * nIter, 0.0);
  // std::vector<double> normL2BDF2(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);
  std::vector<int> numberTimeStep;
  numberTimeStep.resize(nIter);
  for(int iter = 0; iter < nIter; ++iter) {
    nElm[iter] = 1;
    // Maillage
    feMesh0DP0 *mesh = new feMesh0DP0(xa, 1, "M1D");
    // Espaces d'interpolation
    feSpace1DP0 U_M1D = feSpace1DP0(mesh, "U", "M1D", funSol);
    std::vector<feSpace *> fespace = {&U_M1D};
    std::vector<feSpace *> feEssBC = {};
    // Numerotations
    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);
    // Solution
    double t0 = 0.;
    double t1 = 1.;
    int nTimeSteps = 5 * pow(2, iter);
    numberTimeStep[iter] = nTimeSteps;
    feSolution *sol = new feSolution(mesh, fespace, feEssBC, metaNumber);
    // sol->initializeTemporalSolution(t0,t1,nTimeSteps);
    // sol->initializeUnknowns(mesh, metaNumber);
    // sol->initializeEssentialBC(mesh, metaNumber);
    // Formes (bi)lineaires
    int nQuad = 8; // TODO : change to deg
    // std::vector<feSpace*> spaceDiffusion1D_U = {&U_M1D};
    std::vector<feSpace *> spaceSource1D_U = {&U_M1D};
    std::vector<feSpace *> spaceMasse1D_U = {&U_M1D};

    feBilinearForm *source_U_M1D =
      new feBilinearForm(spaceSource1D_U, mesh, nQuad, new feSysElm_0D_Source(1.0, funSource));
    feBilinearForm *masse_U_M1D =
      new feBilinearForm(spaceMasse1D_U, mesh, nQuad, new feSysElm_0D_Masse(1.0, nullptr));

    std::vector<feBilinearForm *> formMatrices = {masse_U_M1D};
    std::vector<feBilinearForm *> formResiduals = {masse_U_M1D, source_U_M1D};
    // Norme de la solution
    feNorm *norm = new feNorm(&U_M1D, mesh, nQuad, funSol);
    std::vector<feNorm *> norms = {norm};
    // Systeme lineaire
    feLinearSystemPETSc *linearSystem =
      new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_PETSC
    linearSystem->initialize();
    // Resolution
    feTolerances tol{1e-10, 1e-12, 10};
    // std::vector<double> normL2BDF1(nTimeSteps,0.0);

    BDF1Solver solver(tol, metaNumber, linearSystem, sol, norms, mesh, t0, t1, nTimeSteps);
    solver.makeSteps(nTimeSteps, fespace);
    std::vector<double> &normL2BDF = solver.getNorm(0);
    // solveBDF2afterBDF1(normL2BDF1,normL2BDF2, tol, metaNumber, linearSystem, formMatrices,
    // formResiduals, sol, norms, mesh, fespace);
    normL2BDF1[2 * iter] = *std::max_element(normL2BDF.begin(), normL2BDF.end());
#endif

    delete norm;
    delete linearSystem;
    delete masse_U_M1D;
    delete source_U_M1D;
    delete sol;
    delete metaNumber;
    delete mesh;
  }
  delete funSource;
  delete funSol;

  // Calcul du taux de convergence
  for(int i = 1; i < nIter; ++i) {
    normL2BDF1[2 * i + 1] = log(normL2BDF1[2 * (i - 1)] / normL2BDF1[2 * i]) / log(2.);
  }
  printf("%12s \t %12s \t %12s \t %12s \n", "nTimeSteps", "nElm", "||E||", "Taux BDF1");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %12.6e \t %12.6e\n", numberTimeStep[i], nElm[i], normL2BDF1[2 * i],
           normL2BDF1[2 * i + 1]);
#ifdef USING_PETSC
  petscFinalize();
#endif
  return 0;
}
