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

double fSol(const double t, const std::vector<double> &x, const std::vector<double> par){
  return pow(t,7);
}

double fSource(const double t, const std::vector<double> &x, const std::vector<double> par){
  return -7.*pow(t,6);
}


int main(int argc, char **argv) {
  #ifdef USING_PETSC
    petscInitialize(argc, argv);
  #endif
  double xa = 0.;


  feFunction *funSol = new feFunction(fSol, {});
  feFunction *funSource = new feFunction(fSource, {});

  int nIter = 7;
  std::vector<double> normL2_BDF1(2 * nIter, 0.0);
  std::vector<double> normL2_DC2F(2 * nIter, 0.0);
  std::vector<double> normL2_DC3F(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);
  std::vector<int> TT;
  TT.resize(nIter);
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
    int nTimeSteps =5 * pow(2, iter);
    TT[iter] = nTimeSteps;
    feSolution *solDC3F = new feSolution(mesh, fespace, feEssBC, metaNumber);
    // Formes (bi)lineaires
    int nQuad = 10; // TODO : change to deg
    std::vector<feSpace *> spaceSource1D_U = {&U_M1D};
    std::vector<feSpace *> spaceMasse1D_U = {&U_M1D};

    feBilinearForm *source_U_M1D =
      new feBilinearForm(spaceSource1D_U, mesh, nQuad, new feSysElm_0D_Source(1.0, funSource));
    feBilinearForm *masse_U_M1D =
      new feBilinearForm(spaceMasse1D_U, mesh, nQuad, new feSysElm_0D_Masse(1.0, nullptr));

    std::vector<feBilinearForm *> formMatrices = {masse_U_M1D};
    std::vector<feBilinearForm *> formResiduals = {masse_U_M1D, source_U_M1D};
    // Norme de la solution
    feNorm *normBDF1 = new feNorm(&U_M1D, mesh, nQuad, funSol);
    feNorm *normDC2F = new feNorm(&U_M1D, mesh, nQuad, funSol);
    feNorm *normDC3F = new feNorm(&U_M1D, mesh, nQuad, funSol);
    std::vector<feNorm *> norms = {normBDF1, normDC2F, normDC3F};
    // Systeme lineaire
    feLinearSystemPETSc *linearSystem =
      new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_PETSC
    linearSystem->initialize();
    // Resolution
    feTolerances tol{1e-10, 1e-9, 10};
    DC3FSolver solver(tol, metaNumber, linearSystem, solDC3F, norms, mesh, t0, t1, nTimeSteps);
    solver.makeSteps(nTimeSteps, fespace);
    std::vector<double> &normL2BDF1 = solver.getNorm(0);
    std::vector<double> &normL2DC2F = solver.getNorm(1);
    std::vector<double> &normL2DC3F = solver.getNorm(2);

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
  printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \n", "nTimeSteps", "nElm", "||E_BDF1||", "Taux BDF1", "||E_DC2F||", "Taux DC2F" , "||E_DC3F||", "Taux DC3F");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \n", TT[i], nElm[i], normL2_BDF1[2 * i], normL2_BDF1[2 * i + 1], normL2_DC2F[2 * i], normL2_DC2F[2 * i + 1], normL2_DC3F[2 * i], normL2_DC3F[2 * i + 1]);
#ifdef USING_PETSC
  petscFinalize();
#endif
  return 0;
}
