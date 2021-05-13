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

double fSol(const double t, const std::vector<double> &x, const std::vector<double> par){
  double c2 = par[1];
  double x1 = x[0];
  return x1*(x1*x1-25.)*pow(t,6) + 1. + x1*x1*x1*pow(t,c2);
}

double fSource(const double t, const std::vector<double> &x, const std::vector<double> par){

  double c1 = par[0];
  double c2 = par[1];
  double beta = (c2 == 0.) ? 0. : c2*pow(t,c2-1.);
  double x1 = x[0]; 
  return -x1*(x1*x1-25.)*6.*pow(t,5.) + c1*6.*x1*pow(t,6.) + c1*6.*x1*pow(t,c2) - x1*x1*x1*beta;
}

int main(int argc, char** argv){
  double xa = 0.;
  double xb = 5.;

  double kd = 0.1;
  std::vector<double> par = {kd, 6.};
  feFunction *funSol = new feFunction(fSol, par);
  feFunction *funSource = new feFunction(fSource, par);

  int nIter = 5;
  std::vector<double> normL2(2*nIter, 0.0);
  std::vector<int> nElm(nIter, 0);

  for(int iter = 0; iter < nIter; ++iter){
    nElm[iter] = 20 * pow(2, iter);
    // Maillage
    feMesh1DP1 *mesh = new feMesh1DP1(xa,xb,nElm[iter],"BXA","BXB","M1D");
    // Espaces d'interpolation
    feSpace1DP0 U_BXA = feSpace1DP0(mesh, "U", "BXA", funSol);
    feSpace1DP0 U_BXB = feSpace1DP0(mesh, "U", "BXB", funSol);
    feSpace1DP2 U_M1D = feSpace1DP2(mesh, "U", "M1D", funSol);
    std::vector<feSpace*> fespace = {&U_BXA, &U_BXB, &U_M1D};
    std::vector<feSpace*> feEssBC = {&U_BXA, &U_BXB};
    // Numerotations
    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);
    // Solution
    double t0 = 0.;
    double t1 = 1.;
    int nTimeSteps = 20 * pow(2, iter);
    feSolution *sol = new feSolution(mesh, fespace, feEssBC, metaNumber);
    sol->initializeTemporalSolution(t0,t1,nTimeSteps);
    sol->initializeUnknowns(mesh, metaNumber);
    sol->initializeEssentialBC(mesh, metaNumber);
    // Formes (bi)lineaires
    int degQuad = 5; 
    std::vector<feSpace*> spaceDiffusion1D_U = {&U_M1D};
    std::vector<feSpace*> spaceSource1D_U = {&U_M1D};
    std::vector<feSpace*> spaceMasse1D_U = {&U_M1D};

    feBilinearForm *diff_U_M1D = new feBilinearForm(spaceDiffusion1D_U, mesh, degQuad, new feSysElm_1D_Diffusion(kd, nullptr));
    feBilinearForm *source_U_M1D = new feBilinearForm(spaceSource1D_U, mesh, degQuad, new feSysElm_1D_Source(1.0, funSource));
    feBilinearForm *masse_U_M1D = new feBilinearForm(spaceMasse1D_U, mesh, degQuad, new feSysElm_1D_Masse(1.0, nullptr));

    std::vector<feBilinearForm*> formMatrices  = {diff_U_M1D, masse_U_M1D};
    std::vector<feBilinearForm*> formResiduals  = {diff_U_M1D, masse_U_M1D, source_U_M1D};
    // Norme de la solution
    feNorm *norm = new feNorm(&U_M1D, mesh, degQuad);
    // Systeme lineaire
    feLinearSystemPETSc *linearSystem = new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_PETSC
    linearSystem->initialize();
    // Resolution
    feTolerances tol{1e-10, 1e-10, 10};
    std::vector<double> normL2BDF(nTimeSteps,0.0);
    solveBDF2(normL2BDF, tol, metaNumber, linearSystem, formMatrices, formResiduals, sol, norm, mesh);
    normL2[2*iter] = *std::max_element(normL2BDF.begin(), normL2BDF.end());
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
  for(int i = 1; i < nIter; ++i){
    normL2[2*i+1] = log(normL2[2*(i-1)]/normL2[2*i])/log(2.);
  }
  printf("%12s \t %12s \t %12s\n", "nElm", "||E||", "Taux BDF2");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12.6e \t %12.6e\n", nElm[i], normL2[2*i], normL2[2*i+1]);

  return 0;
}
