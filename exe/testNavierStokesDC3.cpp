#include <iostream>
#include <fstream>

#include "feNG.h"
#include "feFunction.h"
#include "feMesh.h"
#include "feQuadrature.h"
#include "feSpace.h"
#include "feSpaceTriangle.h"
#include "feCncGeo.h"
#include "feNumber.h"
#include "feSolution.h"
#include "feSysElm.h"
#include "feBilinearForm.h"
#include "feSolver.h"
#include "feExporter.h"

#ifdef HAVE_PETSC
#include "feLinearSystemPETSc.h"
#endif
#ifdef HAVE_MKL
#include "feLinearSystemMklPardiso.h"
#endif


double fSolU(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  double F = exp(-2.0 * t);
  return F * sin(x) * cos(y); // Taylor-Green
  // return pow(x,4) + pow(y,4);
}

double fSolV(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  double F = exp(-2.0 * t);
  return -F * cos(x) * sin(y); // Taylor-Green
  // return -4.0 * pow(x,3) * y;
}

double fSolP(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  // double rho = par[0];
  double F = exp(-2.0 * t);
  return F * F * 1.0 / 4.0 * (cos(2. * x) + cos(2. * y)); // Taylor-Green
  // return pow(x,4) * pow(y,4);
}

double fOne(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  return 1.0;
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  return 0.0;
}

void fSource(const double t, const std::vector<double> &pos, const std::vector<double> &par,
             std::vector<double> &res) {
  // double rho = par[0];
  // double mu = par[1];
  // double nu = mu/rho;
  double x = pos[0];
  double y = pos[1];

  // res[0] = -mu * (12*(x*x*y*y*y*y + x*x*x*x*y*y)) + 2*x*y*y;
  // res[1] = -mu * (-8./5. *x*y*y*y*(10*x*x + 3*y*y)) + 2*y*x*x;

  // res[0] = 4.0*pow(x,3)*(pow(x,4) + pow(y,4)) - 12.0*pow(x,3)*pow(y,4) - 12.0*x*x - 12.0*y*y;
  // res[1] = 4.0*pow(x,4)*pow(y,3) + 24.0*x*y + 16.0*pow(x,6)*y - 12.0*x*x*y*(pow(x,4) + pow(y,4));
  if(fabs(t) <= 1e-10) {
    res[0] = 2.0 * cos(y) * sin(x);
    res[1] = -2.0 * cos(x) * sin(y);
  } else {
    res[0] = 0.0;
    res[1] = 0.0;
  }
}

int main(int argc, char **argv) {
#ifdef USING_PETSC
  petscInitialize(argc, argv);
#endif
  double rho = 1.0;
  double mu = 1.0;
  std::vector<double> stokesParam = {rho, mu};

  feFunction *funSolU = new feFunction(fSolU, {});
  feFunction *funSolV = new feFunction(fSolV, {});
  feFunction *funSolP = new feFunction(fSolP, stokesParam);
  feFunction *funZero = new feFunction(fZero, {});
  feFunction *funOne = new feFunction(fOne, {});
  feVectorFunction *funSource = new feVectorFunction(fSource, stokesParam);

  int nIter = 2;
  std::vector<double> normL2_U(2 * nIter, 0.0);
  std::vector<double> normL2_V(2 * nIter, 0.0);
  std::vector<double> normL2_P(2 * nIter, 0.0);
  std::vector<double> maxNormL2BDF_U(2 * nIter, 0.0);
  std::vector<double> maxNormL2BDF_V(2 * nIter, 0.0);
  std::vector<double> maxNormL2BDF_P(2 * nIter, 0.0);
  std::vector<double> maxNormL2DC3_U(2 * nIter, 0.0);
  std::vector<double> maxNormL2DC3_V(2 * nIter, 0.0);
  std::vector<double> maxNormL2DC3_P(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);

  for(int iter = 0; iter < nIter; ++iter) {
    // std::string meshName = "../../data/squareTaylorGreen2.msh";
    std::string meshName = "../../data/Square/squareNS" + std::to_string(iter+1) + ".msh";
    // std::string meshName = "../../data/squarePression" + std::to_string(iter+1) + ".msh";

    feMesh2DP1 *mesh = new feMesh2DP1(meshName, false);
    nElm[iter] = mesh->getNbInteriorElems();

    feSpace1DP3 U_angle = feSpace1DP3(mesh, "U", "Angle", funSolU);
    feSpace1DP3 U_haut = feSpace1DP3(mesh, "U", "Haut", funSolU);
    feSpace1DP3 U_gauche = feSpace1DP3(mesh, "U", "Gauche", funSolU);
    feSpaceTriP3 U_surface = feSpaceTriP3(mesh, "U", "Surface", funSolU);

    feSpace1DP3 V_angle = feSpace1DP3(mesh, "V", "Angle", funSolV);
    feSpace1DP3 V_haut = feSpace1DP3(mesh, "V", "Haut", funSolV);
    feSpace1DP3 V_gauche = feSpace1DP3(mesh, "V", "Gauche", funSolV);
    feSpaceTriP3 V_surface = feSpaceTriP3(mesh, "V", "Surface", funSolV);

    feSpace1DP0 P_point = feSpace1DP0(mesh, "P", "PointPression", funSolP);
    feSpaceTriP2 P_surface = feSpaceTriP2(mesh, "P", "Surface", funSolP);

    std::vector<feSpace *> fespace = {&U_angle, &U_haut,   &U_gauche,  &U_surface, &V_angle,
                                      &V_haut,  &V_gauche, &V_surface, &P_point,   &P_surface};

    std::vector<feSpace *> feEssBC = {&U_angle, &U_haut,   &U_gauche, &V_angle,
                                      &V_haut,  &V_gauche, &P_point};

    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);
    // metaNumber->printNumberings();

    feSolution *sol = new feSolution(mesh, fespace, feEssBC, metaNumber);
    // sol->initializeUnknowns(mesh, metaNumber);
    // sol->initializeEssentialBC(mesh, metaNumber);

    // Formes (bi)lineaires
    int nQuad = 9;
    std::vector<feSpace *> spacesNS2D = {&U_surface, &V_surface, &P_surface};

    // feBilinearForm *NS2D = new feBilinearForm(spacesNS2D, mesh, nQuad, new
    // feSysElm_2D_NavierStokes(stokesParam, funSource));
    feBilinearForm *NS2D = new feBilinearForm(spacesNS2D, mesh, nQuad,
                                              new feSysElm_2D_NavierStokes(stokesParam, nullptr));

    std::vector<feBilinearForm *> formMatrices = {NS2D};
    std::vector<feBilinearForm *> formResiduals = {NS2D};

    feNorm *normU = new feNorm(&U_surface, mesh, nQuad, funSolU);
    // normU->computeL2Norm(metaNumber, sol, mesh);
    // std::cout<<"fooU "<<normU->getNorm()<<std::endl;
    feNorm *normV = new feNorm(&V_surface, mesh, nQuad, funSolV);
    // normV->computeL2Norm(metaNumber, sol, mesh);
    // std::cout<<"fooV "<<normV->getNorm()<<std::endl;
    feNorm *normP = new feNorm(&P_surface, mesh, nQuad, funSolP);
    // normP->computeL2Norm(metaNumber, sol, mesh);
    // std::cout<<"fooP "<<normP->getNorm()<<std::endl;
    std::vector<feNorm *> norms = {normU, normV, normP};
    // Systeme lineaire

    // feLinearSystemPETSc *linearSystem =
    // new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_MKL    
    feLinearSystemMklPardiso *linearSystem =
    new feLinearSystemMklPardiso( formMatrices, formResiduals, metaNumber, mesh);
  #ifdef HAVE_PETSC
      // linearSystem->initialize();
      // Resolution
      feTolerances tol{1e-9, 1e-9, 20};

      double t0 = 0.;
      double t1 = 1.0;
      int nTimeSteps = 20 * pow(2, iter);
      // sol->initializeTemporalSolution(t0, t1, nTimeSteps);
      std::vector<double> normL2BDF(3 * nTimeSteps, 0.0);
      std::vector<double> normL2DC3(3 * nTimeSteps, 0.0);
      // solveBDF2(normL2BDF, tol, metaNumber, linearSystem, formMatrices, formResiduals, sol, norms,
      // mesh, fespace);
      // solveDC3(normL2BDF, normL2DC3, tol, metaNumber, linearSystem, formMatrices, formResiduals,
      // sol,
      //          norms, mesh, fespace);
      std::string CodeIni = "BDF1/DC"; // Define the way of initialization |"SolEx"->for exact
                                       // solution|  |"BDF1/DCF"->using only initial conditions|
      DC3Solver solver(tol, metaNumber, linearSystem, sol, norms, mesh, t0, t1, nTimeSteps, CodeIni);
      solver.makeSteps(nTimeSteps, fespace);
      // normL2_U[2*iter] = *std::max_element(normL2BDF.begin(), normL2BDF.end());
      // maxNormL2BDF[2*iter] = *std::max_element(normL2BDF.begin(), normL2BDF.end());
      for(int i = 0; i < nTimeSteps; ++i) {
        maxNormL2DC3_U[2 * iter] = fmax(maxNormL2DC3_U[2 * iter], normL2DC3[3 * i]);
        maxNormL2DC3_V[2 * iter] = fmax(maxNormL2DC3_V[2 * iter], normL2DC3[3 * i + 1]);
        maxNormL2DC3_P[2 * iter] = fmax(maxNormL2DC3_P[2 * iter], normL2DC3[3 * i + 2]);
      }

      // std::string vtkFile = "../../data/solutionManufacturee" + std::to_string(iter+1) + ".vtk";
      // std::string vtkFile = "../../data/taylorGreenUnsteady" + std::to_string(iter+1) + ".vtk";
      // feExporterVTK writer(vtkFile, mesh, sol, metaNumber, fespace);

      normU->computeL2Norm(metaNumber, sol, mesh);
      normV->computeL2Norm(metaNumber, sol, mesh);
      normP->computeL2Norm(metaNumber, sol, mesh);
      normL2_U[2 * iter] = normU->getNorm();
      normL2_V[2 * iter] = normV->getNorm();
      normL2_P[2 * iter] = normP->getNorm();
  #endif
    // std::string vtkFile = "../../data/TestNS_DC3.vtk";
    // feExporterVTK writer(vtkFile, mesh, sol, metaNumber, fespace);
    delete normU;
    delete normV;
    delete normP;
    delete linearSystem;
    delete NS2D;
    delete sol;
    delete metaNumber;
    delete mesh;
#endif      
  }
  delete funSource;
  delete funSolU;
  delete funSolV;
  delete funSolP;
  delete funZero;
  delete funOne;

  // Calcul du taux de convergence
  for(int i = 1; i < nIter; ++i) {
    normL2_U[2 * i + 1] =
      -log(normL2_U[2 * i] / normL2_U[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
    normL2_V[2 * i + 1] =
      -log(normL2_V[2 * i] / normL2_V[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
    normL2_P[2 * i + 1] =
      -log(normL2_P[2 * i] / normL2_P[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
  }
  for(int i = 1; i < nIter; ++i) {
    maxNormL2DC3_U[2 * i + 1] = log(maxNormL2DC3_U[2 * (i - 1)] / maxNormL2DC3_U[2 * i]) / log(2.);
    maxNormL2DC3_V[2 * i + 1] = log(maxNormL2DC3_V[2 * (i - 1)] / maxNormL2DC3_V[2 * i]) / log(2.);
    maxNormL2DC3_P[2 * i + 1] = log(maxNormL2DC3_P[2 * (i - 1)] / maxNormL2DC3_P[2 * i]) / log(2.);
  }
  printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s\n", "nElm", "||u-uh||", "tauxU",
         "||v-vh||", "tauxV", "||p-ph||", "tauxP");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n", nElm[i],
           normL2_U[2 * i], normL2_U[2 * i + 1], normL2_V[2 * i], normL2_V[2 * i + 1],
           normL2_P[2 * i], normL2_P[2 * i + 1]);

  printf("%12s \t %12s \t %12s\t %12s \t %12s\t %12s \t %12s\n", "nElm", "||E-DC3_U||",
         "Taux DC3_U", "||E-DC3_V||", "Taux DC3_V", "||E-DC3_P||", "Taux DC3_P");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12.6e \t %12.6e\t %12.6e \t %12.6e\t %12.6e \t %12.6e\n", nElm[i],
           maxNormL2DC3_U[2 * i], maxNormL2DC3_U[2 * i + 1], maxNormL2DC3_V[2 * i],
           maxNormL2DC3_V[2 * i + 1], maxNormL2DC3_P[2 * i], maxNormL2DC3_P[2 * i + 1]);
#ifdef USING_PETSC
  petscFinalize();
#endif
  return 0;
}
