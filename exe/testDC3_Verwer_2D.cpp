#include <iostream>
#include <functional>
#include <vector>
#include <cmath>

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
#include "feLinearSystemMklPardiso.h"
#include "feExporter.h"

#ifdef HAVE_PETSC
#include "feLinearSystemPETSc.h"
#endif
#ifdef HAVE_MKL
#include "feLinearSystemMklPardiso.h"
#endif

// PB 7-71

// double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   return x1*x1*x1* (x1 * x1 - 25.) * pow(t, 6) +  pow(x1, 5) +1.;
// }

// double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double c1 = par[0];
//   double x1 = x[0];
//   return -x1*x1*x1* (x1 * x1 - 25.) * 6. * pow(t, 5.) - c1 * 150. * x1 * pow(t, 6.) +
//          c1 * 20. * x1 *
//          x1*x1* (pow(t, 6) +1. );
// }

// double flambda_A(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double c1 = par[0];
//   return  c1* (5.*pow(x1,4.)*(pow(t, 6)+1) - 75.*x1*x1*pow(t, 6)) ;
// }
// double flambda_B(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double c1 = par[0];
//   return - c1* (5.*pow(x1,4.)*(pow(t, 6)+1) - 75.*x1*x1*pow(t, 6)) ;
// }
// double fSolDot(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double c1 = par[0];
//   // double a = par[1];
//   // double beta = (a == 0.) ? 0. : a * pow(t, a - 1.);
//   double x1 = x[0];
//   return 6.*x1*x1*x1* (x1 * x1 - 25.) * pow(t, 5) ;
// }

// PB 7-73

// double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   return x1*x1*x1* (x1 * x1 - 25.) * pow(t, 6) +  pow(x1, 5)*pow(t, 6) +1.;
// }

// double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double c1 = par[0];
//   double x1 = x[0];
//   return -x1*x1*x1* (x1 * x1 - 25.) * 6. * pow(t, 5.) - 6.* pow(x1, 5)*pow(t, 5) - c1 * 150. * x1
//   * pow(t, 6.) +
//          c1 * 40. * x1 *
//          x1*x1* pow(t, 6);
// }
// double flambda_A(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double c1 = par[0];
//   return  c1* (10.*pow(x1,4.)*pow(t, 6) - 75.*x1*x1*pow(t, 6)) ;
// }
// double flambda_B(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double c1 = par[0];
//   return - c1* (10.*pow(x1,4.)*pow(t, 6) - 75.*x1*x1*pow(t, 6)) ;
// }

// double fSolDot(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   return 6.*x1*x1*x1* (x1 * x1 - 25.) * pow(t, 5) +  6.*pow(x1, 5)*pow(t, 5);
// }

// double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double y1 = x[1];
//   int a = par[1];
//   return pow(t, a) * sin(x1*(2*M_PI))* cos(y1*(2*M_PI));
// }

// double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double y1 = x[1];
//   double c1 = par[0];
//   double a = par[1];
//   double beta = (a == 0.) ? 0. : a * pow(t, a - 1.);
//   return -beta * sin(x1*(2*M_PI))* cos(y1*(2*M_PI)) - c1 *2*pow(t, a) * sin(x1*(2*M_PI))* cos(y1*(2*M_PI))*((2*M_PI)*(2*M_PI));
// }

// double flambda_angle(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double y1 = x[1];
//   double c1 = par[0];
//   double a = par[1];
//   if(y1==0.) return c1*pow(t, a)*sin(x1*(2*M_PI))*sin(y1*(2*M_PI))*(2*M_PI);
//   else return -c1 * pow(t, a) * cos(x1*(2*M_PI))*cos(y1*(2*M_PI))*(2*M_PI);
// }


// double flambda_gauche(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double y1 = x[1];
//   double c1 = par[0];
//   double a = par[1];
//   return c1 * pow(t, a) *cos(x1*(2*M_PI))*cos(y1*(2*M_PI))*(2*M_PI);
// }

// double flambda_haut(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double x1 = x[0];
//   double y1 = x[1];
//   double c1 = par[0];
//   double a = par[1];
//   return c1 * pow(t, a) *sin(x1*(2*M_PI))*sin(y1*(2*M_PI))*(2*M_PI);
// }

// double fSolDot(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   double c1 = par[0];
//   double a = par[1];
//   double beta = (a == 0.) ? 0. : a * pow(t, a - 1.);
//   double x1 = x[0];
//   double y1 = x[1];
//   return beta * sin(x1*(2*M_PI))* cos(y1*(2*M_PI)) ;
// }

double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double y1 = x[1];
  int a = par[1];
  return pow(t, a) * (x1 * x1 + 1)*y1;
}

double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double y1 = x[1];
  double c1 = par[0];
  double a = par[1];
  double beta = (a == 0.) ? 0. : a * pow(t, a - 1.);
  return -beta * (x1 * x1 + 1)*y1  + c1 * 2 * pow(t, a)*y1;
}

double flambda_angle(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double y1 = x[1];
  double c1 = par[0];
  double a = par[1];
  if(y1==0.) return c1*pow(t, a)*(x1 * x1 + 1);
  else return -c1 * pow(t, a) * (2 * x1)* y1;
}


double flambda_haut(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double y1 = x[1];
  double c1 = par[0];
  double a = par[1];
  return -c1*pow(t, a)*(x1 * x1 + 1);
}

double flambda_gauche(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double x1 = x[0];
  double y1 = x[1];
  double c1 = par[0];
  double a = par[1];
  return c1 * pow(t, a) *y1* (2 * x1);
}

double fSolDot(const double t, const std::vector<double> &x, const std::vector<double> par) {
  double c1 = par[0];
  double a = par[1];
  double beta = (a == 0.) ? 0. : a * pow(t, a - 1.);
  double x1 = x[0];
  double y1 = x[1];
  return beta * (x1 * x1 + 1)*y1;
}

double f0(const double t, const std::vector<double> &x, const std::vector<double> par) {
  return 0;
}


int main(int argc, char **argv) {
#ifdef USING_PETSC
  petscInitialize(argc, argv);
#endif
  double exposant = 5.;
  double kd = 0.1;
  std::vector<double> par = {kd, exposant};
  feFunction *funSol = new feFunction(fSol, par);
  feFunction *funSolDot = new feFunction(fSolDot, par);
  feFunction *funSource = new feFunction(fSource, par);
  feFunction *funLambda_angle = new feFunction(flambda_angle, par);
  feFunction *funLambda_gauche = new feFunction(flambda_gauche, par);
  feFunction *funLambda_haut = new feFunction(flambda_haut, par);
  feFunction *fun0 = new feFunction(f0, par);


  int nIter = 1;
  std::vector<double> normL2_BDF2(2 * nIter, 0.0);
  std::vector<double> normL2_DC3(2 * nIter, 0.0);
  std::vector<double> normBDF2_angle(2 * nIter, 0.0);
  std::vector<double> normBDF2_gauche(2 * nIter, 0.0);
  std::vector<double> normDC3_angle(2 * nIter, 0.0);
  std::vector<double> normDC3_gauche(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);
  std::vector<int> TT;
  TT.resize(nIter);

  for(int iter = 0; iter < nIter; ++iter) {
    // std::string meshName = "../../data/Circle/circle.msh";
    // std::string meshName = "../../data/squareTest2.msh";
    std::string meshName = "../../data/Square/square" + std::to_string(iter*2 + 1) + "Msh2.msh";
    // Maillage
     feMesh2DP1 *mesh = new feMesh2DP1(meshName, false);
     nElm[iter] = mesh->getNbInteriorElems();
    // Espaces d'interpolation

    feSpace1DP4 U_angle = feSpace1DP4(mesh, "U", "Angle", fun0);
    feSpace1DP4 U_haut = feSpace1DP4(mesh, "U", "Haut", fun0);
    feSpaceTriP4 U_surface = feSpaceTriP4(mesh, "U", "Surface", funSol);
    feSpace1DP4 U_gauche = feSpace1DP4(mesh, "U", "Gauche", fun0);

    feSpace1DP4 L_angle = feSpace1DP4(mesh, "L", "Angle", funLambda_angle);
    feSpace1DP4 L_haut = feSpace1DP4(mesh, "L", "Haut", funLambda_haut);
    feSpace1DP4 L_gauche = feSpace1DP4(mesh, "L", "Gauche", funLambda_gauche);

    feSpace1DP4 V_angle = feSpace1DP4(mesh, "V", "Angle", fun0);
    feSpace1DP4 V_haut = feSpace1DP4(mesh, "V", "Haut", fun0);
    feSpace1DP4 V_gauche = feSpace1DP4(mesh, "V", "Gauche", fun0);


    std::vector<feSpace *> fespace = {&U_angle, &U_surface, &U_haut, &U_gauche, &L_angle, &L_haut, &L_gauche, &V_angle, &V_haut, &V_gauche};
    std::vector<feSpace *> feEssBC = {};
    // Numerotations
    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);
    // Solution
    double t0 = 0.;
    double t1 = 1.;
    int nTimeSteps = 10 * pow(2, iter);
    TT[iter] = nTimeSteps;
    feSolution *solDC3 = new feSolution(mesh, fespace, feEssBC, metaNumber);
   
    // Formes (bi)lineaires
    int degQuad = 11;
    std::vector<feSpace *> spaceDiffusion2D_U = {&U_surface};
    std::vector<feSpace *> spaceSource2D_U = {&U_surface};
    std::vector<feSpace *> spaceMasse2D_U = {&U_surface};
    std::vector<feSpace *> spaceWeak_angle = {&U_angle, &L_angle, &V_angle};
    std::vector<feSpace *> spaceWeak_haut = {&U_haut, &L_haut, &V_haut};
    std::vector<feSpace *> spaceWeak_gauche = {&U_gauche, &L_gauche, &V_gauche};

    feBilinearForm *diff_U_M2D =
      new feBilinearForm(spaceDiffusion2D_U, mesh, degQuad, new feSysElm_2D_Diffusion(kd, nullptr));
    feBilinearForm *source_U_M2D =
      new feBilinearForm(spaceSource2D_U, mesh, degQuad, new feSysElm_2D_Source(1.0, funSource));
    feBilinearForm *masse_U_M2D =
      new feBilinearForm(spaceMasse2D_U, mesh, degQuad, new feSysElm_2D_Masse(1.0, nullptr));
      feBilinearForm *weakBC_U_angle =
      new feBilinearForm(spaceWeak_angle, mesh, degQuad, new feSysElm_1D_weakBC_edo1(1.0, funSolDot));
    feBilinearForm *weakBC_U_haut =
      new feBilinearForm(spaceWeak_haut, mesh, degQuad, new feSysElm_1D_weakBC_edo1(1.0, funSolDot));
      feBilinearForm *weakBC_U_gauche =
      new feBilinearForm(spaceWeak_gauche, mesh, degQuad, new feSysElm_1D_weakBC_edo1(1.0, funSolDot));

    std::vector<feBilinearForm *> formMatrices = {diff_U_M2D, masse_U_M2D, weakBC_U_angle ,weakBC_U_haut ,weakBC_U_gauche};
    std::vector<feBilinearForm *> formResiduals = {diff_U_M2D, masse_U_M2D, source_U_M2D, weakBC_U_angle ,weakBC_U_haut ,weakBC_U_gauche};

    // Norme de la solution
    feNorm *normBDF2 = new feNorm(&U_surface, mesh, degQuad, funSol);
    feNorm *normDC3 = new feNorm(&U_surface, mesh, degQuad, funSol);
    feNorm *normBDF2Lambda_angle = new feNorm(&L_angle, mesh, degQuad, funLambda_angle);
    feNorm *normBDF2Lambda_gauche = new feNorm(&L_gauche, mesh, degQuad, funLambda_gauche);
    feNorm *normDC3Lambda_angle = new feNorm(&L_angle, mesh, degQuad, funLambda_angle);
    feNorm *normDC3Lambda_gauche = new feNorm(&L_gauche, mesh, degQuad, funLambda_gauche);
    std::vector<feNorm *> norms = {normBDF2, normDC3 , normBDF2Lambda_angle ,normBDF2Lambda_gauche, normDC3Lambda_angle ,normDC3Lambda_gauche};
  
    // Systeme lineaire
    // feLinearSystemPETSc *linearSystem =
    // new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_PETSC
      feLinearSystemMklPardiso *linearSystem =
      new feLinearSystemMklPardiso( formMatrices, formResiduals, metaNumber, mesh);
  #ifdef HAVE_PETSC
      // linearSystem->initialize();
      // Resolution
      feTolerances tol{1e-9, 1e-9, 20};
      std::string CodeIni = "BDF1/DCF"; // Define the way of initialization |"SolEx"->for exact
                                        // solution|  |"BDF1/DCF"->using only initial conditions|
      DC3Solver solver(tol, metaNumber, linearSystem, solDC3, norms, mesh, t0, t1, nTimeSteps,
                       CodeIni);
      solver.makeSteps(nTimeSteps, fespace);
      std::vector<double> &normL2BDF2 = solver.getNorm(0);
      std::vector<double> &normL2DC3 = solver.getNorm(1);
      std::vector<double> &normL2_BDF2_Lambda_angle = solver.getNorm(2);
      std::vector<double> &normL2_BDF2_Lambda_gauche = solver.getNorm(3);
      std::vector<double> &normL2_DC3_Lambda_angle = solver.getNorm(4);
      std::vector<double> &normL2_DC3_Lambda_gauche = solver.getNorm(5);
      normL2_BDF2[2 * iter] = *std::max_element(normL2BDF2.begin(), normL2BDF2.end());
      normL2_DC3[2 * iter] = *std::max_element(normL2DC3.begin(), normL2DC3.end());
      normBDF2_angle[2 * iter] =
        *std::max_element(normL2_BDF2_Lambda_angle.begin(), normL2_BDF2_Lambda_angle.end());
      normBDF2_gauche[2 * iter] =
        *std::max_element(normL2_BDF2_Lambda_gauche.begin(), normL2_BDF2_Lambda_gauche.end());
      normDC3_angle[2 * iter] = *std::max_element(normL2_DC3_Lambda_angle.begin(), normL2_DC3_Lambda_angle.end());
      normDC3_gauche[2 * iter] = *std::max_element(normL2_DC3_Lambda_gauche.begin(), normL2_DC3_Lambda_gauche.end());
  #endif
      std::string vtkFile = "../../data/TestDC3.vtk";
      feExporterVTK writer(vtkFile, mesh, solDC3, metaNumber, fespace);
      delete normBDF2;
      delete normDC3;
      delete normBDF2Lambda_gauche;
      delete normBDF2Lambda_angle;
      delete normDC3Lambda_gauche;
      delete normDC3Lambda_angle;
      delete linearSystem;
      delete masse_U_M2D;
      delete source_U_M2D;
      delete diff_U_M2D;
      delete solDC3;
      delete metaNumber;
      delete mesh;
#endif
  }
  delete funSource;
  delete funSol;

  // Calcul du taux de convergence
  for(int i = 1; i < nIter; ++i) {
    normL2_BDF2[2 * i + 1] = log(normL2_BDF2[2 * (i - 1)] / normL2_BDF2[2 * i]) / log(2.);
    normL2_DC3[2 * i + 1] = log(normL2_DC3[2 * (i - 1)] / normL2_DC3[2 * i]) / log(2.);
    normBDF2_angle[2 * i + 1] = log(normBDF2_angle[2 * (i - 1)] / normBDF2_angle[2 * i]) / log(2.);
    normBDF2_gauche[2 * i + 1] = log(normBDF2_gauche[2 * (i - 1)] / normBDF2_gauche[2 * i]) / log(2.);
    normDC3_angle[2 * i + 1] = log(normDC3_angle[2 * (i - 1)] / normDC3_angle[2 * i]) / log(2.);
    normDC3_gauche[2 * i + 1] = log(normDC3_gauche[2 * (i - 1)] / normDC3_gauche[2 * i]) / log(2.);
  }
  printf("%12s \t %12s \t %12s \t %12s\t %12s \t %12s\n", "nSteps", "nElm", "||E-BDF||",
         "Taux BDF2", "||E-DC3||", "Taux DC3");
  for(int i = 0; i < nIter; ++i)

    printf("%12d \t %12d \t %12.6e \t %12.6e\t %12.6e \t %12.6e\n", TT[i], nElm[i],
           normL2_BDF2[2 * i], normL2_BDF2[2 * i + 1], normL2_DC3[2 * i], normL2_DC3[2 * i + 1]);
  printf("\n");
  printf("\n");
  printf("Lambda based on BDF2");
  printf("\n");
  printf("%12s \t %12s \t %16s \t %16s \t %16s \t %16s \n", "nTimeSteps", "nElm", "||E_lambda_angle||",
         "Taux lambda_angle", "||E_lambda_gauche||", "Taux lambda_gauche");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %16.6e \t %16.6e \t %16.6e \t %16.6e  \n", TT[i], nElm[i],
           normBDF2_angle[2 * i], normBDF2_angle[2 * i + 1], normBDF2_gauche[2 * i], normBDF2_gauche[2 * i + 1]);
  printf("\n");
  printf("\n");
  printf("Lambda based on DC3");
  printf("\n");
  printf("%12s \t %12s \t %16s \t %16s \t %16s \t %16s \n", "nTimeSteps", "nElm", "||E_lambda_angle||",
         "Taux lambda_angle", "||E_lambda_gauche||", "Taux lambda_gauche");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %16.6e \t %16.6e \t %16.6e \t %16.6e  \n", TT[i], nElm[i],
           normDC3_angle[2 * i], normDC3_angle[2 * i + 1], normDC3_gauche[2 * i], normDC3_gauche[2 * i + 1]);
#ifdef USING_PETSC
  petscFinalize();
#endif
  return 0;
}
