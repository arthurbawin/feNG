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
  double q = par[2];
  double p = q -1.0;
  return -pow(y,p) * (1.-pow(x,q));
}

double fSolV(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  double q = par[2];
  double p = q -1.0;
  return pow(x,p) * (1.-pow(y,q));
}
void fVecVeloc(const double t, const std::vector<double> &pos, const std::vector<double> &par,
             std::vector<double> &res) {
  double x = pos[0];
  double y = pos[1];
  double q = par[2];
  double p = q -1.0;
  res[0] = -pow(y,p) * (1.-pow(x,q));       //u
  res[1] = pow(x,p) * (1.-pow(y,q));        //v
  res[2] = pow(y,p)*q*pow(x,p);             //dudx
  res[3] = -p*pow(y,p-1.) * (1-pow(x,q));   //dudy
  res[4] = p*pow(x,p-1.) * (1-pow(y,q));    //dvdx
  res[5] = -pow(x,p)*q*pow(y,p);            //dvdy

}
double fSolP(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  double q = par[2];
  double p = q -1.0;
  // return 2. +x +y;
  return 1;
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  return 0.0;
}

void fSource(const double t, const std::vector<double> &pos, const std::vector<double> &par,
             std::vector<double> &res) {
  double mu = par[1];
  double x = pos[0];
  double y = pos[1];
  double q = par[2];
  double p = q -1.0;
  double u = -pow(y,p) * (1-pow(x,q));
  double v = pow(x,p) * (1-pow(y,q));
  //dx
  double dudx   = pow(y,p)*q*pow(x,p);
  double d2udx2 = pow(y,p)*q*p*pow(x,p-1.);
  double dvdx   = p*pow(x,p-1.) * (1-pow(y,q));
  double d2vdx2 = p*(p-1.)*pow(x,p-2.) * (1-pow(y,q));
  //dy
  double dudy   = -p*pow(y,p-1.) * (1-pow(x,q));
  double d2udy2 = -p*(p-1.)*pow(y,p-2.) * (1-pow(x,q));
  double dvdy   = -pow(x,p)*q*pow(y,p);
  double d2vdy2 = -pow(x,p)*q*p*pow(y,p-1.);

  double dpdx = 0;
  double dpdy = 0;
  double rho = par[0];

  res[0] = rho*(u*dudx + v* dudy)  - mu* (d2udx2 + d2udy2)+ dpdx;
  res[1] = rho*(u*dvdx + v* dvdy) - mu*(d2vdx2 + d2vdy2 ) + dpdy;

  // res[0] = +u*dudx + v* dudy - mu* (d2udx2 + d2udy2) + (1/rho)*dpdx;
  // res[1] = +u*dvdx + v* dvdy - mu*(d2vdx2 + d2vdy2 ) + (1/rho)*dpdy;
}

int main(int argc, char **argv) {
#ifdef USING_PETSC
  petscInitialize(argc, argv);
#endif
  double rho = 1.0;
  double mu = 1.0;
  double q = 2.;
  std::vector<double> stokesParam = {rho, mu, q};

  feFunction *funSolU = new feFunction(fSolU, stokesParam);
  feFunction *funSolV = new feFunction(fSolV, stokesParam);
  feFunction *funSolP = new feFunction(fSolP, stokesParam);
  feFunction *funZero = new feFunction(fZero, {});
  feVectorFunction *funSource = new feVectorFunction(fSource, stokesParam);
  feVectorFunction *funVecVeloc = new feVectorFunction(fVecVeloc, stokesParam);

  int nIter = 2;
  std::vector<double> normL2_U(2 * nIter, 0.0);
  std::vector<double> normL2_V(2 * nIter, 0.0);
  std::vector<double> normL2_Vec(2 * nIter, 0.0);
  std::vector<double> normL2_P(2 * nIter, 0.0);
  std::vector<double> maxNormL2BDF_U(2 * nIter, 0.0);
  std::vector<double> maxNormL2BDF_V(2 * nIter, 0.0);
  std::vector<double> maxNormL2BDF_P(2 * nIter, 0.0);
  std::vector<double> maxNormL2DC3_U(2 * nIter, 0.0);
  std::vector<double> maxNormL2DC3_V(2 * nIter, 0.0);
  std::vector<double> maxNormL2DC3_P(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);

  for(int iter = 0; iter < nIter; ++iter) {
    std::string meshName = "../../data/Square/squareNS" + std::to_string(iter + 1) + ".msh";

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

    // Formes (bi)lineaires
    int nQuad = 16;
    std::vector<feSpace *> spacesNS2D = {&U_surface, &V_surface, &P_surface};

    // feBilinearForm *NS2D = new feBilinearForm(spacesNS2D, mesh, nQuad, new
    // feSysElm_2D_NavierStokes(stokesParam, funSource));
    feBilinearForm *NS2D = new feBilinearForm(spacesNS2D, mesh, nQuad,
                                              new feSysElm_2D_NavierStokes(stokesParam, funSource));

    std::vector<feBilinearForm *> formMatrices = {NS2D};
    std::vector<feBilinearForm *> formResiduals = {NS2D};
    
    
    
    
    feNorm *normU = new feNorm(&U_surface, mesh, nQuad, funSolU);
    feNorm *normV = new feNorm(&V_surface, mesh, nQuad, funSolV);
    std::vector<feSpace *> fespaceNorm = {&U_surface, &V_surface};
    feNorm *normVecVeloc = new feNorm(fespaceNorm, mesh, nQuad,funZero,funVecVeloc, "NormH1");
    feNorm *normP = new feNorm(&P_surface, mesh, nQuad, funSolP);
    std::vector<feNorm *> norms = {normU, normV,normVecVeloc, normP};
    // Systeme lineaire

    // feLinearSystemPETSc *linearSystem =
    // new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_MKL
    feLinearSystemMklPardiso *linearSystem =
      new feLinearSystemMklPardiso(formMatrices, formResiduals, metaNumber, mesh);
#ifdef HAVE_PETSC
    // Resolution
    feTolerances tol{1e-9, 1e-9, 20};

    StationarySolver *solver =
      new StationarySolver(tol, metaNumber, linearSystem, sol, norms, mesh);
    solver->makeSteps(0, fespace);
    
    // maxNormL2DC3_U[2 * iter] = fmax(maxNormL2DC3_U[2 * iter], normL2DC3[3 * i]);
    // maxNormL2DC3_V[2 * iter] = fmax(maxNormL2DC3_V[2 * iter], normL2DC3[3 * i + 1]);
    // maxNormL2DC3_P[2 * iter] = fmax(maxNormL2DC3_P[2 * iter], normL2DC3[3 * i + 2]);
    // std::string vtkFile = "../../data/solutionManufacturee" + std::to_string(iter+1) + ".vtk";
    std::string vtkFile = "../../data/CasTest_NS_steady" + std::to_string(iter+1) + ".vtk";
    feExporterVTK writer(vtkFile, mesh, sol, metaNumber, fespace);

    // normU->computeL2Norm(metaNumber, sol, mesh);
    // normV->computeL2Norm(metaNumber, sol, mesh);
    // normP->computeL2Norm(metaNumber, sol, mesh);
    normL2_U[2 * iter] = normU->getNorm();
    normL2_V[2 * iter] = normV->getNorm();
    normL2_Vec[2 * iter] = normVecVeloc->getNorm();
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

  // Calcul du taux de convergence
  for(int i = 1; i < nIter; ++i) {
    normL2_U[2 * i + 1] =
      -log(normL2_U[2 * i] / normL2_U[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
    normL2_V[2 * i + 1] =
      -log(normL2_V[2 * i] / normL2_V[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
    normL2_Vec[2 * i + 1] =
      -log(normL2_Vec[2 * i] / normL2_Vec[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));  
    normL2_P[2 * i + 1] =
      -log(normL2_P[2 * i] / normL2_P[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
  }
  printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s\n", "nElm", "||u-uh||", "tauxU",
         "||v-vh||", "tauxV", "||Uvec-Uh||", "tauxVec",  "||p-ph||", "tauxP");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n", nElm[i],
           normL2_U[2 * i], normL2_U[2 * i + 1], normL2_V[2 * i], normL2_V[2 * i + 1], normL2_Vec[2 * i], normL2_Vec[2 * i + 1],
           normL2_P[2 * i], normL2_P[2 * i + 1]);

  // printf("%12s \t %12s \t %12s\t %12s \t %12s\t %12s \t %12s\n", "nElm", "||E-DC3_U||",
  //        "Taux DC3_U", "||E-DC3_V||", "Taux DC3_V", "||E-DC3_P||", "Taux DC3_P");
  // for(int i = 0; i < nIter; ++i)
  //   printf("%12d \t %12.6e \t %12.6e\t %12.6e \t %12.6e\t %12.6e \t %12.6e\n", nElm[i],
  //          maxNormL2DC3_U[2 * i], maxNormL2DC3_U[2 * i + 1], maxNormL2DC3_V[2 * i],
  //          maxNormL2DC3_V[2 * i + 1], maxNormL2DC3_P[2 * i], maxNormL2DC3_P[2 * i + 1]);
#ifdef USING_PETSC
  petscFinalize();
#endif
  return 0;
}
