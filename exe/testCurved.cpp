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
#include "feLinearSystemPETSc.h"
#include "feExporter.h"

double fSol(const double t, const std::vector<double> x, const std::vector<double> par) {
  // return pow(x[0],6);
  return pow(x[0], 3);
}

double fSource(const double t, const std::vector<double> x, const std::vector<double> par) {
  double kd = par[0];
  return kd * 30. * pow(x[0], 4);
  // return kd*6*x[0];
}

int main(int argc, char **argv) {
  double kd = 1.0;
  feFunction *funSol = new feFunction(fSol, {kd});
  feFunction *funSource = new feFunction(fSource, {kd});

  int nIter = 1;
  std::vector<double> normL2(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);

  for(int iter = 0; iter < nIter; ++iter) {
    std::string meshName = "../../data/curved4.msh";
    // std::string meshName = "../../data/squareCoarseP2.msh";

    feMesh2DP1 *mesh = new feMesh2DP1(meshName, true);
    nElm[iter] = mesh->getNbInteriorElems();
    mesh->printInfo();

    feSpace1DP2 U_bord = feSpace1DP2(mesh, "U", "Bord", funSol);
    // feSpace1DP2  U_angle   = feSpace1DP2(mesh, "U", "Angle",  funSol);
    // feSpace1DP2  U_haut    = feSpace1DP2(mesh, "U", "Haut",   funSol);
    // feSpace1DP2  U_gauche  = feSpace1DP2(mesh, "U", "Gauche", funSol);
    feSpaceTriP2 U_surface = feSpaceTriP2(mesh, "U", "Surface", funSol);

    // std::vector<feSpace*> fespace = {&U_angle, &U_haut, &U_gauche, &U_surface};
    std::vector<feSpace *> fespace = {&U_bord, &U_surface};
    // std::vector<feSpace*> feEssBC = {&U_angle, &U_haut, &U_gauche};
    std::vector<feSpace *> feEssBC = {&U_bord};

    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);
    // // metaNumber->printNumberings();

    feSolution *sol = new feSolution(mesh, fespace, feEssBC, metaNumber);
    sol->initializeUnknowns(mesh, metaNumber);

    // // Formes (bi)lineaires
    int nQuad = 16; // TODO : change to deg
    std::vector<feSpace *> spaceDiffusion2D_U = {&U_surface};
    feBilinearForm *diffU =
      new feBilinearForm(spaceDiffusion2D_U, mesh, nQuad, new feSysElm_2D_Diffusion(kd, nullptr));
    std::vector<feSpace *> spaceSource2D_U = {&U_surface};
    feBilinearForm *sourceU =
      new feBilinearForm(spaceSource2D_U, mesh, nQuad, new feSysElm_2D_Source(1.0, funSource));

    std::vector<feBilinearForm *> formMatrices = {diffU};
    std::vector<feBilinearForm *> formResiduals = {diffU, sourceU};

    feNorm *norm = new feNorm(&U_surface, mesh, nQuad, funSol);
    std::vector<feNorm *> norms = {norm};
    double integral = 0.0;
    norm->computeArea(metaNumber, sol, mesh);
    integral = norm->getNorm();
    std::cout << "Area = " << integral << std::endl;
    norm->computeIntegral(metaNumber, sol, mesh, funSol);
    integral = norm->getNorm();
    std::cout << "Integral = " << integral << std::endl;

    feLinearSystemPETSc *linearSystem =
      new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
    linearSystem->initialize();
    feTolerances tol{1e-9, 1e-8, 20};
    solveStationary(&normL2[2 * iter], tol, metaNumber, linearSystem, formMatrices, formResiduals,
                    sol, norms, mesh);
    linearSystem->finalize();

    // sol->initializeUnknowns(mesh, metaNumber);

    // std::string vtkFile = "../../data/square" + std::to_string(iter+1) + ".vtk";
    // feExporterVTK writer(vtkFile, mesh, sol, metaNumber, fespace);

    delete linearSystem;
    delete norm;
    delete sourceU;
    delete diffU;
    delete sol;
    delete metaNumber;
    delete mesh;
  }
  delete funSource;
  delete funSol;

  // Calcul du taux de convergence
  for(int i = 1; i < nIter; ++i) {
    normL2[2 * i + 1] =
      -log(normL2[2 * i] / normL2[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
  }
  printf("%12s \t %12s \t %12s\n", "nElm", "||E||", "p");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12.6e \t %12.6e\n", nElm[i], normL2[2 * i], normL2[2 * i + 1]);

  return 0;
}
