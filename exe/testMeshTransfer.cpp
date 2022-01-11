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
  return pow(x[0], 6);
  // return pow(x[0],3);
  // return pow(x[0],2);
}

double fSource(const double t, const std::vector<double> x, const std::vector<double> par) {
  double kd = par[0];
  return kd * 30. * pow(x[0], 4);
  // return kd*6*x[0];
  // return kd*2;
}

double fZero(const double t, const std::vector<double> x, const std::vector<double> par) {
  return 0.;
}

int main(int argc, char **argv) {
  double kd = 1.0;
  feFunction *funSol = new feFunction(fSol, {kd});
  feFunction *funSource = new feFunction(fSource, {kd});
  feFunction *funZero = new feFunction(fZero, {});

  int nIter = 1;
  std::vector<double> normL2(2 * nIter, 0.0);
  std::vector<double> normL2_U(2 * nIter, 0.0);
  std::vector<double> normL2_V(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);

  for(int iter = 0; iter < nIter; ++iter) {
    std::string meshName = "../../data/square" + std::to_string(iter + 1) + ".msh";
    // std::string meshName = "../../data/squareCoarse.msh";

    feMesh2DP1 *mesh = new feMesh2DP1(meshName, false);
    nElm[iter] = mesh->getNbInteriorElems();
    // mesh->printInfo();

    feSpace1DP1 U_angle(mesh, "U", "Angle", funSol);
    feSpace1DP1 U_haut(mesh, "U", "Haut", funSol);
    feSpace1DP1 U_gauche(mesh, "U", "Gauche", funSol);
    feSpaceTriP1 U_surface(mesh, "U", "Surface", funSol);
    feSpace1DP1 V_angle(mesh, "V", "Angle", funSol);
    feSpace1DP1 V_haut(mesh, "V", "Haut", funSol);
    feSpace1DP1 V_gauche(mesh, "V", "Gauche", funSol);
    feSpaceTriP1 V_surface(mesh, "V", "Surface", funSol);

    std::vector<feSpace *> fespace = {&U_angle, &U_haut, &U_gauche, &U_surface,
                                      &V_angle, &V_haut, &V_gauche, &V_surface};
    std::vector<feSpace *> feEssBC = {&U_angle, &U_haut, &U_gauche, &V_angle, &V_haut, &V_gauche};

    feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);

    feSolution *sol = new feSolution(mesh, fespace, feEssBC, metaNumber);

    // Formes (bi)lineaires
    int nQuad = 9; // TODO : change to deg
    std::vector<feSpace *> spaceDiffusion2D_U = {&U_surface};
    std::vector<feSpace *> spaceDiffusion2D_V = {&V_surface};
    std::vector<feSpace *> spaceSource2D_U = {&U_surface};
    std::vector<feSpace *> spaceSource2D_V = {&V_surface};

    feBilinearForm *diffU =
      new feBilinearForm(spaceDiffusion2D_U, mesh, nQuad, new feSysElm_2D_Diffusion(kd, nullptr));
    feBilinearForm *diffV =
      new feBilinearForm(spaceDiffusion2D_V, mesh, nQuad, new feSysElm_2D_Diffusion(kd, nullptr));
    feBilinearForm *sourceU =
      new feBilinearForm(spaceSource2D_U, mesh, nQuad, new feSysElm_2D_Source(1.0, funSource));
    feBilinearForm *sourceV =
      new feBilinearForm(spaceSource2D_V, mesh, nQuad, new feSysElm_2D_Source(1.0, funSource));

    std::vector<feBilinearForm *> formMatrices = {diffU, diffV};
    std::vector<feBilinearForm *> formResiduals = {diffU, diffV, sourceU, sourceV};

    feNorm *normU = new feNorm({&U_surface}, mesh, nQuad, funSol);
    feNorm *normV = new feNorm({&V_surface}, mesh, nQuad, funSol);
    std::vector<feNorm *> norms = {normU, normV};

    feLinearSystemPETSc *linearSystem =
      new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);
    linearSystem->initialize();

    feTolerances tol{1e-9, 1e-8, 3};
    StationarySolver solver(tol, metaNumber, linearSystem, sol, norms, mesh);
    solver.makeSteps(0, fespace);
    normL2[2 * iter] = solver.getNorm(0);
    // solveStationary(&normL2[2*iter], tol, metaNumber, linearSystem, formMatrices, formResiduals,
    // sol, norms, mesh);
    linearSystem->finalize();
    // normL2_U[2*iter] = normU->getNorm();
    // normL2_V[2*iter] = normV->getNorm();

    std::string vtkFile = "transferFrom.vtk";
    feExporterVTK writerFrom(vtkFile, mesh, sol, metaNumber, fespace);

    std::string otherMeshName = "../../data/square" + std::to_string(iter + 3) + ".msh";
    feMesh2DP1 *otherMesh = new feMesh2DP1(otherMeshName, false);
    feSpaceTriP1 U_surface2(otherMesh, "U", "Surface", funZero);
    feSpaceTriP1 V_surface2(otherMesh, "V", "Surface", funZero);
    std::vector<feSpace *> otherSpaces = {&U_surface2, &V_surface2};
    std::vector<feSpace *> otherSpacesEss = {};

    feMetaNumber *otherMetaNumber = new feMetaNumber(otherMesh, otherSpaces, otherSpacesEss);
    feSolution *otherSol = new feSolution(otherMesh, otherSpaces, otherSpacesEss, otherMetaNumber);

    mesh->transfer(otherMesh, metaNumber, otherMetaNumber, solver.getSolutionContainer(), fespace,
                   feEssBC, otherSpaces);

    otherSol->setSolFromContainer(solver.getSolutionContainer(), 0);

    vtkFile = "transferTo.vtk";
    feExporterVTK writerTo(vtkFile, otherMesh, otherSol, otherMetaNumber, otherSpaces);

    delete normU;
    delete normV;
    delete linearSystem;
    delete sourceU;
    delete diffU;
    delete sol;
    delete metaNumber;
    delete mesh;
  }
  delete funSource;
  delete funSol;
  delete funZero;

  // Calcul du taux de convergence
  for(int i = 1; i < nIter; ++i) {
    normL2_U[2 * i + 1] =
      -log(normL2_U[2 * i] / normL2_U[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
    normL2_V[2 * i + 1] =
      -log(normL2_V[2 * i] / normL2_V[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
  }
  printf("%12s \t %12s \t %12s \t %12s \t %12s\n", "nElm", "||E_U||", "p_U", "||E_V||", "p_V");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e\n", nElm[i], normL2_U[2 * i],
           normL2_U[2 * i + 1], normL2_V[2 * i], normL2_V[2 * i + 1]);

  return 0;
}
