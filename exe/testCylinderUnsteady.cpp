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
#include "feRecovery.h"
#include "feMetric.h"

double fInlet(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double y = pos[1];
  double h = 4.0;
  double tEnd = 2.0;
  double tRamp = fmin(1.0, t / tEnd);

  return 1.0 + 2.0 * y / h * (1.0 - tRamp);
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  return 0.0;
}

int main(int argc, char **argv) {
#ifdef HAVE_PETSC
  petscInitialize(argc, argv);
#endif

  // Diamètre du cylindre et nombre de Reynolds
  double r = 0.25;
  double d = 2 * r;
  double U = 1.0;
  double rho = 1.0;
  double Re = 100.0;
  double mu = rho * U * d / Re;

  std::vector<double> stokesParam = {rho, mu};
  feFunction *funZero = new feFunction(fZero, {});
  feFunction *funInlet = new feFunction(fInlet, {});

  std::string root, meshName, metricMeshName, nextMeshName;

  // meshName = "../../data/cyl1Msh2.msh";
  meshName = "../../data/cylLargeBox1Msh2.msh";
  size_t lastindex = meshName.find_last_of(".");
  root = meshName.substr(0, lastindex);

  feMesh2DP1 *mesh = new feMesh2DP1(meshName, false);
  feMesh2DP1 *otherMesh;
  feMesh2DP1::mapType physicalEntities = mesh->getPhysicalEntitiesDescription();

  feSpace *U_bas, *U_haut, *U_gauche, *U_cylindre;
  feSpace *U_bas2, *U_haut2, *U_gauche2, *U_cylindre2;
  feSpace *V_bas, *V_haut, *V_gauche, *V_cylindre;
  feSpace *V_bas2, *V_haut2, *V_gauche2, *V_cylindre2;

  feSpace *U_surface, *V_surface, *U_surface2, *V_surface2, *P_surface, *P_surface2;

  feMetaNumber *metaNumber, *otherMetaNumber;
  feSolution *sol, *otherSol;

  std::vector<feNorm *> norms = {};

  feTolerances tol{1e-6, 1e-6, 50};

  double t0 = 0.;
  double t1 = 20.0;
  int nIter = 20;
  int nStepsPerAdaptationCycle = 4;
  int nTotalSteps = nIter * nStepsPerAdaptationCycle;

  int dQuad = 9;

  BDF2Solver *solver;

  std::string elemType = "P2P1";

  for(int iter = 0; iter < nIter; ++iter) {
    if(elemType == "P2P1") {
      U_bas = new feSpace1DP2(mesh, "U", "Bas", funZero);
      U_haut = new feSpace1DP2(mesh, "U", "Haut", funZero);
      U_gauche = new feSpace1DP2(mesh, "U", "Gauche", funInlet); // Inlet
      U_cylindre = new feSpace1DP2(mesh, "U", "Cylindre", funZero);
      U_surface = new feSpaceTriP2(mesh, "U", "Surface", funZero);

      V_bas = new feSpace1DP2(mesh, "V", "Bas", funZero);
      V_haut = new feSpace1DP2(mesh, "V", "Haut", funZero);
      V_gauche = new feSpace1DP2(mesh, "V", "Gauche", funZero);
      V_cylindre = new feSpace1DP2(mesh, "V", "Cylindre", funZero);
      V_surface = new feSpaceTriP2(mesh, "V", "Surface", funZero);

      P_surface = new feSpaceTriP1(mesh, "P", "Surface", funZero);
    } else if(elemType == "P3P2") {
      U_bas = new feSpace1DP3(mesh, "U", "Bas", funZero);
      U_haut = new feSpace1DP3(mesh, "U", "Haut", funZero);
      U_gauche = new feSpace1DP3(mesh, "U", "Gauche", funInlet); // Inlet
      U_cylindre = new feSpace1DP3(mesh, "U", "Cylindre", funZero);
      U_surface = new feSpaceTriP3(mesh, "U", "Surface", funZero);

      V_bas = new feSpace1DP3(mesh, "V", "Bas", funZero);
      V_haut = new feSpace1DP3(mesh, "V", "Haut", funZero);
      V_gauche = new feSpace1DP3(mesh, "V", "Gauche", funZero);
      V_cylindre = new feSpace1DP3(mesh, "V", "Cylindre", funZero);
      V_surface = new feSpaceTriP3(mesh, "V", "Surface", funZero);

      P_surface = new feSpaceTriP2(mesh, "P", "Surface", funZero);
    } else {
      printf("Unknown finite element type : %s\n", elemType.c_str());
    }

    std::vector<feSpace *> fespace = {U_bas,  U_haut,   U_gauche,   U_cylindre, U_surface, V_bas,
                                      V_haut, V_gauche, V_cylindre, V_surface,  P_surface};

    std::vector<feSpace *> feEssBC = {U_gauche, U_cylindre, V_bas, V_haut, V_gauche, V_cylindre};

    metaNumber = new feMetaNumber(mesh, fespace, feEssBC);
    sol = new feSolution(mesh, fespace, feEssBC, metaNumber);

    std::vector<feSpace *> spacesNS2D = {U_surface, V_surface, P_surface};
    feBilinearForm *NS2D = new feBilinearForm(spacesNS2D, mesh, dQuad,
                                              new feSysElm_2D_NavierStokes(stokesParam, nullptr));
    std::vector<feBilinearForm *> formMatrices = {NS2D};
    std::vector<feBilinearForm *> formResiduals = {NS2D};

    feLinearSystemPETSc *linearSystem =
      new feLinearSystemPETSc(argc, argv, formMatrices, formResiduals, metaNumber, mesh);

    if(iter == 0) {
      // Initialize a new BDF2 integrator
      solver = new BDF2Solver(tol, metaNumber, linearSystem, sol, norms, mesh, t0, t1, nTotalSteps);
    } else {
      // Update the existing integrator
      sol->setSolFromContainer(solver->getSolutionContainer(), 0);

      feNorm *newNorm = new feNorm(U_surface, mesh, dQuad, nullptr);
      newNorm->computeL2Norm(metaNumber, sol, mesh);
      printf("New initial norm of U is %4.4f\n", newNorm->getNorm());
      delete newNorm;

      solver->setMetaNumber(metaNumber);
      solver->setLinearSystem(linearSystem);
      solver->setSolution(sol);
      solver->setMesh(mesh);
    }

    printf("===================================================================================\n");
    printf("                                    Iter %d :                                      \n",
           iter);
    printf("  nDOFs = %10d (%10d inc. + %10d ess.)                                      \n",
           metaNumber->getNbDOFs(), metaNumber->getNbUnknowns(),
           metaNumber->getNbDOFs() - metaNumber->getNbUnknowns());
    printf("  sol.size() = %10d \n", sol->getSolutionReference().size());
    printf("  %1d sols in solutionContainer (%10d dofs) :\n",
           solver->getSolutionContainer()->getNbSol(), solver->getSolutionContainer()->getNbDOFs());
    for(int i = 0; i < solver->getSolutionContainer()->getNbSol(); ++i)
      printf("  sol[%d].size = %10d :\n", i, solver->getSolutionContainer()->getSolution(i).size());
    printf("===================================================================================\n");

    solver->makeSteps(nStepsPerAdaptationCycle, fespace);

    std::vector<double> estErreur(2, 0.);
    feRecovery *recU = new feRecovery(metaNumber, U_surface, mesh, sol, estErreur, funZero);
    // feRecovery *recV = new feRecovery(metaNumber, V_surface, mesh, sol, estErreur, funZero);
    // feRecovery *recP = new feRecovery(metaNumber, P_surface, mesh, sol, estErreur, funZero);

    feMetricOptions metricOptions;
    metricOptions.nTargetVertices = 2500;
    metricOptions.hMin = 0.000001;
    metricOptions.hMax = 100;
    metricOptions.computationMethod = 1;
    metricOptions.LpNorm = 2.0;
    metricOptions.nPhi = 801;

    metricMeshName = root + "SizeField" + std::to_string(iter + 1) + ".msh";
    nextMeshName = root + "Adapted" + std::to_string(iter + 1) + ".msh";

    feMetric *metric = new feMetric(recU, metricOptions);
    metric->writeSizeFieldGmsh(meshName, metricMeshName);

    std::string cmdMMGS = "mmg2d " + metricMeshName + " -o " + nextMeshName + " -hgrad 100";
    std::string cmdGMSH = "gmsh " + nextMeshName + " &";

    system(cmdMMGS.c_str());
    system(cmdGMSH.c_str());

    std::string vtkFile = "transferFrom" + std::to_string(iter + 1) + ".vtk";
    feExporterVTK writerFrom(vtkFile, mesh, sol, metaNumber, fespace);

    std::cout << "Reading " << nextMeshName << std::endl;
    otherMesh = new feMesh2DP1(nextMeshName, false, physicalEntities);

    if(elemType == "P2P1") {
      U_bas2 = new feSpace1DP2(otherMesh, "U", "Bas", funZero);
      U_haut2 = new feSpace1DP2(otherMesh, "U", "Haut", funZero);
      U_gauche2 = new feSpace1DP2(otherMesh, "U", "Gauche", funInlet); // Inlet
      U_cylindre2 = new feSpace1DP2(otherMesh, "U", "Cylindre", funZero);
      U_surface2 = new feSpaceTriP2(otherMesh, "U", "Surface", funZero);

      V_bas2 = new feSpace1DP2(otherMesh, "V", "Bas", funZero);
      V_haut2 = new feSpace1DP2(otherMesh, "V", "Haut", funZero);
      V_gauche2 = new feSpace1DP2(otherMesh, "V", "Gauche", funZero);
      V_cylindre2 = new feSpace1DP2(otherMesh, "V", "Cylindre", funZero);
      V_surface2 = new feSpaceTriP2(otherMesh, "V", "Surface", funZero);

      P_surface2 = new feSpaceTriP1(otherMesh, "P", "Surface", funZero);
    } else if(elemType == "P3P2") {
      U_bas2 = new feSpace1DP3(otherMesh, "U", "Bas", funZero);
      U_haut2 = new feSpace1DP3(otherMesh, "U", "Haut", funZero);
      U_gauche2 = new feSpace1DP3(otherMesh, "U", "Gauche", funInlet); // Inlet
      U_cylindre2 = new feSpace1DP3(otherMesh, "U", "Cylindre", funZero);
      U_surface2 = new feSpaceTriP3(otherMesh, "U", "Surface", funZero);

      V_bas2 = new feSpace1DP3(otherMesh, "V", "Bas", funZero);
      V_haut2 = new feSpace1DP3(otherMesh, "V", "Haut", funZero);
      V_gauche2 = new feSpace1DP3(otherMesh, "V", "Gauche", funZero);
      V_cylindre2 = new feSpace1DP3(otherMesh, "V", "Cylindre", funZero);
      V_surface2 = new feSpaceTriP3(otherMesh, "V", "Surface", funZero);

      P_surface2 = new feSpaceTriP2(otherMesh, "P", "Surface", funZero);
    } else {
      printf("Unknown finite element type : %s\n", elemType.c_str());
    }

    std::vector<feSpace *> otherSpaces = {U_bas2,      U_haut2,    U_gauche2, U_cylindre2,
                                          U_surface2,  V_bas2,     V_haut2,   V_gauche2,
                                          V_cylindre2, V_surface2, P_surface2};

    std::vector<feSpace *> otherSpacesEss = {U_gauche2, U_cylindre2, V_bas2,
                                             V_haut2,   V_gauche2,   V_cylindre2};

    otherMetaNumber = new feMetaNumber(otherMesh, otherSpaces, otherSpacesEss);
    otherSol = new feSolution(otherMesh, otherSpaces, otherSpacesEss, otherMetaNumber);

    // Compute norm of U on the current mesh
    feNorm *normBeforeTransfer = new feNorm(U_surface, mesh, dQuad, nullptr);
    normBeforeTransfer->computeL2Norm(metaNumber, sol, mesh);
    printf("Initial norm of U is %4.4f\n", normBeforeTransfer->getNorm());
    delete normBeforeTransfer;

    mesh->transfer(otherMesh, metaNumber, otherMetaNumber, solver->getSolutionContainer(), fespace,
                   feEssBC, otherSpaces);

    otherSol->setSolFromContainer(solver->getSolutionContainer(), 0);

    // Compute norm of U on the new mesh
    feNorm *normAfterTransfer = new feNorm(U_surface2, otherMesh, dQuad, nullptr);
    normAfterTransfer->computeL2Norm(otherMetaNumber, otherSol, otherMesh);
    printf("Final   norm of U is %4.4f\n", normAfterTransfer->getNorm());
    delete normAfterTransfer;

    vtkFile = "transferTo" + std::to_string(iter + 1) + ".vtk";
    feExporterVTK writerTo(vtkFile, otherMesh, otherSol, otherMetaNumber, otherSpaces);

    mesh = otherMesh;
    meshName = nextMeshName;

    delete otherSol;
    delete otherMetaNumber;

    for(feSpace *fS : otherSpaces) delete fS;
    delete metric;
    delete recU;
    delete linearSystem;
    delete NS2D;
    delete sol;
    delete metaNumber;
    for(feSpace *fS : fespace) delete fS;
  }
  delete solver;
  delete mesh;

  delete funZero;
  delete funInlet;

#ifdef HAVE_PETSC
  petscFinalize();
#endif

  return 0;
}
