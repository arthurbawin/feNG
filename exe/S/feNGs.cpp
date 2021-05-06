#include <iostream>
#include <functional>
#include <vector>

#include "feMesh.h"
#include "feQuadrature.h"
#include "feSpace.h"
#include "feCncGeo.h"
#include "feNumber.h"
#include "feSolution.h"
#include "feSysElm.h"
#include "feBilinearForm.h"
#include "feLinearSystem.h"
#include "feSolver.h"

double f1(const double t, const std::vector<double> x){
	return x[0];
}

double f2(const double t, const std::vector<double> x){
	return 0.0;
}

void test1(int argc, char** argv){

	double xa = 1.;
	double xb = 4.;
	int nElm = 10;

	feMesh1DP1 *mesh = new feMesh1DP1(xa,xb,nElm,"BXA","BXB","M1D");

  feSpace1DP0 U_BXA = feSpace1DP0(mesh, "U", "BXA", f1);
  feSpace1DP0 U_BXB = feSpace1DP0(mesh, "U", "BXB", f1);
  feSpace1DP3 U_M1D = feSpace1DP3(mesh, "U", "M1D", f1);
  feSpace1DP0 V_BXA = feSpace1DP0(mesh, "V", "BXA", f1);
  feSpace1DP0 V_BXB = feSpace1DP0(mesh, "V", "BXB", f1);
  feSpace1DP3 V_M1D = feSpace1DP3(mesh, "V", "M1D", f1);
  std::vector<feSpace*> fespace = {&U_BXA, &U_BXB, &U_M1D, &V_BXA, &V_BXB, &V_M1D};
  std::vector<feSpace*> feEssBC = {&U_BXA, &U_BXB, &V_BXA, &V_BXB};

	// std::vector<feSpace*> fespace(6);
	// fespace[0] = new feSpace1DP0(mesh, "U", "BXA", f1);
	// fespace[1] = new feSpace1DP0(mesh, "U", "BXB", f1);
	// fespace[2] = new feSpace1DP3(mesh, "U", "M1D", f1);
	// fespace[3] = new feSpace1DP0(mesh, "V", "BXA", f1);
	// fespace[4] = new feSpace1DP0(mesh, "V", "BXB", f1);
	// fespace[5] = new feSpace1DP3(mesh, "V", "M1D", f1);

	// std::vector<feSpace*> feEssBC(4);
	// feEssBC[0] = new feSpace1DP0(mesh, "U", "BXA", f1);
	// feEssBC[1] = new feSpace1DP0(mesh, "U", "BXB", f1);
	// feEssBC[2] = new feSpace1DP0(mesh, "V", "BXA", f1);
	// feEssBC[3] = new feSpace1DP0(mesh, "V", "BXB", f1);

	feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC); // Le déjà célèbre META_NUMER

	feSolution *sol = new feSolution(mesh, fespace, feEssBC, metaNumber);
	sol->initializeUnknowns(mesh, metaNumber);
	sol->initializeEssentialBC(mesh, metaNumber);

  // feSysElm_1D_Source *source1D = new feSysElm_1D_Source(0.0, f_source);
  // std::vector<feSpace*> spaceSource1D(1, fespace[2]);
  // feBilinearForm *b = new feBilinearForm(spaceSource1D, mesh, 3, source1D);
  // b->initialize(metaNumber, mesh, sol, 1);
  // b->computeMatrix(metaNumber, mesh, sol, 0);
  // b->computeRHS(metaNumber, mesh, sol, 0);
  // b->printInfo();

  feSysElm_1D_Diffusion *diffusion1D = new feSysElm_1D_Diffusion(1.0, nullptr);
  std::vector<feSpace*> spaceDiffusion1D_U = {&U_M1D};
  std::vector<feSpace*> spaceDiffusion1D_V = {&V_M1D};

  int nQuadraturePoints = 3; // TODO : change to deg
  feBilinearForm *bDiff_U_M1D = new feBilinearForm(spaceDiffusion1D_U, mesh, nQuadraturePoints, diffusion1D);
  feBilinearForm *bDiff_V_M1D = new feBilinearForm(spaceDiffusion1D_V, mesh, nQuadraturePoints, diffusion1D);

  std::vector<feBilinearForm*> formMatrices  = {bDiff_U_M1D, bDiff_V_M1D};
  std::vector<feBilinearForm*> formResiduals = {bDiff_U_M1D, bDiff_V_M1D};
  feLinearSystem *linearSystem = new feLinearSystem(formMatrices, formResiduals, metaNumber, mesh);
  linearSystem->initialize(argc,argv);
  linearSystem->assembleMatrices(sol);
  linearSystem->assembleResiduals(sol);
  double normDx, normRes;
  linearSystem->solve(&normDx, &normRes);
  linearSystem->check(); // Prints the solution
  linearSystem->finalize();

  delete linearSystem;
  delete bDiff_U_M1D;
  delete bDiff_V_M1D;
  delete diffusion1D;
  // delete b;
  // delete source1D;
	delete sol;
	delete metaNumber;

	// for(auto fS : feEssBC) delete fS;
	// for(auto fS : fespace) delete fS;
	delete mesh;
}

double kd = 0.1;

double fSol(const double t, const std::vector<double> x){
  return pow(x[0],6);
}

double fSource(const double t, const std::vector<double> x){
  return kd*30.*pow(x[0],4);
}

void testStationnaire(int argc, char** argv){
  double xa = 0.;
  double xb = 5.;
  int nElm = 20;

  // Maillage
  feMesh1DP1 *mesh = new feMesh1DP1(xa,xb,nElm,"BXA","BXB","M1D");
  // Espaces d'interpolation
  feSpace1DP0 U_BXA = feSpace1DP0(mesh, "U", "BXA", fSol);
  feSpace1DP0 U_BXB = feSpace1DP0(mesh, "U", "BXB", fSol);
  feSpace1DP3 U_M1D = feSpace1DP3(mesh, "U", "M1D", fSol);
  std::vector<feSpace*> fespace = {&U_BXA, &U_BXB, &U_M1D};
  std::vector<feSpace*> feEssBC = {&U_BXA, &U_BXB};
  // Numerotations
  feMetaNumber *metaNumber = new feMetaNumber(mesh, fespace, feEssBC);
  // Solution
  feSolution *sol = new feSolution(mesh, fespace, feEssBC, metaNumber);
  sol->initializeUnknowns(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber);
  // Formes (bi)lineaires
  int nQuad = 3; // TODO : change to deg
  std::vector<feSpace*> spaceDiffusion1D_U = {&U_M1D};
  feBilinearForm *diff_U_M1D = new feBilinearForm(spaceDiffusion1D_U, mesh, nQuad, new feSysElm_1D_Diffusion(kd, nullptr));
  std::vector<feSpace*> spaceSource1D_U = {&U_M1D};
  feBilinearForm *source_U_M1D = new feBilinearForm(spaceSource1D_U, mesh, nQuad, new feSysElm_1D_Source(1.0, fSource));

  std::vector<feBilinearForm*> formMatrices  = {diff_U_M1D};
  std::vector<feBilinearForm*> formResiduals  = {diff_U_M1D, source_U_M1D};

  feLinearSystem *linearSystem = new feLinearSystem(formMatrices, formResiduals, metaNumber, mesh);
  linearSystem->initialize(argc,argv);

  feTolerances tol{1e-9, 1e-8, 10};
  feNorm *norm = new feNorm();
  solveStationary(tol, metaNumber, linearSystem, formMatrices, formResiduals, sol, norm, mesh);

  // TODO : tester la convergence

  delete norm;
  delete linearSystem;
  delete source_U_M1D;
  delete diff_U_M1D;
  delete sol;
  delete metaNumber;
  delete mesh;
}

int main(int argc, char** argv){
	// test1(argc,argv);
  testStationnaire(argc,argv);
  return 0;
}
