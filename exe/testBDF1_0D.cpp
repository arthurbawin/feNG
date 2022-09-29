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
  return pow(t, 5);
}

double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {
  return -5 * pow(t, 4);
}

int main(int argc, char **argv) {

  petscInitialize(argc, argv);

  double xa = 0.;

  feFunction *funSol = new feFunction(fSol, {});
  feFunction *funSource = new feFunction(fSource, {});

  int nIter = 1;

  std::vector<double> normL2BDF1(2 * nIter, 0.0);

  std::vector<int> nElm(nIter, 0);
  std::vector<int> numberTimeStep;
  numberTimeStep.resize(nIter);
  

  int degreeQuadrature = 10;
  for(int iter = 0; iter < nIter; ++iter) {
    nElm[iter] = 1;
    
    // Maillage
    feMesh0DP0 mesh(xa, 1, "M1D");
    
    // Espaces d'interpolation
    int dim;
    feSpace *U_M1D;
    feCheck(createFiniteElementSpace(U_M1D, &mesh, dim=0, POINT, 0, "U", "M1D", degreeQuadrature, funSol));
    
    std::vector<feSpace *> fespace = {U_M1D};
    std::vector<feSpace *> feEssBC = {};
    
    
    // Numerotations
    feMetaNumber metaNumber(&mesh, fespace, feEssBC);
    

    // Solution
    feSolution sol(&mesh, fespace, feEssBC, &metaNumber);

    
    // Formes (bi)lineaires
    feBilinearForm source({U_M1D}, &mesh, degreeQuadrature, new feSysElm_0D_Source(1.0, funSource));
    feBilinearForm masse({U_M1D}, &mesh, degreeQuadrature, new feSysElm_0D_Masse(1.0, nullptr));
        

    // Systeme lineaire
    feLinearSystem *system;
    // feCheck(createLinearSystem(system, MKLPARDISO, spaces, {&NS2D, &weakBC_UCylinder, &weakBC_VCylinder},&metaNumber, &mesh, argc, argv));
    feCheck(createLinearSystem(system, PETSC, fespace, {&source, &masse}, &metaNumber, &mesh, argc, argv));
    
    // Norme de la solution
    feComputer *norm = new feComputer(U_M1D, &mesh, &metaNumber, "L2Norm_1Field", funSol);
    std::vector<feComputer *> comput = {norm};


    //exporter
    feExporter *exporter;
    feCheck(createVisualizationExporter(exporter, VTK, &metaNumber, &sol, &mesh, fespace));
    int exportEveryNSteps = 1;
    std::string vtkFileRoot = "../../data/ResultatVK/VonKarman";
    feExportData exportData = {nullptr, exportEveryNSteps, vtkFileRoot};


    // Resolution
    TimeIntegrator *solver;
    feTolerances tol{1e-8, 1e-8, 10};
    double t0 = 0.;
    double t1 = 1.;
    std::string CodeIni = "";
    int nTimeSteps = 10 * pow(2, iter);
    numberTimeStep[iter] = nTimeSteps;
    
    feCheck(createTimeIntegrator(solver, BDF1, tol, system, &metaNumber, &sol, &mesh, comput, exportData, t0, t1, nTimeSteps, CodeIni));
    feCheck(solver->makeSteps(nTimeSteps));

    std::vector<double> &normL2_norm = norm->getResult();
    normL2BDF1[2*iter] = *std::max_element(normL2_norm.begin(), normL2_norm.end());


    delete U_M1D; 
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

  petscFinalize();

  return 0;
}
