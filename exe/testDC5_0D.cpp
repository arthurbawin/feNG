#include "feAPI.h"
#define USING_PETSC

double fSol(const double t, const std::vector<double> &x, const std::vector<double> par)
{
  return sin(2 * M_PI * t) + 1;
}

double fSource(const double t, const std::vector<double> &x, const std::vector<double> par)
{
  return -2 * M_PI * cos(2 * M_PI * t);
}

// double fSol(const double t, const std::vector<double> &x, const std::vector<double> par) {
//   return pow(t, 6.) +1;
// }

// double fSource(const double t, const std::vector<double> &x, const std::vector<double> par) {

//   return -6. * pow(t, 5.);
// }

int main(int argc, char **argv)
{
#ifdef USING_PETSC
  petscInitialize(argc, argv);
#endif
  double xa = 0.;

  feFunction *funSol = new feFunction(fSol, {});
  feFunction *funSource = new feFunction(fSource, {});

  int nIter = 5;
  std::vector<double> normL2_BDF2(2 * nIter, 0.0);
  std::vector<double> normL2_DC3(2 * nIter, 0.0);
  std::vector<double> normL2_DC4(2 * nIter, 0.0);
  std::vector<double> normL2_DC5(2 * nIter, 0.0);
  std::vector<int> nElm(nIter, 0);
  std::vector<int> TT;
  TT.resize(nIter);

  // Boucle sur le nombre de pas de temps
  for(int iter = 0; iter < nIter; ++iter) {
    nElm[iter] = 1;

    // Maillage
    feMesh0DP0 mesh(xa, 1, "M1D");

    // Parameters
    int verbosity = 2;
    int order = 1;
    int degreeQuadrature = 10;

    // Espaces d'interpolation
    int dim;
    feSpace *U_M1D;
    feCheck(createFiniteElementSpace(U_M1D, &mesh, dim = 0, POINT, order, "U", "M1D",
                                     degreeQuadrature, funSol));
    std::vector<feSpace *> spaces = {U_M1D};
    std::vector<feSpace *> essentialSpaces = {};

    // Create a "meta" structure containing the numbering of each unknown field.
    feMetaNumber metaNumber(&mesh, spaces, essentialSpaces);

    // Create the solution structure (contains essentially the solution at all DOFs).
    feSolution sol(&mesh, spaces, essentialSpaces, &metaNumber);

    // Formes (bi)lineaires
    std::vector<feSpace *> spaceSource1D_U = {U_M1D};
    std::vector<feSpace *> spaceMasse1D_U = {U_M1D};
    feBilinearForm source_U_M1D(spaceSource1D_U, &mesh, degreeQuadrature,
                                new feSysElm_0D_Source(1.0, funSource));
    feBilinearForm masse_U_M1D(spaceMasse1D_U, &mesh, degreeQuadrature,
                               new feSysElm_0D_Masse(1.0, nullptr));

    // Initialize the linear system. Assembly of the elementary matrices and RHS is
    // performed in the solve step. Two linear solvers are available :
    // MKL Pardiso (direct solver) and PETSc (collection of iterative solvers).
    feLinearSystem *system;
    feCheck(createLinearSystem(system, PETSC, {&masse_U_M1D, &source_U_M1D}, &metaNumber, &mesh,
                               argc, argv));

    // Norme de la solution
    // Option available ;
    //    - L2Norm_1Field
    //    - L2Norm_MultipleFields
    //    - H1Norm_MultipleFields
    //    - SemiH1Norm_MultipleFields
    //    - Int
    //    - IntNormalx
    //    - IntNormaly
    //    - NormLambda
    //    - NormLambdaNS
    feComputer *normBDF2 = new feComputer(U_M1D, &mesh, &metaNumber, "L2Norm_1Field", funSol);
    feComputer *normDC3 = new feComputer(U_M1D, &mesh, &metaNumber, "L2Norm_1Field", funSol);
    feComputer *normDC4 = new feComputer(U_M1D, &mesh, &metaNumber, "L2Norm_1Field", funSol);
    feComputer *normDC5 = new feComputer(U_M1D, &mesh, &metaNumber, "L2Norm_1Field", funSol);
    std::vector<feComputer *> comput = {normBDF2, normDC3, normDC4, normDC5};

    // Create an exporter structure to write the solution for visualization. Currently the only
    // supported format for visualization is VTK.
    feExporter *exporter;
    feCheck(createVisualizationExporter(exporter, VTK, &metaNumber, &sol, &mesh, spaces));
    int exportEveryNSteps = 1;
    std::string vtkFileRoot = "../../data/ResultatVK/VonKarman";
    feExportData exportData = {exporter, exportEveryNSteps, vtkFileRoot};

    // Solve the discrete problem. Initialize a TimeIntegrator object and tolerances on the
    // Newton-Raphson nonlinear solver (tolerance on the solution correction dx, tolerance on the
    // residual, max number of iterations). Here the PDE is linear : the nonlinear solver should
    // converge in 2 iterations. The TimeIntegrator can be STATIONARY, BDF1, BDF2 or a
    // deferred-correction method (DC2F, DC3, DC3F). The solution will be exported for visualization
    // according to the exportData structure. The linear system is assembled and solved in the
    // "makeSteps()" call.
    TimeIntegrator *solver;
    double t0 = 0.;
    double t1 = 1.;
    int nTimeSteps = 5 * pow(2, iter);
    TT[iter] = nTimeSteps;
    feTolerances tol{1e-10, 1e-10, 10};
    feCheck(createTimeIntegrator(solver, DC5, tol, system, &metaNumber, &sol, &mesh, comput,
                                 exportData, t0, t1, nTimeSteps));
    feCheck(solver->makeSteps(nTimeSteps));
    // Post process the results
    std::vector<double> &normL2BDF2 = normBDF2->getResult();
    std::vector<double> &normL2DC3 = normDC3->getResult();
    std::vector<double> &normL2DC4 = normDC4->getResult();
    std::vector<double> &normL2DC5 = normDC5->getResult();
    normL2_BDF2[2 * iter] = *std::max_element(normL2BDF2.begin(), normL2BDF2.end());
    normL2_DC3[2 * iter] = *std::max_element(normL2DC3.begin(), normL2DC3.end());
    normL2_DC4[2 * iter] = *std::max_element(normL2DC4.begin(), normL2DC4.end());
    normL2_DC5[2 * iter] = *std::max_element(normL2DC5.begin(), normL2DC5.end());

    delete solver;
    delete exporter;
    delete U_M1D;
  }
  delete funSource;
  delete funSol;

  // Calcul du taux de convergence
  for(int i = 1; i < nIter; ++i) {
    normL2_BDF2[2 * i + 1] = log(normL2_BDF2[2 * (i - 1)] / normL2_BDF2[2 * i]) / log(2.);
    normL2_DC3[2 * i + 1] = log(normL2_DC3[2 * (i - 1)] / normL2_DC3[2 * i]) / log(2.);
    normL2_DC4[2 * i + 1] = log(normL2_DC4[2 * (i - 1)] / normL2_DC4[2 * i]) / log(2.);
    normL2_DC5[2 * i + 1] = log(normL2_DC5[2 * (i - 1)] / normL2_DC5[2 * i]) / log(2.);
  }
  printf("%12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \t %12s \n",
         "nTimeSteps", "nElm", "||E_BDF2||", "Taux BDF2", "||E_DC3||", "Taux DC3", "||E_DC4||",
         "Taux DC4", "||E_DC5||", "Taux DC5");
  for(int i = 0; i < nIter; ++i)
    printf("%12d \t %12d \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t %12.6e \t "
           "%12.6e \n",
           TT[i], nElm[i], normL2_BDF2[2 * i], normL2_BDF2[2 * i + 1], normL2_DC3[2 * i],
           normL2_DC3[2 * i + 1], normL2_DC4[2 * i], normL2_DC4[2 * i + 1], normL2_DC5[2 * i],
           normL2_DC5[2 * i + 1]);

  // #ifdef USING_PETSC
  //   petscFinalize();
  // #endif
  return 0;
}
