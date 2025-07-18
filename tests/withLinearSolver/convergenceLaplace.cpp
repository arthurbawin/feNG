
#include "feAPI.h"
#include "../tests.h"

#include <gtest/gtest.h>

/*

  Convergence study for the heat equation with homogeneous Dirichlet BC
  on the unit square and with continuous Lagrange FE.

*/

static int my_argc;
static char** my_argv;

namespace {

  double fSol(const feFunctionArguments &args, const std::vector<double> & /* par */)
  {
    double x = args.pos[0];
    double y = args.pos[1];
    return pow(x, 6) + pow(y, 6);
  }

  double fSource(const feFunctionArguments &args, const std::vector<double> &par)
  {
    double x = args.pos[0];
    double y = args.pos[1];
    double k = par[0];
    return + k * 30. * (pow(x, 4) + pow(y, 4));
  }

  double round_to(double value, double precision)
  {
    return std::round(value / precision) * precision;
  }

  int meshConvergence(std::stringstream &resultBuffer, int order, int numMeshes, int degreeQuadrature)
  {
    setVerbose(1);  

    double k = 1.0;
    feFunction funSol(fSol);
    feFunction funSource(fSource, {k});
    feConstantFunction funZero(0);
    feConstantFunction diffusivity(k);

    std::vector<int> nElm(numMeshes, 0);
    std::vector<double> L2errorU(2 * numMeshes, 0.0);

    //
    // Compute the solution on each mesh
    //
    for(int i = 0; i < numMeshes; ++i)
    {
      std::string meshFile = "../../../data/square" + std::to_string(i+1) + ".msh";

      feMesh2DP1 mesh(meshFile);
      nElm[i] = mesh.getNumInteriorElements();

      feSpace *uBord, *u;
      feCheck(createFiniteElementSpace(uBord, &mesh, elementType::LAGRANGE, order, "U", "Bord", degreeQuadrature, &funSol));
      feCheck(createFiniteElementSpace(    u, &mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &funZero));

      std::vector<feSpace *> spaces = {u, uBord};
      std::vector<feSpace *> essentialSpaces = {uBord};

      feMetaNumber numbering(&mesh, spaces, essentialSpaces);
      feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);

      feBilinearForm *diff = nullptr, *source = nullptr;
      feCheck(createBilinearForm(diff, {u}, new feSysElm_Diffusion<2>(&diffusivity)));
      feCheck(createBilinearForm(source, {u}, new feSysElm_Source(&funSource)));

      std::vector<feBilinearForm*> forms = {diff, source};

  #if defined(HAVE_MKL) || defined(HAVE_PETSC)
      feLinearSystem *system;
      #if defined(HAVE_MKL)
      // Decide if the Pardiso linear system is stored on a single process (-1)
      // or distributed on the processes (>= 0)
      int ownershipSplit = -1;
      feCheck(createLinearSystem(system, MKLPARDISO, {diff, source}, &numbering, ownershipSplit));
      #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
      feCheck(createLinearSystem(system, PETSC_MUMPS, {diff, source}, &numbering));
      #else
      feCheck(createLinearSystem(system, PETSC, {diff, source}, &numbering));
      system->setRelativeTol(1e-10);
      system->setAbsoluteTol(1e-20);
      #endif

      feNorm normU(L2_ERROR, {u}, &sol, &funSol);
      std::vector<feNorm*> norms = {&normU};

      TimeIntegrator *solver;
      feNLSolverOptions NLoptions{1e-13, 1e-13, 1e4, 10, 3, 1e-1};
      feExportData exportData = {nullptr, 1, ""};
      feCheck(createTimeIntegrator(solver, timeIntegratorScheme::STATIONARY, NLoptions, system, &sol, &mesh, norms, exportData));
      feCheck(solver->makeSteps(0));

      // feExporter *exporter;
      // feCheck(createVisualizationExporter(exporter, VTK, &numbering, &sol, &mesh, spaces));
      // feCheck(exporter->writeStep("foo.vtk"));

      // Compute L2 error
      feNorm norm(L2_ERROR, {u}, &sol, &funSol);
      L2errorU[2 * i] = round_to(norm.compute(), 1e-14);

      delete solver;
      delete system;
  #endif

      delete diff;
      delete source;
      delete uBord;
      delete u;
    }

    // Compute the convergence rate
    for(int i = 1; i < numMeshes; ++i) {
      L2errorU[2 * i + 1] =
        -log(L2errorU[2 * i] / L2errorU[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
    }

    resultBuffer << "Lagrange elements P" << order << std::endl;
    printf("%12s \t %12s \t %12s \n", "nElm", "||E_U||", "rate");
    resultBuffer
      << std::setw(16) << std::right << "nElm"
      << std::setw(16) << std::right << "L2 error"
      << std::setw(16) << std::right << "rate" << std::endl;
    for(int i = 0; i < numMeshes; ++i) {
      printf("%12d \t %12.6e \t %12.6e\n", nElm[i], L2errorU[2 * i], L2errorU[2 * i + 1]);
      resultBuffer
        << std::scientific
        << std::setw(16) << std::right
        << std::setprecision(6) << nElm[i] << std::setw(16) << std::right
        << std::setprecision(6) << L2errorU[2 * i] << std::setw(16) << std::right
        << std::setprecision(6) << L2errorU[2 * i + 1]
        << std::endl;
    }
    resultBuffer << std::endl;

    return 0;
  }
}

TEST(Convergence, LaplaceSquare)
{
  initialize(my_argc, my_argv);
  
  std::string testRoot = "../../../tests/withLinearSolver/convergenceLaplace";

  std::stringstream resultBuffer;

  int degreeQuadrature = 12;
  meshConvergence(resultBuffer, 1, 4, degreeQuadrature);
  meshConvergence(resultBuffer, 2, 4, degreeQuadrature);
  meshConvergence(resultBuffer, 3, 4, degreeQuadrature);
  meshConvergence(resultBuffer, 4, 4, degreeQuadrature);

  // Compare the resultBuffer with the file testRoot.output
  // If output file does not exist, a temporary file .tmp_output
  // is created, which may be used as the new .output file for this test.
  EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);

  finalize();
}

//
// Custom main to capture argc and argv
//
int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  my_argc = argc;
  my_argv = argv;
  return RUN_ALL_TESTS();
}