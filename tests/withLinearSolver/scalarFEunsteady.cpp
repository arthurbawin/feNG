
#include "feAPI.h"
#include "../tests.h"

#include <gtest/gtest.h>

/*
  Tests the implementation of scalar FE spaces and weak forms
*/

static int my_argc;
static char** my_argv;

namespace diffusion {

  double uSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double t = args.t;
    return exp(-t) * sin(M_PI*x)*sin(M_PI*y);
  }

  double uSrc_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double x = args.pos[0];
    const double y = args.pos[1];
    const double t = args.t;
    const double dudt = - exp(-t) * sin(M_PI*x)*sin(M_PI*y);
    const double lapu = - 2. * exp(-t) * M_PI*M_PI * sin(M_PI*x)*sin(M_PI*y);
    return - (dudt - lapu);
  }

  feStatus solve(const std::string &meshFile,
                 const timeIntegratorScheme scheme,
                 const double t0,
                 const double t1,
                 const int nTimeSteps,
                 const int order,
                 const int degreeQuadrature,
                 const feNLSolverOptions &NLoptions,
                 int &numInteriorElements,
                 double &L2Error)
  {
    feStatus s;
    feFunction uSol(uSol_f);
    feFunction uSrc(uSrc_f);

    feMesh2DP1 mesh(meshFile);
    numInteriorElements = mesh.getNumInteriorElements();

    feSpace *u = nullptr, *uBord = nullptr;
    s = createFiniteElementSpace(uBord, &mesh, elementType::LAGRANGE, order, "U",    "Bord", degreeQuadrature, &uSol);
    if(s != FE_STATUS_OK) return s;
    s = createFiniteElementSpace(    u, &mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &scalarConstant::zero);
    if(s != FE_STATUS_OK) return s;

    std::vector<feSpace*> spaces = {uBord, u};
    std::vector<feSpace*> essentialSpaces = {uBord};

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
    
    feBilinearForm *dudt = nullptr, *diff = nullptr, *src = nullptr;
    s = createBilinearForm(dudt, {u}, new feSysElm_TransientMass(&scalarConstant::one));
    if(s != FE_STATUS_OK) return s;
    s = createBilinearForm(diff, {u}, new feSysElm_Diffusion<2>(&scalarConstant::one));
    if(s != FE_STATUS_OK) return s;
    s = createBilinearForm( src, {u}, new feSysElm_Source(&uSrc));
    if(s != FE_STATUS_OK) return s;
    std::vector<feBilinearForm*> forms = {dudt, diff, src};

    feLinearSystem *system;
    #if defined(HAVE_MKL)
      s = createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns());
    #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
      s = createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns());
    #else
      s = createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns());
    #endif
    if(s != FE_STATUS_OK) return s;

    feNorm *errorU_L2 = nullptr;
    s = createNorm(errorU_L2, L2_ERROR, {u}, &sol, &uSol);
    if(s != FE_STATUS_OK) return s;
    std::vector<feNorm *> norms = {};

    TimeIntegrator *solver;
    s = createTimeIntegrator(solver, scheme, NLoptions,
      system, &sol, &mesh, norms, {nullptr, 1, ""}, t0, t1, nTimeSteps);
    if(s != FE_STATUS_OK) return s;
    s = solver->makeSteps(nTimeSteps);
    if(s != FE_STATUS_OK) return s;

    L2Error = errorU_L2->compute();

    delete errorU_L2;
    delete solver;
    delete system;
    for(feBilinearForm* f : forms)
      delete f;
    for(feSpace *space : spaces)
      delete space;

    return FE_STATUS_OK;
  }

  feStatus meshConvergence(std::stringstream &resultBuffer,
    timeIntegratorScheme scheme, double t0, double t1, 
    int initial_nTimeSteps, int order, int numMeshes, int degreeQuadrature)
  {
    feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 10, 4, 1e-1};
    {
      std::vector<int> nElm(numMeshes);
      std::vector<double> err(numMeshes, 0.);
      int nTimeSteps = initial_nTimeSteps;

      for(int i = 0; i < numMeshes; ++i, nTimeSteps *= 2)
      {
        std::string meshFile = "../../../data/square" + std::to_string(i+1) + ".msh";
        feStatus s = solve(meshFile, scheme, t0, t1, nTimeSteps, order, degreeQuadrature, NLoptions, nElm[i], err[i]);
        if(s != FE_STATUS_OK) return s;
      }
      resultBuffer << "Unsteady diffusion - Error on u - Lagrange elements P" << order << std::endl;
      computeAndPrintConvergence(2, numMeshes, err, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
    }
    return FE_STATUS_OK;
  }
}

TEST(ScalarFEUnsteady, DiffusionBDF1)
{
  initialize(my_argc, my_argv);
  setVerbose(0);
  std::string testRoot = "../../../tests/withLinearSolver/scalarFEunsteady_bdf1";
  std::stringstream resultBuffer;
  int degreeQuadrature = 8;
  double t0 = 0.;
  double t1 = 1.;
  int initial_nTimeSteps = 10;
  timeIntegratorScheme scheme = timeIntegratorScheme::BDF1;
  ASSERT_TRUE(diffusion::meshConvergence(resultBuffer, scheme, t0, t1, initial_nTimeSteps, 1, 4, degreeQuadrature) == FE_STATUS_OK);

  EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
  finalize();
}

TEST(ScalarFEUnsteady, DiffusionBDF2)
{
  initialize(my_argc, my_argv);
  setVerbose(0);
  std::string testRoot = "../../../tests/withLinearSolver/scalarFEunsteady_bdf2";
  std::stringstream resultBuffer;
  int degreeQuadrature = 8;
  double t0 = 0.;
  double t1 = 1.;
  int initial_nTimeSteps = 10;
  timeIntegratorScheme scheme = timeIntegratorScheme::BDF2;
  ASSERT_TRUE(diffusion::meshConvergence(resultBuffer, scheme, t0, t1, initial_nTimeSteps, 1, 4, degreeQuadrature) == FE_STATUS_OK);
  // initial_nTimeSteps = 20;
  // ASSERT_TRUE(diffusion::meshConvergence(resultBuffer, scheme, t0, t1, initial_nTimeSteps, 2, 4, degreeQuadrature) == FE_STATUS_OK);

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