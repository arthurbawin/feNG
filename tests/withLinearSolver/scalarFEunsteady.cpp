
#include "feAPI.h"
#include "../tests.h"

#include <gtest/gtest.h>

/*
  Tests the implementation of scalar FE spaces and weak forms
*/

static int my_argc;
static char** my_argv;

static std::vector<timeIntegratorScheme> schemes = {timeIntegratorScheme::BDF1,
                                                    timeIntegratorScheme::BDF2,
                                                    timeIntegratorScheme::DC2BDF1,
                                                    timeIntegratorScheme::DC3BDF1};

static std::vector<std::string> schemesStr = {"BDF1", "BDF2", "DC2BDF1", "DC3BDF1"};

//
// Solve the linear ODE : u'(t) = f(t)
//                         u(0) = u0.
struct linearODE {

  static double FT(const double t) { return t*t*t*t; }
  static double DFDT(const double t) { return 4.*t*t*t; }

  // static double FT(const double t) { return sin(t); }
  // static double DFDT(const double t) { return  cos(t); }

  // static double FT(const double t) { return exp(-t); }
  // static double DFDT(const double t) { return -exp(-t); }

  static double uSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double t = args.t;
    return FT(t);
  }

  static double uSolDot_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double t = args.t;
    return DFDT(t);
  }

  static double uSrc_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double t = args.t;                      
    return - DFDT(t);
  }

  static feStatus solve(const std::string &meshFile,
                 const timeIntegratorScheme scheme,
                 const double t0,
                 const double t1,
                 const int nTimeSteps,
                 const int order,
                 const int degreeQuadrature,
                 const feNLSolverOptions &NLoptions,
                 double &error)
  {
    feFunction uSol(uSol_f), uSolDot(uSolDot_f), uSrc(uSrc_f);

    feMesh2DP1 mesh(meshFile);

    feSpace *u = nullptr;
    feCheckReturn(createFiniteElementSpace(u, &mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &uSol));
    std::vector<feSpace*> spaces = {u}, essentialSpaces = {};

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
    
    feBilinearForm *dudt = nullptr, *src = nullptr;
    feCheckReturn(createBilinearForm(dudt, {u}, new feSysElm_TransientMass(&scalarConstant::one)));
    feCheckReturn(createBilinearForm( src, {u}, new feSysElm_Source(&uSrc)));
    std::vector<feBilinearForm*> forms = {dudt, src};

    feLinearSystem *system;
    #if defined(HAVE_MKL)
      feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
    #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
      feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
    #else
      feCheckReturn(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
    #endif

    feNorm *errorU = nullptr, *errorUdot = nullptr;
    feCheckReturn(createNorm(errorU, LINF_ERROR, {u}, &sol, &uSol));
    feCheckReturn(createNorm(errorUdot, LINF_UDOT_ERROR, {u}, &sol, &uSolDot));
    std::vector<feNorm *> norms = {errorU, errorUdot};

    TimeIntegrator *solver;
    feCheckReturn(createTimeIntegrator(solver, scheme, NLoptions,
      system, &sol, &mesh, norms, {nullptr, 1, ""}, t0, t1, nTimeSteps));
    feCheckReturn(solver->makeSteps(nTimeSteps));

    std::vector<std::vector<double>> postProc = solver->getPostProcessingData();
    error = 0.;
    for(int j = 0; j < nTimeSteps+1; ++j) {
      feInfo("Error at step %2d (t = %1.3e) = %1.6e (U) = %1.6e (Udot)", j, postProc[0][j], postProc[1][j], postProc[2][j]);
      error = fmax(error, postProc[1][j]);
    }

    delete errorU;
    delete errorUdot;
    delete solver;
    delete system;
    for(feBilinearForm* f : forms)
      delete f;
    for(feSpace *space : spaces)
      delete space;

    return FE_STATUS_OK;
  }
};

//
// A nonlinear ODE : u'(t) = -u^2
//                    u(0) = 1,
//
// with solution u(t) = 1 / (1+t).
//
struct nonlinearODE {

  static double uSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double t = args.t;
    return 1. / (1. + t);
  }

  static feStatus solve(const std::string &meshFile,
                 const timeIntegratorScheme scheme,
                 const double t0,
                 const double t1,
                 const int nTimeSteps,
                 const int order,
                 const int degreeQuadrature,
                 const feNLSolverOptions &NLoptions,
                 double &error)
  {
    feFunction uSol(uSol_f);

    feMesh2DP1 mesh(meshFile);

    feSpace *u = nullptr;
    feCheckReturn(createFiniteElementSpace(u, &mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &uSol));
    std::vector<feSpace*> spaces = {u}, essentialSpaces = {};

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
    
    feBilinearForm *dudt = nullptr, *uSquared = nullptr;
    feCheckReturn(createBilinearForm(    dudt, {u}, new feSysElm_TransientMass(&scalarConstant::one)));
    feCheckReturn(createBilinearForm(uSquared, {u}, new feSysElm_MassPower(&scalarConstant::one, 2.)));
    std::vector<feBilinearForm*> forms = {dudt, uSquared};

    feLinearSystem *system;
    #if defined(HAVE_MKL)
      feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
    #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
      feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
    #else
      feCheckReturn(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
    #endif

    feNorm *errorU = nullptr;
    feCheckReturn(createNorm(errorU, LINF_ERROR, {u}, &sol, &uSol));
    std::vector<feNorm *> norms = {errorU};

    TimeIntegrator *solver;
    feCheckReturn(createTimeIntegrator(solver, scheme, NLoptions,
      system, &sol, &mesh, norms, {nullptr, 1, ""}, t0, t1, nTimeSteps));
    feCheckReturn(solver->makeSteps(nTimeSteps));

    std::vector<std::vector<double>> postProc = solver->getPostProcessingData();
    error = 0.;
    for(int j = 0; j < nTimeSteps+1; ++j) {
      feInfo("Error at step %2d (t = %1.3e) = %1.6e", j, postProc[0][j], postProc[1][j]);
      error = fmax(error, postProc[1][j]);
    }

    delete errorU;
    delete solver;
    delete system;
    for(feBilinearForm* f : forms)
      delete f;
    for(feSpace *s : spaces)
      delete s;

    return FE_STATUS_OK;
  }
};

//
// Solve dudt - lap(u) + f = 0 for manufactured solution.
//
// The "spatial" part of u is captured by Lagrange polynomials.
//
struct diffusion {

  // static double FT(const double t) { return 1. + t*t*t*t*t*t; }
  // static double DFDT(const double t) { return 6.*t*t*t*t*t; }

  // static double FT(const double t) { return exp(-t); }
  // static double DFDT(const double t) { return -exp(-t); }

  // static double FT(const double t) { return exp(t); }
  // static double DFDT(const double t) { return exp(t); }

  // static double FT(const double t) { return cos(t); }
  // static double DFDT(const double t) { return -sin(t); }

  static double FT(const double t) { return 2./(1.+t); }
  static double DFDT(const double t) { return -2./((1.+t)*(1.+t)); }

  static double uSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    // const double x = args.pos[0];
    // const double y = args.pos[1];
    const double t = args.t;
    return FT(t);
    // return FT(t) * x;
  }

  static double duSoldt_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    // const double x = args.pos[0];
    // const double y = args.pos[1];
    const double t = args.t;
    return DFDT(t);
    // return DFDT(t) * x;
  }

  static double uSrc_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    // const double x = args.pos[0];
    // const double y = args.pos[1];
    const double t = args.t;
    // const double dudt = DFDT(t) * x;
    const double dudt = DFDT(t);
    const double lapu = 0.;
    return - (dudt - lapu);
  }

  static feStatus solve(const std::string &meshFile,
                 const timeIntegratorScheme scheme,
                 const double t0,
                 const double t1,
                 const int nTimeSteps,
                 const int order,
                 const int degreeQuadrature,
                 const feNLSolverOptions &NLoptions,
                 double &error)
  {
    feFunction uSol(uSol_f), duSoldt(duSoldt_f), uSrc(uSrc_f);

    feMesh2DP1 mesh(meshFile);

    feSpace *u = nullptr, *uBord = nullptr;
    feCheckReturn(createFiniteElementSpace(uBord, &mesh, elementType::LAGRANGE, order, "U",    "Bord", degreeQuadrature, &uSol));
    feCheckReturn(createFiniteElementSpace(    u, &mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &uSol));

    std::vector<feSpace*> spaces = {u, uBord}, essentialSpaces = {uBord};

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
    
    feBilinearForm *dudt = nullptr, *diff = nullptr, *src = nullptr;
    feCheckReturn(createBilinearForm(dudt, {u}, new feSysElm_TransientMass(&scalarConstant::one)));
    feCheckReturn(createBilinearForm(diff, {u}, new feSysElm_Diffusion<2>(&scalarConstant::one)));
    feCheckReturn(createBilinearForm( src, {u}, new feSysElm_Source(&uSrc)));
    // std::vector<feBilinearForm*> forms = {dudt, diff, src};
    std::vector<feBilinearForm*> forms = {dudt, src};

    feLinearSystem *system;
    #if defined(HAVE_MKL)
      feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
    #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
      feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
    #else
      feCheckReturn(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
    #endif

    feNorm *errorU = nullptr, *errorUDot = nullptr;
    feCheckReturn(createNorm(errorU, LINF_ERROR, {u}, &sol, &uSol));
    feCheckReturn(createNorm(errorUDot, LINF_UDOT_ERROR, {u}, &sol, &duSoldt));
    std::vector<feNorm *> norms = {errorU, errorUDot};

    TimeIntegrator *solver;
    feCheckReturn(createTimeIntegrator(solver, scheme, NLoptions,
      system, &sol, &mesh, norms, {nullptr, 1, ""}, t0, t1, nTimeSteps));
    feCheckReturn(solver->makeSteps(nTimeSteps));

    std::vector<std::vector<double>> postProc = solver->getPostProcessingData();
    error = 0.;
    for(int j = 0; j < nTimeSteps+1; ++j) {
      feInfo("Error at step %2d (t = %1.3e) = %1.16e (U) = %1.16e (Udot)", j, postProc[0][j], postProc[1][j], postProc[2][j]);
      error = fmax(error, postProc[1][j]);
    }

    delete errorU;
    delete errorUDot;
    delete solver;
    delete system;
    for(feBilinearForm* f : forms)
      delete f;
    for(feSpace *s : spaces)
      delete s;

    return FE_STATUS_OK;
  }
};

namespace nonLinearDiffusion {

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
    const double cpx = cos(M_PI*x);
    const double cpy = cos(M_PI*y);
    const double spx = sin(M_PI*x);
    const double spy = sin(M_PI*y);
    const double dudt = - exp(-t) * sin(M_PI*x)*sin(M_PI*y);
    const double lapu = 2.*M_PI*M_PI* (exp(-3*t)*cpx*cpx*spx*spy*spy*spy - exp(-t)*spx*spy*(exp(-2*t)*spx*spx*spy*spy + 1) + exp(-3*t)*cpy*cpy*spx*spx*spx*spy);
    return - (dudt - lapu);
  }

  double conductivity_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double u = args.u;
    return 1. + u*u;
  }

  double ddu_conductivity_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
  {
    const double u = args.u;
    return 2.*u;
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
    feFunction uSol(uSol_f);
    feFunction uSrc(uSrc_f);
    feFunction conductivity(conductivity_f);
    feFunction ddu_conductivity(ddu_conductivity_f);

    feMesh2DP1 mesh(meshFile);
    numInteriorElements = mesh.getNumInteriorElements();

    feSpace *u = nullptr, *uBord = nullptr;
    feCheckReturn(createFiniteElementSpace(uBord, &mesh, elementType::LAGRANGE, order, "U",    "Bord", degreeQuadrature, &uSol));
    feCheckReturn(createFiniteElementSpace(    u, &mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &scalarConstant::zero));

    std::vector<feSpace*> spaces = {uBord, u};
    std::vector<feSpace*> essentialSpaces = {uBord};

    feMetaNumber numbering(&mesh, spaces, essentialSpaces);
    feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);
    
    feBilinearForm *dudt = nullptr, *diff = nullptr, *src = nullptr;
    feCheckReturn(createBilinearForm(dudt, {u}, new feSysElm_TransientMass(&scalarConstant::one)));
    feCheckReturn(createBilinearForm(diff, {u}, new feSysElm_NonlinearDiffusion<2>(&conductivity, &ddu_conductivity)));
    feCheckReturn(createBilinearForm( src, {u}, new feSysElm_Source(&uSrc)));
    std::vector<feBilinearForm*> forms = {dudt, diff, src};

    feLinearSystem *system;
    #if defined(HAVE_MKL)
      feCheckReturn(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
    #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
      feCheckReturn(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
    #else
      feCheckReturn(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
    #endif

    feNorm *errorU_L2 = nullptr;
    feCheckReturn(createNorm(errorU_L2, L2_ERROR, {u}, &sol, &uSol));
    std::vector<feNorm *> norms = {};

    TimeIntegrator *solver;
    feCheckReturn(createTimeIntegrator(solver, scheme, NLoptions,
      system, &sol, &mesh, norms, {nullptr, 1, ""}, t0, t1, nTimeSteps));
    feCheckReturn(solver->makeSteps(nTimeSteps));

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

//
// Convergence study - Refine only the time step
//
template<typename problemType>
feStatus timeStepConvergence(std::stringstream &resultBuffer,
                             const timeIntegratorScheme scheme,
                             const double t0,
                             const double t1, 
                             const int initial_nTimeSteps,
                             const int order,
                             const int N,
                             const int degreeQuadrature,
                             const std::string &message)
{
  feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 10, 4, 1e-1};
  {
    std::vector<int> timeSteps(N);
    std::vector<double> err(N, 0.);
    int nTimeSteps = initial_nTimeSteps;

    for(int i = 0; i < N; ++i, nTimeSteps *= 2)
    {
      timeSteps[i] = nTimeSteps;
      // std::string meshFile = "../../../data/square1.msh";
      std::string meshFile = "../../../data/dummy.msh";
      feCheckReturn(problemType::solve(meshFile, scheme, t0, t1, nTimeSteps, order, degreeQuadrature, NLoptions, err[i]));
    }
    // resultBuffer << "Unsteady diffusion - Error on u - Lagrange elements P" << order << std::endl;
    resultBuffer << message << std::endl;
    computeAndPrintConvergence(1, N, err, timeSteps, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
  }
  return FE_STATUS_OK;
}

TEST(ScalarFEUnsteady, LinearODE)
{
  initialize(my_argc, my_argv);
  // setVerbose(0);
  std::string testRoot = "../../../tests/withLinearSolver/scalarFEunsteady_linearODE";
  std::stringstream resultBuffer;
  int degreeQuadrature = 8;
  double t0 = 0.;
  double t1 = 1.;
  int initial_nTimeSteps = 10;

  for(size_t i = 0; i < schemes.size(); ++i)
  {
    std::string message = "Linear ODE - Error on scalar-valued solution u - Time integrator: " + schemesStr[i];
    ASSERT_TRUE(timeStepConvergence<linearODE>(resultBuffer, schemes[i], 
      t0, t1, initial_nTimeSteps, 1, 5, degreeQuadrature, message) == FE_STATUS_OK);
  }

  EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
  finalize();
}

TEST(ScalarFEUnsteady, NonLinearODE)
{
  initialize(my_argc, my_argv);
  // setVerbose(0);
  std::string testRoot = "../../../tests/withLinearSolver/scalarFEunsteady_nonlinearODE";
  std::stringstream resultBuffer;
  int degreeQuadrature = 8;
  double t0 = 0.;
  double t1 = 1.;
  int initial_nTimeSteps = 10;

  for(size_t i = 0; i < schemes.size(); ++i)
  {
    std::string message = "Nonlinear ODE - Error on scalar-valued solution u - Time integrator: " + schemesStr[i];
    ASSERT_TRUE(timeStepConvergence<nonlinearODE>(resultBuffer, schemes[i], 
      t0, t1, initial_nTimeSteps, 1, 5, degreeQuadrature, message) == FE_STATUS_OK);
  }

  EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
  finalize();
}

// TEST(ScalarFEUnsteady, Diffusion)
// {
//   initialize(my_argc, my_argv);
//   // setVerbose(0);
//   std::string testRoot = "../../../tests/withLinearSolver/scalarFEunsteady_diffusion";
//   std::stringstream resultBuffer;
//   int degreeQuadrature = 4;
//   double t0 = 0.;
//   double t1 = 1.;
//   int initial_nTimeSteps = 10;

//   // for(size_t i = 0; i < schemes.size(); ++i)
//   // for(size_t i = 0; i < 1; ++i)
//   for(size_t i = 2; i < 3; ++i)
//   // for(size_t i = 3; i < 4; ++i)
//   {
//     std::string message = "Diffusion PDE w/ strong BC - Error on sclar-valued solution u - Time integrator: " + schemesStr[i];
//     ASSERT_TRUE(timeStepConvergence<diffusion>(resultBuffer, schemes[i], t0, t1, initial_nTimeSteps, 1, 1, degreeQuadrature, message) == FE_STATUS_OK);
//   }

//   EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
//   finalize();
// }

// TEST(ScalarFEUnsteady, NonLinearDiffusion)
// {
//   initialize(my_argc, my_argv);
//   setVerbose(0);
//   std::string testRoot = "../../../tests/withLinearSolver/scalarFEunsteady_nonlinear_bdf2";
//   std::stringstream resultBuffer;
//   int degreeQuadrature = 8;
//   double t0 = 0.;
//   double t1 = 1.;
//   int initial_nTimeSteps = 10;
//   timeIntegratorScheme scheme = timeIntegratorScheme::BDF2;
//   ASSERT_TRUE(nonLinearDiffusion::meshConvergence(resultBuffer, scheme, t0, t1, initial_nTimeSteps, 1, 4, degreeQuadrature) == FE_STATUS_OK);
//   // initial_nTimeSteps = 20;
//   // ASSERT_TRUE(diffusion::meshConvergence(resultBuffer, scheme, t0, t1, initial_nTimeSteps, 2, 4, degreeQuadrature) == FE_STATUS_OK);

//   EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
//   finalize();
// }

//
// Custom main to capture argc and argv
//
int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  my_argc = argc;
  my_argv = argv;
  return RUN_ALL_TESTS();
}