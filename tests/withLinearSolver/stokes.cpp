
#include "feAPI.h"
#include "../tests.h"

#include <gtest/gtest.h>

/*
  Tests the implementation of mixed vector FE spaces and the Stokes equations
  Manufactured solution (polynomial) and Poiseuille flow in rectangular channel
*/

static int my_argc;
static char** my_argv;

namespace stokes {

  namespace MMS {

    void uSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/, std::vector<double> &res)
    {
      const double x = args.pos[0];
      const double y = args.pos[1];
      res[0] = pow(x, 4) * pow(y, 4);
      res[1] = -4./5. * pow(x, 3) * pow(y, 5);
    }

    double pSol_f(const feFunctionArguments &args, const std::vector<double> &/*par*/)
    {
      const double x = args.pos[0];
      const double y = args.pos[1];
      return x*x*y*y;
    }

    void uSrc_f(const feFunctionArguments &args, const std::vector<double> &/*par*/, std::vector<double> &res)
    {
      const double x = args.pos[0];
      const double y = args.pos[1];
      const double minus_dpdx[2] = {-2.*x*y*y, -2.*x*x*y};
      const double lap_u[2] = {12. * (x*x*y*y*y*y + x*x*x*x*y*y),
                              -4./5. * (6. *x*y*y*y*y*y + 20. * x*x*x*y*y*y)};
      res[0] = - (minus_dpdx[0] + lap_u[0]);
      res[1] = - (minus_dpdx[1] + lap_u[1]);
    }

    int solve(const std::string &meshFile,
              const int order,
              const int degreeQuadrature,
              const feNLSolverOptions &NLoptions,
              int &numInteriorElements,
              double &L2ErrorVelocity,
              double &L2ErrorPressure,
              const bool divergenceFormulation)
    {
      feVectorFunction uSol(uSol_f);
      feVectorFunction uSrc(uSrc_f);
      feFunction pSol(pSol_f);

      feMesh2DP1 mesh(meshFile);
      numInteriorElements = mesh.getNumInteriorElements();

      feSpace *u = nullptr, *uB = nullptr, *p = nullptr, *pB = nullptr;
      feCheck(createFiniteElementSpace( u, &mesh, elementType::VECTOR_LAGRANGE,   order, "U", "Domaine", degreeQuadrature, &uSol));
      feCheck(createFiniteElementSpace(uB, &mesh, elementType::VECTOR_LAGRANGE,   order, "U", "Bord"   , degreeQuadrature, &uSol));
      feCheck(createFiniteElementSpace( p, &mesh, elementType::LAGRANGE       , order-1, "P", "Domaine", degreeQuadrature, &pSol));
      feCheck(createFiniteElementSpace(pB, &mesh, elementType::LAGRANGE       , order-1, "P", "Bord"   , degreeQuadrature, &pSol));
      std::vector<feSpace*> spaces = {u, uB, p, pB};
      std::vector<feSpace*> essentialSpaces = {uB, pB};

      feMetaNumber numbering(&mesh, spaces, essentialSpaces);
      feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);

      feBilinearForm *divSigma = nullptr, *diffU = nullptr, *gradP = nullptr, *divU = nullptr, *source = nullptr;
      feCheck(createBilinearForm( divU, {p, u}, new feSysElm_MixedDivergence(&scalarConstant::one)));
      feCheck(createBilinearForm(source,   {u}, new feSysElm_VectorSource(&uSrc)));
      std::vector<feBilinearForm*> forms = {divU, source};

      if(divergenceFormulation) {
        UNUSED(diffU, gradP);
        feCheck(createBilinearForm(divSigma, {u, p}, new feSysElm_DivergenceNewtonianStress(&scalarConstant::one, &scalarConstant::one)));
        forms.push_back(divSigma);
      } else {
        // Laplacian formulation
        UNUSED(divSigma);
        feCheck(createBilinearForm(gradP, {u, p}, new feSysElm_MixedGradient(&scalarConstant::minusOne)));
        feCheck(createBilinearForm(diffU,    {u}, new feSysElm_VectorDiffusion(&scalarConstant::minusOne, &scalarConstant::one)));
        forms.push_back(diffU);
        forms.push_back(gradP);
      }

      feLinearSystem *system;
      #if defined(HAVE_MKL)
        feCheck(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
      #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
        feCheck(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
      #else
        feCheck(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
      #endif

      std::vector<feNorm *> norms = {};

      // feExporter *exporter;
      // feCheck(createVisualizationExporter(exporter, VTK, &numbering, &sol, &mesh, spaces));
      // feExportData exportData = {exporter, 1, "vec"};

      TimeIntegrator *solver;
      feCheck(createTimeIntegrator(solver, timeIntegratorScheme::STATIONARY, NLoptions, system, &sol, &mesh, norms, {nullptr, 1, ""}));
      feCheck(solver->makeSteps(1));

      feNorm *errorU_L2 = nullptr, *errorP_L2 = nullptr;
      feCheck(createNorm(errorU_L2, VECTOR_L2_ERROR, {u}, &sol, nullptr, &uSol));
      feCheck(createNorm(errorP_L2, L2_ERROR, {p}, &sol, &pSol));
      L2ErrorVelocity = errorU_L2->compute();
      L2ErrorPressure = errorP_L2->compute();

      delete errorU_L2;
      delete errorP_L2;
      delete solver;
      delete system;
      for(feBilinearForm* f : forms)
        delete f;
      for(feSpace *s : spaces)
        delete s;

      return 0;
    }

    int meshConvergence(std::stringstream &resultBuffer, int order, int numMeshes, int degreeQuadrature)
    {
      feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 10, 4, 1e-1};

      std::vector<int> nElm(numMeshes);
      std::vector<double> errU(numMeshes, 0.);
      std::vector<double> errP(numMeshes, 0.);

      bool divergenceFormulation = true;
      for(int i = 0; i < numMeshes; ++i)
      {
        std::string meshFile = "../../../data/square" + std::to_string(i+1) + ".msh";
        stokes::MMS::solve(meshFile, order, degreeQuadrature, NLoptions, nElm[i], errU[i], errP[i], divergenceFormulation);
      }
      resultBuffer << "Taylor-Hood elements P" << order << "-P" << order-1 << std::endl; 
      resultBuffer << "Stokes MMS - Divergence formulation - Error on u" << std::endl;
      computeAndPrintConvergence(2, numMeshes, errU, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
      resultBuffer << "Stokes MMS - Divergence formulation - Error on p" << std::endl;
      computeAndPrintConvergence(2, numMeshes, errP, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);

      divergenceFormulation = false;
      for(int i = 0; i < numMeshes; ++i)
      {
        std::string meshFile = "../../../data/square" + std::to_string(i+1) + ".msh";
        stokes::MMS::solve(meshFile, order, degreeQuadrature, NLoptions, nElm[i], errU[i], errP[i], divergenceFormulation);
      }
      resultBuffer << "Taylor-Hood elements P" << order << "-P" << order-1 << std::endl; 
      resultBuffer << "Stokes MMS - Laplacian formulation - Error on u" << std::endl;
      computeAndPrintConvergence(2, numMeshes, errU, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);
      resultBuffer << "Stokes MMS - Laplacian formulation - Error on p" << std::endl;
      computeAndPrintConvergence(2, numMeshes, errP, nElm, DEFAULT_SIGNIFICANT_DIGITS, resultBuffer);

      return 0;
    }

  } // namespace MMS

  namespace poiseuille {

    void uSol_f(const feFunctionArguments &args, const std::vector<double> &par, std::vector<double> &res)
    {
      const double y = args.pos[1];
      const double H = par[0];
      const double dpdx = par[1];
      res[0] = -dpdx/2. * y * (H - y);
      res[1] = 0.;
    }

    double pSol_f(const feFunctionArguments &args, const std::vector<double> &par)
    {
      const double x = args.pos[0];
      const double L = par[0];
      const double dpdx = par[1];
      return -dpdx * (L-x);
    }

    int solve(const std::string &meshFile,
              const int order,
              const int degreeQuadrature,
              const feNLSolverOptions &NLoptions,
              int &numInteriorElements,
              double &L2ErrorVelocity,
              double &L2ErrorPressure,
              const bool divergenceFormulation)
    {
      double H = 1.; // Height of the channel
      double L = 5.; // Length of the channel
      double dpdx = -1.0;
      feVectorFunction uSol(uSol_f, {H, dpdx});
      feFunction pSol(pSol_f, {L, dpdx});

      feMesh2DP1 mesh(meshFile);
      numInteriorElements = mesh.getNumInteriorElements();

      feSpace *u = nullptr, *p = nullptr;
      feSpace *uInlet = nullptr, *uOutlet = nullptr, *uNoSlip = nullptr;
      feCheck(createFiniteElementSpace(      u, &mesh, elementType::VECTOR_LAGRANGE,   order, "U", "Domain", degreeQuadrature, &vectorConstant::zero));
      feCheck(createFiniteElementSpace( uInlet, &mesh, elementType::VECTOR_LAGRANGE,   order, "U",  "Inlet", degreeQuadrature, &uSol));
      feCheck(createFiniteElementSpace(uOutlet, &mesh, elementType::VECTOR_LAGRANGE,   order, "U", "Outlet", degreeQuadrature, &vectorConstant::zero));
      feCheck(createFiniteElementSpace(uNoSlip, &mesh, elementType::VECTOR_LAGRANGE,   order, "U", "NoSlip", degreeQuadrature, &vectorConstant::zero));
      feCheck(createFiniteElementSpace(      p, &mesh, elementType::LAGRANGE,        order-1, "P", "Domain", degreeQuadrature, &scalarConstant::zero));

      std::vector<feSpace*> spaces = {uInlet, uNoSlip, u, p};
      std::vector<feSpace*> essentialSpaces = {uNoSlip, uInlet};

      if(divergenceFormulation) {
        // Solve for u-velocity on outlet but set v-velocity to zero
        // This in turn sets p = 0 on the outlet
        spaces.push_back(uOutlet);
        uOutlet->setEssentialComponent(1, true);
      } else {
        // Free BC on the outlet :
        // effectively sets p = 0 and dvdx = 0 on the outlet
      }

      feMetaNumber numbering(&mesh, spaces, essentialSpaces);
      feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);

      feBilinearForm *divSigma = nullptr, *diffU = nullptr, *gradP = nullptr, *divU = nullptr;
      feCheck(createBilinearForm( divU, {p, u}, new feSysElm_MixedDivergence(&scalarConstant::one)));
      std::vector<feBilinearForm*> forms = {divU};

      if(divergenceFormulation) {
        UNUSED(diffU, gradP);
        feCheck(createBilinearForm(divSigma, {u, p}, new feSysElm_DivergenceNewtonianStress(&scalarConstant::one, &scalarConstant::one)));
        forms.push_back(divSigma);
      } else {
        // Laplacian formulation
        UNUSED(divSigma);
        feCheck(createBilinearForm(gradP, {u, p}, new feSysElm_MixedGradient(&scalarConstant::minusOne)));
        feCheck(createBilinearForm(diffU,    {u}, new feSysElm_VectorDiffusion(&scalarConstant::minusOne, &scalarConstant::one)));
        forms.push_back(diffU);
        forms.push_back(gradP);
      }

      feLinearSystem *system;
      #if defined(HAVE_MKL)
        feCheck(createLinearSystem(system, MKLPARDISO, forms, numbering.getNbUnknowns()));
      #elif defined(HAVE_PETSC) && defined(PETSC_HAVE_MUMPS)
        feCheck(createLinearSystem(system, PETSC_MUMPS, forms, numbering.getNbUnknowns()));
      #else
        feCheck(createLinearSystem(system, PETSC, forms, numbering.getNbUnknowns()));
      #endif

      std::vector<feNorm *> norms = {};

      // feExporter *exporter;
      // feCheck(createVisualizationExporter(exporter, VTK, &numbering, &sol, &mesh, spaces));
      // feExportData exportData = {exporter, 1, "vec"};

      TimeIntegrator *solver;
      feCheck(createTimeIntegrator(solver, timeIntegratorScheme::STATIONARY, NLoptions, system, &sol, &mesh, norms, {nullptr, 1, ""}));
      feCheck(solver->makeSteps(1));

      feNorm *errorU_L2 = nullptr, *errorP_L2 = nullptr;
      feCheck(createNorm(errorU_L2, VECTOR_L2_ERROR, {u}, &sol, nullptr, &uSol));
      feCheck(createNorm(errorP_L2, L2_ERROR, {p}, &sol, &pSol));
      L2ErrorVelocity = errorU_L2->compute();
      L2ErrorPressure = errorP_L2->compute();

      delete errorU_L2;
      delete errorP_L2;
      delete solver;
      delete system;
      for(feBilinearForm* f : forms)
        delete f;
      for(feSpace *s : spaces)
        delete s;

      return 0;
    }

    int meshConvergence(std::stringstream &resultBuffer, int order, int degreeQuadrature, double errors[4])
    {
      UNUSED(resultBuffer);
      feNLSolverOptions NLoptions{1e-10, 1e-10, 1e4, 10, 4, 1e-1};

      int nElm;
      double errU, errP;
      std::string meshFile = "../../../data/poiseuille0.msh";

      bool divergenceFormulation = true;
      stokes::poiseuille::solve(meshFile, order, degreeQuadrature, NLoptions, nElm, errU, errP, divergenceFormulation);
      std::cout << "Taylor-Hood elements P" << order << "-P" << order-1 << std::endl;
      std::cout << "Stokes Poiseuille - Divergence formulation - Error on u" << std::endl;
      std::cout
        << std::scientific
        << std::setw(16) << std::right
        << "nElm = " << std::setprecision(6) << nElm << std::setw(16) << std::right
        << "errU = " << std::setprecision(6) << errU << std::setw(16) << std::right
        << std::endl;
      std::cout << "Stokes Poiseuille - Divergence formulation - Error on p" << std::endl;
      std::cout
        << std::scientific
        << std::setw(16) << std::right
        << "nElm = " << std::setprecision(6) << nElm << std::setw(16) << std::right
        << "errP = " << std::setprecision(6) << errP << std::setw(16) << std::right
        << std::endl;

      errors[0] = errU;
      errors[1] = errP;

      divergenceFormulation = true;
      stokes::poiseuille::solve(meshFile, order, degreeQuadrature, NLoptions, nElm, errU, errP, divergenceFormulation);
      std::cout << "Taylor-Hood elements P" << order << "-P" << order-1 << std::endl;
      std::cout << "Stokes Poiseuille - Laplacian formulation - Error on u" << std::endl;
      std::cout
        << std::scientific
        << std::setw(16) << std::right
        << "nElm = " << std::setprecision(6) << nElm << std::setw(16) << std::right
        << "errU = " << std::setprecision(6) << errU << std::setw(16) << std::right
        << std::endl;
      std::cout << "Stokes Poiseuille - Laplacian formulation - Error on p" << std::endl;
      std::cout
        << std::scientific
        << std::setw(16) << std::right
        << "nElm = " << std::setprecision(6) << nElm << std::setw(16) << std::right
        << "errP = " << std::setprecision(6) << errP << std::setw(16) << std::right
        << std::endl;

      errors[2] = errU;
      errors[3] = errP;

      return 0;
    }
  } // namespace poiseuille
}


TEST(Stokes, MMS)
{
  initialize(my_argc, my_argv);
  setVerbose(0);
  
  std::string testRoot = "../../../tests/withLinearSolver/stokes_MMS";

  std::stringstream resultBuffer;

  int degreeQuadrature = 8;
  stokes::MMS::meshConvergence(resultBuffer, 2, 4, degreeQuadrature);
  stokes::MMS::meshConvergence(resultBuffer, 3, 4, degreeQuadrature);
  stokes::MMS::meshConvergence(resultBuffer, 4, 4, degreeQuadrature);

  // Compare the resultBuffer with the file testRoot.output
  // If output file does not exist, a temporary file .tmp_output
  // is created, which may be used as the new .output file for this test.
  EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);

  finalize();
}

TEST(Stokes, Poiseuille)
{
  initialize(my_argc, my_argv);
  setVerbose(0);
  
  std::stringstream resultBuffer;

  double error_tolerance = 1e-13;

  int degreeQuadrature = 8;
  double errors[4];
  stokes::poiseuille::meshConvergence(resultBuffer, 2, degreeQuadrature, errors);
  for(int i = 0; i < 4; ++i) ASSERT_TRUE(errors[i] < error_tolerance);
  stokes::poiseuille::meshConvergence(resultBuffer, 3, degreeQuadrature, errors);
  for(int i = 0; i < 4; ++i) ASSERT_TRUE(errors[i] < error_tolerance);
  stokes::poiseuille::meshConvergence(resultBuffer, 4, degreeQuadrature, errors);
  for(int i = 0; i < 4; ++i) ASSERT_TRUE(errors[i] < error_tolerance);

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