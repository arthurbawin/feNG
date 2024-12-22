
#include "feAPI.h"
#include "tests.h"

#include <gtest/gtest.h>

/*

  Evaluate the convergence of the interpolation of x^6+y^6 on the unit
  square and cube with Lagrange shape functions.

*/

static int my_argc;
static char** my_argv;

namespace {

  double fSol(const double /* t */, const std::vector<double> &pos, const std::vector<double> & /* par */)
  {
    double x = pos[0];
    double y = pos[1];
    return pow(x, 6) + pow(y, 6);
  }

  int meshConvergence(std::stringstream &resultBuffer, int dim, int order, int numMeshes, int degreeQuadrature)
  {
    if(dim != 2 && dim != 3)
      dim = 3;
      
    setVerbose(2);

    feFunction funSol(fSol);

    std::vector<int> nElm(numMeshes, 0);
    std::vector<double> L2errorU(2 * numMeshes, 0.0);

    std::string meshFile;

    //
    // Interpolate the solution on each mesh
    //
    for(int i = 0; i < numMeshes; ++i)
    {
      if(dim == 2) {
        meshFile = "../../data/square" + std::to_string(i + 1) + ".msh";
      } else {
        meshFile = "../../data/cube" + std::to_string(i + 1) + ".msh";
      }

      feMesh2DP1 mesh(meshFile);
      nElm[i] = mesh.getNumInteriorElements();

      feSpace *uDomaine;
      feCheck(createFiniteElementSpace(uDomaine, &mesh, elementType::LAGRANGE, order, "U", "Domaine", degreeQuadrature, &funSol));

      std::vector<feSpace *> spaces = {uDomaine};
      std::vector<feSpace *> essentialSpaces = {uDomaine};

      feMetaNumber numbering(&mesh, spaces, essentialSpaces);
      feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);

      sol.initializeUnknowns(&mesh);
      sol.initializeEssentialBC(&mesh);

      // Compute L2 error
      feNorm norm(L2_ERROR, {uDomaine}, &sol, &funSol);
      L2errorU[2 * i] = norm.compute();

      delete uDomaine;
    }

    // Compute the convergence rate
    for(int i = 1; i < numMeshes; ++i) {

      if(dim == 2) {
        L2errorU[2 * i + 1] = -log(L2errorU[2 * i] / L2errorU[2 * (i - 1)]) / log(sqrt(nElm[i]) / sqrt(nElm[i - 1]));
      } else {
        L2errorU[2 * i + 1] = -log(L2errorU[2 * i] / L2errorU[2 * (i - 1)]) /  log(pow((double) nElm[i] / nElm[i - 1], 1./3.));
      }
    }

    resultBuffer << "Dimension " << dim << " - Lagrange elements P" << order << std::endl;
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

TEST(Interpolation, ConvergenceInterpolation)
{
  initialize(my_argc, my_argv);
  
  std::string testRoot = "../../tests/interpolation";

  std::stringstream resultBuffer;

  int degreeQuadrature = 12;

  int dim = 2;
  meshConvergence(resultBuffer, dim, 1, 4, degreeQuadrature);
  meshConvergence(resultBuffer, dim, 2, 4, degreeQuadrature);
  meshConvergence(resultBuffer, dim, 3, 4, degreeQuadrature);
  meshConvergence(resultBuffer, dim, 4, 4, degreeQuadrature);

  dim = 3;
  meshConvergence(resultBuffer, dim, 1, 3, degreeQuadrature);
  meshConvergence(resultBuffer, dim, 2, 3, degreeQuadrature);
  meshConvergence(resultBuffer, dim, 3, 3, degreeQuadrature);
  // meshConvergence(resultBuffer, dim, 4, 3, degreeQuadrature); // P4 not yet fixed (face numbering)

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