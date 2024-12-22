
#include "feAPI.h"
#include "../tests.h"

#include <gtest/gtest.h>

/*

  Number the degrees of freedom for a scalar field on the unit square
  with Lagrange finite elements of degree 4.
  Print the code (DOF_ESSENTIAL = -1, DOF_UNKNOWN = -2, DOF_NOT_ASSIGNED = -3)
  of each dof and the final numbering for this field.

*/

namespace {

  int createDOFNumbering(std::stringstream &resultBuffer)
  {
    setVerbose(3);

    const char *meshFile = "../../../data/square1.msh";
    feMesh2DP1 mesh(meshFile);

    // Define a finite element space on each subdomain of the computational domain.
    int order = 4;
    int degreeQuadrature = 1;
    feConstantFunction funZero(0.);
    feSpace *uBord, *uDomaine;
    feCheck(createFiniteElementSpace(uBord, &mesh, elementType::LAGRANGE, order, "U", "Bord",
                                     degreeQuadrature, &funZero));
    feCheck(createFiniteElementSpace(uDomaine, &mesh, elementType::LAGRANGE, order, "U", "Domaine",
                                     degreeQuadrature, &funZero));

    std::vector<feSpace *> spaces = {uBord, uDomaine};
    std::vector<feSpace *> essentialSpaces = {uBord};

    // Create the DOFs numbering for each field. Here there is only one field U.
    feMetaNumber numbering(&mesh, spaces, essentialSpaces);

    // Write the output
    numbering.printCodes(resultBuffer);
    numbering.printNumberings(resultBuffer);

    // Free the used memory
    delete uDomaine;
    delete uBord;

    return 0;
  }

}

TEST(DOF_Numbering, LagrangeP4Square) {

  std::string testRoot = "../../../tests/withoutLinearSolver/dofNumbering";

  std::stringstream resultBuffer;

  createDOFNumbering(resultBuffer);

  // Compare the resultBuffer with the file testRoot.output
  // If output file does not exist, a temporary file .tmp_output
  // is created, which may be used as the new .output file for this test.
  EXPECT_EQ(compareOutputFiles(testRoot, resultBuffer), 0);
}

