#include "feSolver.h"
#include "feSolutionContainer.h"

void solveQNBDF(feSolutionContainer *solDot, feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem, 
  std::vector<feBilinearForm*> &formMatrices, std::vector<feBilinearForm*> &formResiduals,
  feSolution *sol, feMesh *mesh)
{
  bool newton = true;
  bool status = linearSystem->getRecomputeStatus();;
  int iter = 0;
  double normDx, normResidual;

  while(newton){

    linearSystem->setToZero();
    solDot->computeSolTimeDerivative(sol);
    linearSystem->assemble(sol);
    linearSystem->solve(&normDx, &normResidual);
    linearSystem->correctSolution(sol);
    solDot->setSol(0, sol->getSolutionCopy());

    printf("iter %2d : ||dx|| = %10.10e \t ||res|| = %10.10e\n", ++iter, normDx, normResidual);

    newton = !((normDx <= tol.tolDx && normResidual <= tol.tolResidual) || iter > tol.maxIter);

    linearSystem->setRecomputeStatus(false);
  }
  linearSystem->setRecomputeStatus(status);
}

double solveStationary(feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem, 
  std::vector<feBilinearForm*> &formMatrices, std::vector<feBilinearForm*> &formResiduals,
  feSolution *sol, feNorm *norm, feMesh *mesh)
{
	linearSystem->setRecomputeStatus(true);

  sol->initializeUnknowns(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber);

  int nSol = 2;
  feStationarySolution *solSta = new feStationarySolution(nSol, sol->getCurrentTime(), metaNumber);
  solSta->initialize(sol, mesh, metaNumber);

  // FESolveSTATNL
  printf("\n");
  printf("Étape 1 - recomputeMatrix = %s : Solution stationnaire\n", linearSystem->getRecomputeStatus() ? "true" : "false");
  solveQNBDF(solSta, tol, metaNumber, linearSystem, formMatrices, formResiduals, sol, mesh);

  norm->computeL2Norm(metaNumber, sol, mesh);

  return norm->getNorm();
}