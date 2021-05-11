#include "feSolver.h"
#include "feSolutionContainer.h"

void solveQNBDF(feSolutionContainer *solDot, feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem, 
  std::vector<feBilinearForm*> &formMatrices, std::vector<feBilinearForm*> &formResiduals,
  feSolution *sol, feMesh *mesh)
{
  bool newton = true;
  bool status = linearSystem->getRecomputeStatus();
  int iter = 0;
  double normDx, normResidual;

  while(newton){

    linearSystem->setToZero();
    solDot->computeSolTimeDerivative(sol);
    linearSystem->assemble(sol);
    linearSystem->solve(&normDx, &normResidual);
    //FECORRECTIONSOLUTION
    linearSystem->correctSolution(sol);
    solDot->setSol(0, sol->getSolutionCopy());

    printf("iter %2d : ||dx|| = %10.10e \t ||res|| = %10.10e\n", ++iter, normDx, normResidual);

    newton = !((normDx <= tol.tolDx && normResidual <= tol.tolResidual) || iter > tol.maxIter);

    linearSystem->setRecomputeStatus(false);
  }
  if(iter > tol.maxIter){
    printf("=== ! === Not converged at iter %2d : ||dx|| = %10.10e \t ||res|| = %10.10e\n", iter, normDx, normResidual);
  } else{
    printf("Converged in %2d iterations : ||dx|| = %10.10e \t ||res|| = %10.10e\n", iter, normDx, normResidual);
  }
  linearSystem->setRecomputeStatus(status);
}

void solveStationary(double *normL2, feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem, 
  std::vector<feBilinearForm*> &formMatrices, std::vector<feBilinearForm*> &formResiduals,
  feSolution *sol, feNorm *norm, feMesh *mesh)
{
	linearSystem->setRecomputeStatus(true);
  sol->initializeUnknowns(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber);
  // Initialize container
  int nSol = 2;
  feStationarySolution *solSta = new feStationarySolution(nSol, sol->getCurrentTime(), metaNumber);
  solSta->initialize(sol, mesh, metaNumber);
  // Solve
  printf("\nÉtape 1 - recomputeMatrix = %s : Solution stationnaire\n", linearSystem->getRecomputeStatus() ? "true" : "false");
  solveQNBDF(solSta, tol, metaNumber, linearSystem, formMatrices, formResiduals, sol, mesh);
  // Compute L2 norm of solution
  norm->computeL2Norm(metaNumber, sol, mesh);
  *normL2 = norm->getNorm();
}

double f = 0.0;
bool K1K2 = false;

void solveBDF2NL(std::vector<double> &normL2, feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem, 
  std::vector<feBilinearForm*> &formMatrices, std::vector<feBilinearForm*> &formResiduals,
  feSolution *sol, feNorm *norm, feMesh *mesh)
{
  linearSystem->setRecomputeStatus(true);
  double dt = sol->getTimeStep();
  double nSteps = sol->getNbTimeSteps();
  double t0 = sol->getInitialTime();

  sol->initializeUnknowns(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber);

  int nSol = 5;
  feSolutionBDF2 *solBDF2 = new feSolutionBDF2(nSol, sol->getCurrentTime(), metaNumber);
  solBDF2->initialize(sol, mesh, metaNumber);

  std::vector<double> tK1K2(nSteps+1, 0.0);
  for(int i = 0; i < nSteps+1; ++i)
    tK1K2[i] = t0 + i*dt;
  for(int i = 0; i < nSteps; i+=2)
    if(i+2 < nSteps+1)
      tK1K2[i+1] = tK1K2[i] + f*(tK1K2[i+2] - tK1K2[i]);

  // Initialization
  if(K1K2)
    dt = tK1K2[1] - tK1K2[0];
  solBDF2->rotate(dt);
  solBDF2->initialize(sol, mesh, metaNumber);
  sol->setSolFromContainer(solBDF2);
  norm->computeL2Norm(metaNumber, sol, mesh);
  normL2[0] = norm->getNorm();

  // bool status = true;
  for(int iTime = 1; iTime < nSteps; ++iTime){
    if(K1K2)
      dt = tK1K2[iTime+1] - tK1K2[iTime];
    solBDF2->rotate(dt);
    if(iTime==1)
      linearSystem->setRecomputeStatus(true);
    else
      linearSystem->setRecomputeStatus(false);
    // FESOLVEBDF2NL
    initializeBDF2(sol, metaNumber, mesh, solBDF2);
    printf("\n");
    printf("Étape 1 - recomputeMatrix = %s : Solution BDF2 - t = %6.6e\n", linearSystem->getRecomputeStatus() ? "true" : "false", sol->getCurrentTime());
    solveQNBDF(solBDF2, tol, metaNumber, linearSystem, formMatrices, formResiduals, sol, mesh);
    // FEPSTCLC
    sol->setSolDotToZero();
    linearSystem->setResidualToZero();
    linearSystem->assembleResiduals(sol);
    linearSystem->assignResidualToDCResidual(solBDF2);
    // Compute L2 norm of the solution
    sol->setSolFromContainer(solBDF2);
    norm->computeL2Norm(metaNumber, sol, mesh);
    normL2[iTime] = norm->getNorm();
  }
}