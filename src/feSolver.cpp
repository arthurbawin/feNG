#include "feSolver.h"
#include "feSolutionContainer.h"
#include "feExporter.h"

void solveQNBDF(feSolutionContainer *solDot, feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem, 
  std::vector<feBilinearForm*> &formMatrices, std::vector<feBilinearForm*> &formResiduals,
  feSolution *sol, feMesh *mesh)
{
  bool newton = true;
  bool status = linearSystem->getRecomputeStatus();
  int iter = 0, linearSystemIter;
  double normDx, normResidual, normAxb;

  while(newton){

    linearSystem->setToZero();
    solDot->computeSolTimeDerivative(sol, linearSystem);
    linearSystem->assemble(sol);
    linearSystem->solve(&normDx, &normResidual, &normAxb, &linearSystemIter);
    //FECORRECTIONSOLUTION
    linearSystem->correctSolution(sol);
    solDot->setSol(0, sol->getSolutionCopy());

    printf("iter %2d : ||A*dx-res|| = %10.10e (%4d iter.) \t ||dx|| = %10.10e \t ||res|| = %10.10e\n", ++iter, normAxb, linearSystemIter, normDx, normResidual);

    newton = !((normDx <= tol.tolDx && normResidual <= tol.tolResidual) || iter > tol.maxIter);

    linearSystem->setRecomputeStatus(true);
  }
  if(iter > tol.maxIter){
    printf("=== ! === Not converged at iter %2d : ||A*dx-res|| = %10.10e  (%4d iter.) \t ||dx|| = %10.10e \t ||res|| = %10.10e\n", iter, normAxb, linearSystemIter, normDx, normResidual);
  } else{
    printf("Converged in %2d iterations : ||dx|| = %10.10e \t ||res|| = %10.10e\n", iter, normDx, normResidual);
  }
  linearSystem->setRecomputeStatus(status);
}

void solveStationary(double *normL2, feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem, 
  std::vector<feBilinearForm*> &formMatrices, std::vector<feBilinearForm*> &formResiduals,
  feSolution *sol, std::vector<feNorm*> &norms, feMesh *mesh)
{
	linearSystem->setRecomputeStatus(true);
  sol->initializeUnknowns(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber);
  // Initialize container
  int nSol = 1;
  feStationarySolution *solSta = new feStationarySolution(nSol, sol->getCurrentTime(), metaNumber);
  solSta->initialize(sol, mesh, metaNumber);
  // Solve
  printf("\nSteady state solution - recomputeMatrix = %s\n", linearSystem->getRecomputeStatus() ? "true" : "false");
  solveQNBDF(solSta, tol, metaNumber, linearSystem, formMatrices, formResiduals, sol, mesh);
  // Compute L2 norm of solution(s)
  for(auto *norm : norms){
    norm->computeL2Norm(metaNumber, sol, mesh);
    *normL2 = norm->getNorm(); // TODO : fix this
  }
}

/* Nom cryptique : 
  - reset sol->soldot
  - recalcule le residu
  - assigne le residu a solContainer->F */
void fePstClc(feSolution *sol, feLinearSystem *linearSystem, feSolutionContainer *solContainer){
  sol->setSolDotToZero();
  linearSystem->setResidualToZero();
  linearSystem->assembleResiduals(sol);
  linearSystem->assignResidualToDCResidual(solContainer);
}

double f = 0.5;
bool K1K2 = false;

void solveBDF2(std::vector<double> &normL2, feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem, 
  std::vector<feBilinearForm*> &formMatrices, std::vector<feBilinearForm*> &formResiduals,
  feSolution *sol, std::vector<feNorm*> &norms, feMesh *mesh, std::vector<feSpace*> &spaces)
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
  // norm->computeL2Norm(metaNumber, sol, mesh);
  // normL2[0] = norm->getNorm();
  for(auto *norm : norms){
    norm->computeL2Norm(metaNumber, sol, mesh);
    normL2[0] = norm->getNorm(); // TODO : fix this
  }

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
    fePstClc(sol, linearSystem, solBDF2);
    // Compute L2 norm of the solution
    sol->setSolFromContainer(solBDF2);
    // norm->computeL2Norm(metaNumber, sol, mesh);
    // normL2[iTime] = norm->getNorm();
    // int cnt = 0;
    for(auto *norm : norms){
      norm->computeL2Norm(metaNumber, sol, mesh);
      normL2[iTime] = norm->getNorm(); // TODO : fix this
      // std::cout<<"Norme "<<cnt++<<" = "<<normL2[iTime]<<std::endl;
    }

    // std::string vtkFile = "../../data/taylorGreenUnsteady" + std::to_string(iTime) + ".vtk";
    // feExporterVTK writer(vtkFile, mesh, sol, metaNumber, spaces);
  }

  delete solBDF2;
}

void solveDC3(std::vector<double> &normL2BDF2, std::vector<double> &normL2DC3,feTolerances tol, feMetaNumber *metaNumber, 
  feLinearSystem *linearSystem, std::vector<feBilinearForm*> &formMatrices, std::vector<feBilinearForm*> &formResiduals,
  feSolution *sol, std::vector<feNorm*> &norms, feMesh *mesh, std::vector<feSpace*> &spaces)
{
  linearSystem->setRecomputeStatus(true);
  double dt = sol->getTimeStep();
  double nSteps = sol->getNbTimeSteps();
  double t0 = sol->getInitialTime();

  sol->initializeUnknowns(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber);

  int nSol = 5;
  feSolutionBDF2 *solBDF2 = new feSolutionBDF2(nSol, sol->getCurrentTime(), metaNumber);
  feSolutionDCF *solDC3 = new feSolutionDCF(nSol, sol->getCurrentTime(), metaNumber);
  solBDF2->initialize(sol, mesh, metaNumber);
  solDC3->initialize(sol, mesh, metaNumber);

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
  solDC3->rotate(dt);
  solBDF2->initialize(sol, mesh, metaNumber);
  solDC3->initialize(sol, mesh, metaNumber);

  sol->setSolFromContainer(solBDF2);
  // norm->computeL2Norm(metaNumber, sol, mesh);
  // normL2BDF2[0] = norm->getNorm();

  sol->setSolFromContainer(solDC3);
  // norm->computeL2Norm(metaNumber, sol, mesh);
  // normL2DC3[0] = norm->getNorm();

  // std::string vtkFile = "../../data/cylFine" + std::to_string(0) + ".vtk";
  // feExporterVTK writer(vtkFile, mesh, sol, metaNumber, spaces);
  // int cnt = 0;
  for(auto *norm : norms){
    norm->computeL2Norm(metaNumber, sol, mesh);
    normL2BDF2[0] = norm->getNorm();
    normL2DC3[0] = norm->getNorm();
    // normL2DC3[0+cnt] = norm->getNorm();
    // normL2[iTime] = norm->getNorm(); // TODO : fix this
    // std::cout<<"Norme "<<cnt++<<" = "<<normL2DC3[0+cnt]<<std::endl;
  }

  // Start integration with a BDF2
  for(int iTime = 1; iTime <= 2; ++iTime){
    if(K1K2)
      dt = tK1K2[iTime+1] - tK1K2[iTime];
    solBDF2->rotate(dt);
    solDC3->rotate(dt);
    solDC3->initialize(sol, mesh, metaNumber);
    if(iTime==1)
      linearSystem->setRecomputeStatus(true);
    else
      linearSystem->setRecomputeStatus(false);
    // FESOLVEBDF2NL
    initializeBDF2(sol, metaNumber, mesh, solBDF2);
    printf("\n");
    printf("Étape 1 - recomputeMatrix = %s : Solution BDF2 - t = %6.6e\n", linearSystem->getRecomputeStatus() ? "true" : "false", sol->getCurrentTime());
    solveQNBDF(solBDF2, tol, metaNumber, linearSystem, formMatrices, formResiduals, sol, mesh);
    fePstClc(sol, linearSystem, solBDF2);
    // Compute L2 norm of the solution
    sol->setSolFromContainer(solBDF2);
    // norm->computeL2Norm(metaNumber, sol, mesh);
    // normL2BDF2[iTime] = norm->getNorm();
    sol->setSolFromContainer(solDC3);
    // norm->computeL2Norm(metaNumber, sol, mesh);
    // normL2DC3[iTime] = norm->getNorm();

    // std::string vtkFile = "../../data/cylFine" + std::to_string(iTime) + ".vtk";
    // feExporterVTK writer(vtkFile, mesh, sol, metaNumber, spaces);
    // int cnt = 0;
    for(auto *norm : norms){
      norm->computeL2Norm(metaNumber, sol, mesh);
      normL2BDF2[iTime] = norm->getNorm();
      normL2DC3[iTime] = norm->getNorm();
      // normL2DC3[3*iTime+cnt] = norm->getNorm();
      // normL2[iTime] = norm->getNorm(); // TODO : fix this
      // std::cout<<"Norme "<<cnt++<<" = "<<normL2DC3[3*iTime+cnt]<<std::endl;
    }
  }
  // Continue with DC3
  for(int iTime = 3; iTime < nSteps; ++iTime){
  // for(int iTime = 3; iTime < 4; ++iTime){
    if(K1K2)
      dt = tK1K2[iTime+1] - tK1K2[iTime];
    solBDF2->rotate(dt);
    solDC3->rotate(dt);
    // if(iTime==3)
    //   linearSystem->setRecomputeStatus(true);
    // else
    //   linearSystem->setRecomputeStatus(false);
    // FESOLVEDC3F
    initializeBDF2(sol, metaNumber, mesh, solBDF2);
    printf("\n");
    printf("Étape 2 - recomputeMatrix = %s : Solution DC3F - t = %6.6e\n", linearSystem->getRecomputeStatus() ? "true" : "false", sol->getCurrentTime());
    solveQNBDF(solBDF2, tol, metaNumber, linearSystem, formMatrices, formResiduals, sol, mesh);
    fePstClc(sol, linearSystem, solBDF2);
    linearSystem->setRecomputeStatus(false);
    initializeDC3F(sol, metaNumber, mesh, solBDF2, solDC3);
    solveQNBDF(solDC3, tol, metaNumber, linearSystem, formMatrices, formResiduals, sol, mesh);
    fePstClc(sol, linearSystem, solDC3);
    // // Compute L2 norm of the solution
    sol->setSolFromContainer(solBDF2);
    // norm->computeL2Norm(metaNumber, sol, mesh);
    // normL2BDF2[iTime] = norm->getNorm();
    sol->setSolFromContainer(solDC3);
    // norm->computeL2Norm(metaNumber, sol, mesh);
    // normL2DC3[iTime] = norm->getNorm();

    // std::string vtkFile = "../../data/cylFine" + std::to_string(iTime) + ".vtk";
    // feExporterVTK writer(vtkFile, mesh, sol, metaNumber, spaces);
    // int cnt = 0;
    for(auto *norm : norms){
      norm->computeL2Norm(metaNumber, sol, mesh);
      normL2BDF2[iTime] = norm->getNorm();
      normL2DC3[iTime] = norm->getNorm();
      // normL2DC3[3*iTime+cnt] = norm->getNorm();
      // normL2[iTime] = norm->getNorm(); // TODO : fix this
      // std::cout<<"Norme "<<cnt++<<" = "<<normL2DC3[3*iTime+cnt]<<std::endl;
    }
  }

  delete solDC3;
  delete solBDF2;
}