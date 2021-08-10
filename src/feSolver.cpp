#include "feSolver.h"
#include "feSolutionContainer.h"
#include "feExporter.h"

void solveQNBDF(feSolutionContainer *solDot, feTolerances tol, feMetaNumber *metaNumber,
                feLinearSystem *linearSystem, feSolution *sol, feMesh *mesh) {
  bool newton = true;
  bool status = linearSystem->getRecomputeStatus();
  int iter = 0, linearSystemIter;
  double normDx, normResidual, normAxb;

  while(newton) {
    linearSystem->setToZero();
    solDot->computeSolTimeDerivative(sol, linearSystem);
    linearSystem->assemble(sol);
    linearSystem->solve(&normDx, &normResidual, &normAxb, &linearSystemIter);
    // FECORRECTIONSOLUTION
    linearSystem->correctSolution(sol);
    solDot->setSol(0, sol->getSolutionCopy());

    printf(
      "iter %2d : ||A*dx-res|| = %10.10e (%4d iter.) \t ||dx|| = %10.10e \t ||res|| = %10.10e\n",
      ++iter, normAxb, linearSystemIter, normDx, normResidual);

    newton = !((normDx <= tol.tolDx && normResidual <= tol.tolResidual) || iter > tol.maxIter);

    linearSystem->setRecomputeStatus(true);
  }
  if(iter > tol.maxIter) {
    printf("=== ! === Not converged at iter %2d : ||A*dx-res|| = %10.10e  (%4d iter.) \t ||dx|| = "
           "%10.10e \t ||res|| = %10.10e\n",
           iter, normAxb, linearSystemIter, normDx, normResidual);
  } else {
    printf("Converged in %2d iterations : ||dx|| = %10.10e \t ||res|| = %10.10e\n", iter, normDx,
           normResidual);
  }
  linearSystem->setRecomputeStatus(status);
}

void solveStationary(double *normL2, feTolerances tol, feMetaNumber *metaNumber,
                     feLinearSystem *linearSystem, std::vector<feBilinearForm *> &formMatrices,
                     std::vector<feBilinearForm *> &formResiduals, feSolution *sol,
                     std::vector<feNorm *> &norms, feMesh *mesh) {
  linearSystem->setRecomputeStatus(true);
  sol->initializeUnknowns(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber);
  // Initialize container
  int nSol = 1;
  feStationarySolution *solSta = new feStationarySolution(nSol, sol->getCurrentTime(), metaNumber);
  solSta->initialize(sol, mesh, metaNumber);
  // Solve
  printf("\nSteady state solution - recomputeMatrix = %s\n",
         linearSystem->getRecomputeStatus() ? "true" : "false");
  solveQNBDF(solSta, tol, metaNumber, linearSystem, sol, mesh);
  // Compute L2 norm of solution(s)
  for(auto *norm : norms) {
    norm->computeL2Norm(metaNumber, sol, mesh);
    *normL2 = norm->getNorm(); // TODO : fix this
  }
}

StationarySolver::StationarySolver(feTolerances tol, feMetaNumber *metaNumber,
                                   feLinearSystem *linearSystem, feSolution *sol,
                                   std::vector<feNorm *> &norms, feMesh *mesh)
  : TimeIntegrator(tol, metaNumber, linearSystem, sol, norms, mesh) {
  // Initialize the solution container with a single solution
  int nSol = 1;
  _solutionContainer = new feStationarySolution(nSol, _sol->getCurrentTime(), _metaNumber);
  _solutionContainer->initialize(_sol, _mesh, _metaNumber);

  _normL2.resize(norms.size());
  for(auto &n : _normL2) n.resize(1, 0.);
}

void StationarySolver::makeSteps(int nSteps, std::vector<feSpace *> &spaces) {
  printf("StationarySolver : Advancing %d steps\n", nSteps);
  _linearSystem->setRecomputeStatus(true);
  _sol->initializeUnknowns(_mesh, _metaNumber);
  _sol->initializeEssentialBC(_mesh, _metaNumber);
  // Solve
  printf("\nSteady state solution - recomputeMatrix = %s\n",
         _linearSystem->getRecomputeStatus() ? "true" : "false");
  solveQNBDF(_solutionContainer, _tol, _metaNumber, _linearSystem, _sol, _mesh);
  // Compute L2 norm of solution(s)
  for(size_t i = 0; i < _norms.size(); ++i) {
    _norms[i]->computeL2Norm(_metaNumber, _sol, _mesh);
    _normL2[i][0] = _norms[i]->getNorm();
  }
}

/* Nom cryptique :
  - reset sol->soldot
  - recalcule le residu
  - assigne le residu a solContainer->F */
void fePstClc(feSolution *sol, feLinearSystem *linearSystem, feSolutionContainer *solContainer) {
  sol->setSolDotToZero();
  linearSystem->setResidualToZero();
  linearSystem->assembleResiduals(sol);
  linearSystem->assignResidualToDCResidual(solContainer);
}

double f = 0.5;
bool K1K2 = false;

void solveBDF2(std::vector<double> &normL2, feTolerances tol, feMetaNumber *metaNumber,
               feLinearSystem *linearSystem, std::vector<feBilinearForm *> &formMatrices,
               std::vector<feBilinearForm *> &formResiduals, feSolution *sol,
               std::vector<feNorm *> &norms, feMesh *mesh, std::vector<feSpace *> &spaces) {
  linearSystem->setRecomputeStatus(true);
  double dt = sol->getTimeStep();
  double nSteps = sol->getNbTimeSteps();
  double t0 = sol->getInitialTime();

  sol->initializeUnknowns(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber);

  int nSol = 5;
  feSolutionBDF2 *solBDF2 = new feSolutionBDF2(nSol, sol->getCurrentTime(), metaNumber);
  solBDF2->initialize(sol, mesh, metaNumber);

  std::vector<double> tK1K2(nSteps + 1, 0.0);
  for(int i = 0; i < nSteps + 1; ++i) tK1K2[i] = t0 + i * dt;
  for(int i = 0; i < nSteps; i += 2)
    if(i + 2 < nSteps + 1) tK1K2[i + 1] = tK1K2[i] + f * (tK1K2[i + 2] - tK1K2[i]);

  // Initialization
  if(K1K2) dt = tK1K2[1] - tK1K2[0];
  solBDF2->rotate(dt);
  solBDF2->initialize(sol, mesh, metaNumber);
  sol->setSolFromContainer(solBDF2);
  // norm->computeL2Norm(metaNumber, sol, mesh);
  // normL2[0] = norm->getNorm();
  for(auto *norm : norms) {
    norm->computeL2Norm(metaNumber, sol, mesh);
    normL2[0] = norm->getNorm(); // TODO : fix this
  }

  // bool status = true;
  for(int iTime = 1; iTime < nSteps; ++iTime) {
    if(K1K2) dt = tK1K2[iTime + 1] - tK1K2[iTime];
    solBDF2->rotate(dt);
    if(iTime == 1)
      linearSystem->setRecomputeStatus(true);
    else
      linearSystem->setRecomputeStatus(false);
    // FESOLVEBDF2NL
    initializeBDF2(sol, metaNumber, mesh, solBDF2);
    printf("\n");
    printf("Étape 1 - recomputeMatrix = %s : Solution BDF2 - t = %6.6e\n",
           linearSystem->getRecomputeStatus() ? "true" : "false", sol->getCurrentTime());
    solveQNBDF(solBDF2, tol, metaNumber, linearSystem, sol, mesh);
    fePstClc(sol, linearSystem, solBDF2);
    // Compute L2 norm of the solution
    sol->setSolFromContainer(solBDF2);
    // norm->computeL2Norm(metaNumber, sol, mesh);
    // normL2[iTime] = norm->getNorm();
    // int cnt = 0;
    for(auto *norm : norms) {
      norm->computeL2Norm(metaNumber, sol, mesh);
      normL2[iTime] = norm->getNorm(); // TODO : fix this
      // std::cout<<"Norme "<<cnt++<<" = "<<normL2[iTime]<<std::endl;
    }

    // std::string vtkFile = "../../data/cylindre" + std::to_string(iTime) + ".vtk";
    // feExporterVTK writer(vtkFile, mesh, sol, metaNumber, spaces);
  }

  delete solBDF2;
}

BDF2Solver::BDF2Solver(feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem,
                       feSolution *sol, std::vector<feNorm *> &norms, feMesh *mesh, double t0,
                       double tEnd, int nTimeSteps)
  : TimeIntegrator(tol, metaNumber, linearSystem, sol, norms, mesh, t0, tEnd, nTimeSteps) {
  // Initialize the solution container
  int nSol = 3;
  _solutionContainer = new feSolutionBDF2(nSol, _sol->getCurrentTime(), _metaNumber);
  _solutionContainer->initialize(_sol, _mesh, _metaNumber);

  _normL2.resize(norms.size());
  for(auto &n : _normL2) n.resize(_nTimeSteps, 0.);

  printf("Initializing BDF2 solver : integrating from t0 = %f to tEnd = %f in %d steps\n", _t0,
         _tEnd, _nTimeSteps);
}

void BDF2Solver::makeSteps(int nSteps, std::vector<feSpace *> &spaces) {
  printf("BDF2 : Advancing %d steps from t = %f to t = %f\n", nSteps, _tCurrent,
         _tCurrent + nSteps * _dt);

  if(_currentStep == 0) {
    // Initialization and first step
    _linearSystem->setRecomputeStatus(true);
    _solutionContainer->initialize(_sol, _mesh, _metaNumber);
    _solutionContainer->rotate(_dt);
    _sol->setSolFromContainer(_solutionContainer);
    _tCurrent += _dt;
    ++_currentStep;

    for(size_t i = 0; i < _norms.size(); ++i) {
      _norms[i]->computeL2Norm(_metaNumber, _sol, _mesh);
      _normL2[i][0] = _norms[i]->getNorm();
    }

    printf("Current step = %d/%d : t = %f\n", _currentStep, nSteps, _tCurrent);

    --nSteps; // To advance the same number of steps than if currentStep != 0
  }

  for(int i = 0; i < nSteps; ++i) {
    _solutionContainer->rotate(_dt);
    if(i == 0)
      _linearSystem->setRecomputeStatus(true);
    else
      _linearSystem->setRecomputeStatus(false);
    initializeBDF2(_sol, _metaNumber, _mesh, dynamic_cast<feSolutionBDF2 *>(_solutionContainer));
    printf("\n");
    printf("Étape 1 - recomputeMatrix = %s : Solution BDF2 - t = %6.6e\n",
           _linearSystem->getRecomputeStatus() ? "true" : "false", _sol->getCurrentTime());
    solveQNBDF(_solutionContainer, _tol, _metaNumber, _linearSystem, _sol, _mesh);
    fePstClc(_sol, _linearSystem, _solutionContainer);

    // Compute L2 norm of the solution
    _sol->setSolFromContainer(_solutionContainer);
    for(size_t i = 0; i < _norms.size(); ++i) {
      _norms[i]->computeL2Norm(_metaNumber, _sol, _mesh);
      _normL2[i][_currentStep] = _norms[i]->getNorm();
    }

    _tCurrent += _dt;
    ++_currentStep;
    printf("Current step = %d : t = %f\n", _currentStep, _tCurrent);

    std::string vtkFile = "../../data/cylindreJF/cylindreAdapt" + std::to_string(_currentStep) + ".vtk";
    feExporterVTK writer(vtkFile, _mesh, _sol, _metaNumber, spaces);
  }
}

BDF1Solver::BDF1Solver(feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem,
                       feSolution *sol, std::vector<feNorm *> &norms, feMesh *mesh, double t0,
                       double tEnd, int nTimeSteps)
  : TimeIntegrator(tol, metaNumber, linearSystem, sol, norms, mesh, t0, tEnd, nTimeSteps) {
  // Initialize the solution container
  int nSol = 3;
  _solutionContainer = new feSolutionBDF1(nSol, _sol->getCurrentTime(), _metaNumber);
  _solutionContainer->initialize(_sol, _mesh, _metaNumber);

  _normL2.resize(norms.size());
  for(auto &n : _normL2) n.resize(_nTimeSteps, 0.);

  printf("Initializing BDF1 solver : integrating from t0 = %f to tEnd = %f in %d steps\n", _t0,
         _tEnd, _nTimeSteps);
}

void BDF1Solver::makeSteps(int nSteps, std::vector<feSpace *> &spaces) {
  printf("BDF1 : Advancing %d steps from t = %f to t = %f\n", nSteps, _tCurrent,
         _tCurrent + nSteps * _dt);

  if(_currentStep == 0) {
    // Initialization and first step
    _linearSystem->setRecomputeStatus(true);
    _solutionContainer->initialize(_sol, _mesh, _metaNumber);
    _solutionContainer->rotate(_dt);
    _sol->setSolFromContainer(_solutionContainer);
    _tCurrent += _dt;
    ++_currentStep;

    for(size_t i = 0; i < _norms.size(); ++i) {
      _norms[i]->computeL2Norm(_metaNumber, _sol, _mesh);
      _normL2[i][0] = _norms[i]->getNorm();
    }

    printf("Current step = %d/%d : t = %f\n", _currentStep, nSteps, _tCurrent);

    --nSteps; // To advance the same number of steps than if currentStep != 0
  }

  for(int i = 0; i < nSteps; ++i) {
    _solutionContainer->rotate(_dt);
    if(i == 0)
      _linearSystem->setRecomputeStatus(true);
    else
      _linearSystem->setRecomputeStatus(false);
    initializeBDF1(_sol, _metaNumber, _mesh, dynamic_cast<feSolutionBDF1 *>(_solutionContainer));
    printf("\n");
    printf("Étape 1 - recomputeMatrix = %s : Solution BDF1 - t = %6.6e\n",
           _linearSystem->getRecomputeStatus() ? "true" : "false", _sol->getCurrentTime());
    solveQNBDF(_solutionContainer, _tol, _metaNumber, _linearSystem, _sol, _mesh);
    fePstClc(_sol, _linearSystem, _solutionContainer);

    // Compute L2 norm of the solution
    _sol->setSolFromContainer(_solutionContainer);
    for(size_t i = 0; i < _norms.size(); ++i) {
      _norms[i]->computeL2Norm(_metaNumber, _sol, _mesh);
      _normL2[i][_currentStep] = _norms[i]->getNorm();
    }

    _tCurrent += _dt;
    ++_currentStep;
    printf("Current step = %d : t = %f\n", _currentStep, _tCurrent);

    std::string vtkFile = "../../data/cylindreJF/cylindreAdapt" + std::to_string(_currentStep) + ".vtk";
    feExporterVTK writer(vtkFile, _mesh, _sol, _metaNumber, spaces);
  }
}

void solveBDF1(std::vector<double> &normL2, feTolerances tol, feMetaNumber *metaNumber,
               feLinearSystem *linearSystem, std::vector<feBilinearForm *> &formMatrices,
               std::vector<feBilinearForm *> &formResiduals, feSolution *sol,
               std::vector<feNorm *> &norms, feMesh *mesh, std::vector<feSpace *> &spaces) {
  linearSystem->setRecomputeStatus(true);
  double dt = sol->getTimeStep();
  double nSteps = sol->getNbTimeSteps();
  // double t0 = sol->getInitialTime();

  sol->initializeUnknowns(mesh, metaNumber);
  sol->initializeEssentialBC(mesh, metaNumber);

  int nSol = 5;
  feSolutionBDF1 *solBDF1 = new feSolutionBDF1(nSol, sol->getCurrentTime(), metaNumber);
  solBDF1->initialize(sol, mesh, metaNumber);

  // std::vector<double> tK1K2(nSteps+1, 0.0);
  // for(int i = 0; i < nSteps+1; ++i)
  //   tK1K2[i] = t0 + i*dt;
  // for(int i = 0; i < nSteps; i+=2)
  //   if(i+2 < nSteps+1)
  //     tK1K2[i+1] = tK1K2[i] + f*(tK1K2[i+2] - tK1K2[i]);

  // Initialization
  // if(K1K2)
  //   dt = tK1K2[1] - tK1K2[0];
  solBDF1->rotate(dt);
  solBDF1->initialize(sol, mesh, metaNumber);
  sol->setSolFromContainer(solBDF1);
  // norm->computeL2Norm(metaNumber, sol, mesh);
  // normL2[0] = norm->getNorm();
  for(auto *norm : norms) {
    norm->computeL2Norm(metaNumber, sol, mesh);
    normL2[0] = norm->getNorm(); // TODO : fix this
  }

  // bool status = true;
  for(int iTime = 1; iTime < nSteps; ++iTime) {
    // if(K1K2)
    //   dt = tK1K2[iTime+1] - tK1K2[iTime];
    solBDF1->rotate(dt);
    if(iTime == 1)
      linearSystem->setRecomputeStatus(true);
    else
      linearSystem->setRecomputeStatus(false);
    // FESOLVEBDF2NL
    initializeBDF1(sol, metaNumber, mesh, solBDF1);
    printf("\n");
    printf("Étape 1 - recomputeMatrix = %s : Solution BDF1 - t = %6.6e\n",
           linearSystem->getRecomputeStatus() ? "true" : "false", sol->getCurrentTime());
    solveQNBDF(solBDF1, tol, metaNumber, linearSystem, sol, mesh);
    fePstClc(sol, linearSystem, solBDF1);
    // Compute L2 norm of the solution
    sol->setSolFromContainer(solBDF1);
    // norm->computeL2Norm(metaNumber, sol, mesh);
    // normL2[iTime] = norm->getNorm();
    // int cnt = 0;
    for(auto *norm : norms) {
      norm->computeL2Norm(metaNumber, sol, mesh);
      normL2[iTime] = norm->getNorm(); // TODO : fix this
      // std::cout<<"Norme "<<cnt++<<" = "<<normL2[iTime]<<std::endl;
    }

    // std::string vtkFile = "../../data/BDF1__" + std::to_string(iTime) + ".vtk";
    // feExporterVTK writer(vtkFile, mesh, sol, metaNumber, spaces);
    // int NumberOfDofs=81;
    // getNDofsFromFile(vtkFile, NumberOfDofs);
    // std::cout<<"Nb de Dofs"<<NumberOfDofs<<std::endl;

    // std::vector<double> lastSol;
    // lastSol.resize(NumberOfDofs);
    // readLastSol(vtkFile, lastSol);
    // std::cout<<"le fichier en question est le"<<vtkFile<<std::endl;
    // std::cout<<"le nombre de Dofs est de "<<NumberOfDofs<<std::endl;
    // for(int i=0; i<NumberOfDofs;i++){
    //   std::cout<<"la valeur du champ BDF1 vaut"<<lastSol[i]<<std::endl;
    // }
  }

  delete solBDF1;
}

void solveDC3(std::vector<double> &normL2BDF2, std::vector<double> &normL2DC3, feTolerances tol,
              feMetaNumber *metaNumber, feLinearSystem *linearSystem,
              std::vector<feBilinearForm *> &formMatrices,
              std::vector<feBilinearForm *> &formResiduals, feSolution *sol,
              std::vector<feNorm *> &norms, feMesh *mesh, std::vector<feSpace *> &spaces) {
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

  std::vector<double> tK1K2(nSteps + 1, 0.0);
  for(int i = 0; i < nSteps + 1; ++i) tK1K2[i] = t0 + i * dt;
  for(int i = 0; i < nSteps; i += 2)
    if(i + 2 < nSteps + 1) tK1K2[i + 1] = tK1K2[i] + f * (tK1K2[i + 2] - tK1K2[i]);

  // Initialization
  if(K1K2) dt = tK1K2[1] - tK1K2[0];
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

  std::string vtkFile = "../../data/VK/cylindre" + std::to_string(0) + ".vtk";
  feExporterVTK writer(vtkFile, mesh, sol, metaNumber, spaces);
  // int cnt = 0;
  for(auto *norm : norms) {
    norm->computeL2Norm(metaNumber, sol, mesh);
    normL2BDF2[0] = norm->getNorm();
    normL2DC3[0] = norm->getNorm();
    // normL2DC3[0+cnt] = norm->getNorm();
    // normL2[iTime] = norm->getNorm(); // TODO : fix this
    // std::cout<<"Norme "<<cnt++<<" = "<<normL2DC3[0+cnt]<<std::endl;
  }

  // Start integration with a BDF2
  for(int iTime = 1; iTime <= 2; ++iTime) {
    if(K1K2) dt = tK1K2[iTime + 1] - tK1K2[iTime];
    solBDF2->rotate(dt);
    solDC3->rotate(dt);
    solDC3->initialize(sol, mesh, metaNumber);
    if(iTime == 1)
      linearSystem->setRecomputeStatus(true);
    else
      linearSystem->setRecomputeStatus(false);
    // FESOLVEBDF2NL
    initializeBDF2(sol, metaNumber, mesh, solBDF2);
    printf("\n");
    printf("Étape 1 - recomputeMatrix = %s : Solution BDF2 - t = %6.6e\n",
           linearSystem->getRecomputeStatus() ? "true" : "false", sol->getCurrentTime());
    solveQNBDF(solBDF2, tol, metaNumber, linearSystem, sol, mesh);
    fePstClc(sol, linearSystem, solBDF2);
    // Compute L2 norm of the solution
    sol->setSolFromContainer(solBDF2);
    // norm->computeL2Norm(metaNumber, sol, mesh);
    // normL2BDF2[iTime] = norm->getNorm();
    sol->setSolFromContainer(solDC3);
    // norm->computeL2Norm(metaNumber, sol, mesh);
    // normL2DC3[iTime] = norm->getNorm();

    std::string vtkFile = "../../data/VK/cylindre" + std::to_string(iTime) + ".vtk";
    feExporterVTK writer(vtkFile, mesh, sol, metaNumber, spaces);

    // int cnt = 0;
    for(auto *norm : norms) {
      norm->computeL2Norm(metaNumber, sol, mesh);
      normL2BDF2[iTime] = norm->getNorm();
      normL2DC3[iTime] = norm->getNorm();
      // normL2DC3[3*iTime+cnt] = norm->getNorm();
      // normL2[iTime] = norm->getNorm(); // TODO : fix this
      // std::cout<<"Norme "<<cnt++<<" = "<<normL2DC3[3*iTime+cnt]<<std::endl;
    }
  }
  // Continue with DC3
  for(int iTime = 3; iTime < nSteps; ++iTime) {
    // for(int iTime = 3; iTime < 4; ++iTime){
    if(K1K2) dt = tK1K2[iTime + 1] - tK1K2[iTime];
    solBDF2->rotate(dt);
    solDC3->rotate(dt);
    // if(iTime==3)
    //   linearSystem->setRecomputeStatus(true);
    // else
    //   linearSystem->setRecomputeStatus(false);
    // FESOLVEDC3F
    initializeBDF2(sol, metaNumber, mesh, solBDF2);
    printf("\n");
    printf("Étape 2 - recomputeMatrix = %s : Solution DC3F - t = %6.6e\n",
           linearSystem->getRecomputeStatus() ? "true" : "false", sol->getCurrentTime());
    solveQNBDF(solBDF2, tol, metaNumber, linearSystem, sol, mesh);
    fePstClc(sol, linearSystem, solBDF2);
    linearSystem->setRecomputeStatus(false);
    initializeDC3F(sol, metaNumber, mesh, solBDF2, solDC3);
    solveQNBDF(solDC3, tol, metaNumber, linearSystem, sol, mesh);
    fePstClc(sol, linearSystem, solDC3);
    // // Compute L2 norm of the solution
    sol->setSolFromContainer(solBDF2);
    // norm->computeL2Norm(metaNumber, sol, mesh);
    // normL2BDF2[iTime] = norm->getNorm();
    sol->setSolFromContainer(solDC3);
    // norm->computeL2Norm(metaNumber, sol, mesh);
    // normL2DC3[iTime] = norm->getNorm();

    std::string vtkFile = "../../data/VK/cylindre" + std::to_string(iTime) + ".vtk";
    feExporterVTK writer(vtkFile, mesh, sol, metaNumber, spaces);

    // int cnt = 0;
    for(auto *norm : norms) {
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