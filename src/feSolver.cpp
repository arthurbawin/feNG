#include "feSolver.h"
#include "feSolutionContainer.h"
#include "feExporter.h"
#include "feExporter.h"
#include "feNG.h"

double _normR0, _normFirstR0;

extern int FE_VERBOSE;

feStatus createTimeIntegrator(TimeIntegrator *&solver, timeIntegratorScheme scheme,
                              feTolerances tolerances, feLinearSystem *system,
                              feMetaNumber *metaNumber, feSolution *solution, feMesh *mesh,
                              std::vector<feNorm *> &norms, feExportData exportData, double tBegin,
                              double tEnd, int nTimeSteps, std::string initializationCode)
{
  switch(scheme) {
    case STATIONARY:
      solver =
        new StationarySolver(tolerances, metaNumber, system, solution, norms, mesh, exportData);
      break;
    case BDF1:
      solver = new BDF1Solver(tolerances, metaNumber, system, solution, norms, mesh, exportData,
                              tBegin, tEnd, nTimeSteps);
      break;
    case BDF2:
      solver = new BDF2Solver(tolerances, metaNumber, system, solution, norms, mesh, exportData,
                              tBegin, tEnd, nTimeSteps, initializationCode);
      break;
    case DC2F:
      solver = new DC2FSolver(tolerances, metaNumber, system, solution, norms, mesh, exportData,
                              tBegin, tEnd, nTimeSteps);
      break;
    case DC3:
      solver = new DC3Solver(tolerances, metaNumber, system, solution, norms, mesh, exportData,
                             tBegin, tEnd, nTimeSteps, initializationCode);
      break;
    case DC3F:
      solver = new DC3FSolver(tolerances, metaNumber, system, solution, norms, mesh, exportData,
                              tBegin, tEnd, nTimeSteps);
      break;
    default:
      return feErrorMsg(FE_STATUS_ERROR, "Unsupported time integration scheme.");
  }
  return FE_STATUS_OK;
}

void solveQNBDF(feSolutionContainer *solDot, feTolerances tol, feMetaNumber *metaNumber,
                feLinearSystem *linearSystem, feSolution *sol, feMesh *mesh)
{
  bool newton = true;
  bool status = linearSystem->getRecomputeStatus();
  bool prev_status = status;
  int iter = 0, linearSystemIter;
  double normDx, normResidual, normAxb;
  int cnt = 0; // nb of iteration without recompute the matrix
  // std::cout<<"linearSystem == "<<status<<std::endl;
  // if(status == false) prev_status =
  while(newton) {
    // linearSystem->setRecomputeStatus(false); //By default, we don't reassemble the matrix
    // tic();
    linearSystem->setToZero();
    solDot->computeSolTimeDerivative(sol, linearSystem);
    linearSystem->assemble(sol);
    linearSystem->solve(&normDx, &normResidual, &normAxb, &linearSystemIter);
    if(iter == 0 && status == 1) _normR0 = normResidual;
    if(iter == 0 && status == 0)
      _normFirstR0 = normResidual; // in case we have to recalculate the matrix and update the R0
    // FECORRECTIONSOLUTION
    linearSystem->correctSolution(sol);
    solDot->setSol(0, sol->getSolutionCopy());

    feInfoCond(
      FE_VERBOSE > 0,
      "    iter %2d : ||A*dx-res|| = %10.10e (%4d iter.) \t ||dx|| = %10.10e \t ||res|| = %10.10e",
      ++iter, normAxb, linearSystemIter, normDx, normResidual);

    newton = !((normDx <= tol.tolDx && normResidual <= tol.tolResidual) || iter > tol.maxIter);

    // linearSystem->setRecomputeStatus(true);
    // std::cout<<"R0 == "<<_normR0<<std::endl;
    // std::cout<<"normResidual / _normR0 == "<<normResidual / _normR0<<std::endl;
    // std::cout<<"cnt == "<<cnt<<std::endl;
    if((iter > 2 || prev_status == 0) && ((normResidual / _normR0) < 0.001) && cnt < 8) {
      linearSystem->setRecomputeStatus(false);

      // std::cout<<"On garde la matrice"<<std::endl;
    } else if(iter == 1 && prev_status == 0)
      linearSystem->setRecomputeStatus(false); // we decide to try with the old matrix at the first
                                               // iteration just to see if it's enough
    else {
      if(iter == 2 && prev_status == 0)
        _normR0 = _normFirstR0; // if we have to recompute the matrix we update the R0
      linearSystem->setRecomputeStatus(true);
      // std::cout<<"On recalcule la matrice"<<std::endl;
      cnt = 0;
    }
    cnt++;
  }

  if(iter > tol.maxIter) {
    feWarning("Not converged at iter %2d : ||A*dx-res|| = %10.10e  (%4d iter.) \t ||dx|| = "
              "%10.10e \t ||res|| = %10.10e",
              iter, normAxb, linearSystemIter, normDx, normResidual);
  } else {
    feInfoCond(FE_VERBOSE > 0,
               "Converged in %2d Newton iterations : ||dx|| = %10.10e \t ||res|| = %10.10e", iter,
               normDx, normResidual);
  }
  linearSystem->setRecomputeStatus(status);
}

// Deprecated
void solveStationary(double *normL2, feTolerances tol, feMetaNumber *metaNumber,
                     feLinearSystem *linearSystem, std::vector<feBilinearForm *> &formMatrices,
                     std::vector<feBilinearForm *> &formResiduals, feSolution *sol,
                     std::vector<feNorm *> &norms, feMesh *mesh)
{
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
                                   std::vector<feNorm *> &norms, feMesh *mesh,
                                   feExportData exportData)
  : TimeIntegrator(tol, metaNumber, linearSystem, sol, norms, mesh, exportData)
{
  // Initialize the solution container with a single solution
  int nSol = 1;
  _solutionContainer = new feStationarySolution(nSol, _sol->getCurrentTime(), _metaNumber);
  _solutionContainer->initialize(_sol, _mesh, _metaNumber);

  _normL2.resize(norms.size());
  for(auto &n : _normL2) n.resize(1, 0.);
}

feStatus StationarySolver::makeSteps(int nSteps)
{
  feInfoCond(FE_VERBOSE > 0, "Solving for steady state solution - recomputeMatrix = %s",
             _linearSystem->getRecomputeStatus() ? "true" : "false");
  _linearSystem->setRecomputeStatus(true);
  _sol->initializeUnknowns(_mesh, _metaNumber);
  _sol->initializeEssentialBC(_mesh, _metaNumber);
  // Solve
  solveQNBDF(_solutionContainer, _tol, _metaNumber, _linearSystem, _sol, _mesh);
  // Compute L2 norm of solution(s)
  for(int k = 0; k < _norms.size(); ++k) {
    if(_norms[k]->getNbFields() > 1) {
      if(_norms[k]->getTypeNorm() == "NormH1") {
        _norms[k]->computeH1NormVec(_metaNumber, _sol, _mesh);
        _normL2[k][0] = _norms[k]->getNorm();
      } else if(_norms[k]->getTypeNorm() == "NormSemiH1") {
        _norms[k]->computeSemiH1NormVec(_metaNumber, _sol, _mesh);
        _normL2[k][0] = _norms[k]->getNorm();
      } else if(_norms[k]->getTypeNorm() == "NormL2") {
        _norms[k]->computeL2NormVec(_metaNumber, _sol, _mesh);
        _normL2[k][0] = _norms[k]->getNorm();
      } else {
        printf("%s\n", "Norm doesn't exist");
      }

    } else {
      if(_norms[k]->getTypeNorm() == "IntFluxNx") {
        _norms[k]->computeIntFluxNx(_metaNumber, _sol, _mesh);
        _normL2[k][0] = _norms[k]->getNorm();
      } else if(_norms[k]->getTypeNorm() == "IntFluxNy") {
        _norms[k]->computeIntFluxNy(_metaNumber, _sol, _mesh);
        _normL2[k][0] = _norms[k]->getNorm();
      } else {
        _norms[k]->computeL2Norm(_metaNumber, _sol, _mesh);
        _normL2[k][0] = _norms[k]->getNorm();
      }
    }
  }
  // Export for visualization
  if(_exportData.exporter != nullptr) {
    std::string fileName = _exportData.fileNameRoot + "Stationary.vtk";
    feCheck(_exportData.exporter->writeStep(fileName));
  }

  return FE_STATUS_OK;
}

/* Nom cryptique :
  - reset sol->soldot
  - recalcule le residu
  - assigne le residu a solContainer->F */
void fePstClc(feSolution *sol, feLinearSystem *linearSystem, feSolutionContainer *solContainer)
{
  sol->setSolDotToZero();
  linearSystem->setResidualToZero();
  linearSystem->assembleResiduals(sol);
  linearSystem->assignResidualToDCResidual(solContainer);
}

// Parameter of varaible time step
double _f = 0.20; //_f =0.2 means dt1/dt2 = 4  _f=0.25 means dt1/dt2 = 3
bool K1K2 = false;

BDF2Solver::BDF2Solver(feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem,
                       feSolution *sol, std::vector<feNorm *> &norms, feMesh *mesh,
                       feExportData exportData, double t0, double tEnd, int nTimeSteps,
                       std::string CodeIni)
  : TimeIntegrator(tol, metaNumber, linearSystem, sol, norms, mesh, exportData, t0, tEnd,
                   nTimeSteps, CodeIni)
{
  // Initialize the solution container
  int nSol = 3;
  _solutionContainer = new feSolutionBDF2(nSol, _sol->getCurrentTime(), _metaNumber);
  _solutionContainer->initialize(_sol, _mesh, _metaNumber);

  _normL2.resize(norms.size());
  for(auto &n : _normL2) n.resize(_nTimeSteps, 0.);

  // printf("Initializing BDF2 solver : integrating from t0 = %f to tEnd = %f in %d steps\n", _t0,
  //        _tEnd, _nTimeSteps);
}

feStatus BDF2Solver::makeSteps(int nSteps)
{
  printf("BDF2 : Advancing %d steps from t = %f to t = %f\n", nSteps, _tCurrent,
         _tCurrent + nSteps * _dt);

  if(_currentStep == 0) {
    // Initialization and first step

    if(_CodeIni == "BDF1/DCF") {
      printf("\n");
      printf(" ----------------------------- \n");
      printf("Using DC2F to initialize BDF2 \n");
      printf(" ----------------------------- \n");
      printf("\n");
      double _t_ini = _t0 + 2 * _dt;
      std::vector<feNorm *> norms2 = {
        _norms[0], _norms[0]}; // Il faut un vecteur norme de taille 2 sinon pb //!!!A changer!!!
      DC2FSolver solver(_tol, _metaNumber, _linearSystem, _sol, norms2, _mesh, _exportData, _t0,
                        _t_ini, 2);
      solver.makeSteps(2);
      feSolutionContainer *_solutionContainerDC2F = solver.getSolutionContainer();
      _solutionContainer->setSol(0, _solutionContainerDC2F->getSolution(0));
      _solutionContainer->setSol(1, _solutionContainerDC2F->getSolution(1));
      _solutionContainer->setSol(2, _solutionContainerDC2F->getSolution(2));
      delete _solutionContainerDC2F;
      _solutionContainer->rotate(_dt);
    } else {
      _solutionContainer->rotate(_dt);
      _sol->setSolFromContainer(_solutionContainer);
      _solutionContainer->initialize(_sol, _mesh, _metaNumber);
    }

    _linearSystem->setRecomputeStatus(true);
    _sol->setSolFromContainer(_solutionContainer);
    _tCurrent += _dt;
    ++_currentStep;

    for(size_t i = 0; i < _norms.size(); ++i) {
      _norms[i]->computeL2Norm(_metaNumber, _sol, _mesh);
      _normL2[i][0] = _norms[i]->getNorm();
    }

    feInfoCond(FE_VERBOSE > 0, "Current step = %d/%d : t = %f\n", _currentStep, nSteps, _tCurrent);

    if(_exportData.exporter != nullptr && (_currentStep % _exportData.exportEveryNSteps) == 0) {
      std::string fileName = _exportData.fileNameRoot + std::to_string(_currentStep) + ".vtk";
      feCheck(_exportData.exporter->writeStep(fileName));
    }

    --nSteps; // To advance the same number of steps than if currentStep != 0
  }

  std::vector<double> tK1K2(nSteps + 1, 0.0);
  for(int i = 0; i < nSteps + 1; ++i) tK1K2[i] = _t0 + i * _dt;
  for(int i = 0; i < nSteps; i += 2)
    if(i + 2 < nSteps + 1) tK1K2[i + 1] = tK1K2[i] + _f * (tK1K2[i + 2] - tK1K2[i]);

  for(int i = 0; i < nSteps; ++i) {
    if(i == 0) {
      _linearSystem->setRecomputeStatus(true);
      if(_CodeIni == "BDF1/DCF") {
        _solutionContainer->setSol(1, _solutionContainer->getSolution(2));
        printf("\n");
        printf(" ----------------------------- \n");
        printf("Continue with BDF2 \n");
        printf(" ----------------------------- \n");
        printf("\n");
      }
    } else {
      _linearSystem->setRecomputeStatus(false);
    }
    _solutionContainer->rotate(_dt);
    initializeBDF2(_sol, _metaNumber, _mesh, dynamic_cast<feSolutionBDF2 *>(_solutionContainer));
    printf("\n");
    printf("Étape 1 - recomputeMatrix = %s : Solution BDF2 - t = %6.6e\n",
           _linearSystem->getRecomputeStatus() ? "true" : "false", _sol->getCurrentTime());
    solveQNBDF(_solutionContainer, _tol, _metaNumber, _linearSystem, _sol, _mesh);
    fePstClc(_sol, _linearSystem, _solutionContainer);

    // Calcul of norms
    for(int k = 0; k < _norms.size(); ++k) {
      if(_norms[k]->getNbFields() > 1) {
        if(_norms[k]->getTypeNorm() == "NormH1") {
          _norms[k]->computeH1NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "NormSemiH1") {
          _norms[k]->computeSemiH1NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "NormL2") {
          _norms[k]->computeL2NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else {
          printf("%s\n", "Norm doesn't exist");
        }

      } else {
        if(_norms[k]->getTypeNorm() == "IntFlux") {
          _norms[k]->computeNormLambda(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "IntForce") {
          _norms[k]->computeNormLambdaNS(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
          std::cout << "Norm L-DC3" << _normL2[k][_currentStep] << std::endl;
        } else if(_norms[k]->getTypeNorm() == "IntLambda") {
          _norms[k]->computeIntegralNum(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else {
          _norms[k]->computeL2Norm(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        }
      }
    }

    _tCurrent += _dt;
    ++_currentStep;
    if(K1K2) _dt = tK1K2[i + 1] - tK1K2[i];
    printf("\n");
    printf("Current step = %d : t = %f\n", _currentStep, _tCurrent);

    if(_exportData.exporter != nullptr && (_currentStep % _exportData.exportEveryNSteps) == 0) {
      std::string fileName = _exportData.fileNameRoot + std::to_string(_currentStep) + ".vtk";
      _exportData.exporter->writeStep(fileName);
    }
  }

  return FE_STATUS_OK;
}

BDF1Solver::BDF1Solver(feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem,
                       feSolution *sol, std::vector<feNorm *> &norms, feMesh *mesh,
                       feExportData exportData, double t0, double tEnd, int nTimeSteps)
  : TimeIntegrator(tol, metaNumber, linearSystem, sol, norms, mesh, exportData, t0, tEnd,
                   nTimeSteps)
{
  // Initialize the solution container
  int nSol = 2;
  _solutionContainer = new feSolutionBDF1(nSol, _sol->getCurrentTime(), _metaNumber);
  _solutionContainer->initialize(_sol, _mesh, _metaNumber);

  _normL2.resize(norms.size());
  for(auto &n : _normL2) n.resize(_nTimeSteps, 0.);

  printf("Initializing BDF1 solver : integrating from t0 = %f to tEnd = %f in %d steps\n", _t0,
         _tEnd, _nTimeSteps);
}

feStatus BDF1Solver::makeSteps(int nSteps)
{
  printf("BDF1 : Advancing %d steps from t = %f to t = %f\n", nSteps, _tCurrent,
         _tCurrent + nSteps * _dt);

  if(_currentStep == 0) {
    // Initialization and first step
    _linearSystem->setRecomputeStatus(true);
    _solutionContainer->initialize(_sol, _mesh, _metaNumber);
    _sol->setSolFromContainer(_solutionContainer);

    for(size_t i = 0; i < _norms.size(); ++i) {
      _norms[i]->computeL2Norm(_metaNumber, _sol, _mesh);
      _normL2[i][_currentStep] = _norms[i]->getNorm();
    }

    printf("Current step = %d/%d : t = %f\n", _currentStep, nSteps, _tCurrent);
  }

  std::vector<double> tK1K2(nSteps + 1, 0.0);
  for(int i = 0; i < nSteps + 1; ++i) tK1K2[i] = _t0 + i * _dt;
  for(int i = 0; i < nSteps; i += 2)
    if(i + 2 < nSteps + 1) tK1K2[i + 1] = tK1K2[i] + _f * (tK1K2[i + 2] - tK1K2[i]);

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
    _sol->setSolFromContainer(_solutionContainer);
    // Calcul of norms
    for(int k = 0; k < _norms.size(); ++k) {
      if(_norms[k]->getNbFields() > 1) {
        if(_norms[k]->getTypeNorm() == "NormH1") {
          _norms[k]->computeH1NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "NormSemiH1") {
          _norms[k]->computeSemiH1NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "NormL2") {
          _norms[k]->computeL2NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else {
          printf("%s\n", "Norm doesn't exist");
        }

      } else {
        if(_norms[k]->getTypeNorm() == "IntFlux") {
          _norms[k]->computeNormLambda(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "IntForce") {
          _norms[k]->computeNormLambdaNS(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
          std::cout << "Norm L-DC3" << _normL2[k][_currentStep] << std::endl;
        } else if(_norms[k]->getTypeNorm() == "IntLambda") {
          _norms[k]->computeIntegralNum(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else {
          _norms[k]->computeL2Norm(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        }
      }
    }

    _tCurrent += _dt;
    ++_currentStep;
    if(K1K2) _dt = tK1K2[i + 1] - tK1K2[i];
    printf("Current step = %d : t = %f\n", _currentStep, _tCurrent);

    // std::string vtkFile =
    //   "../../data/VTK_VonKarman/dataTestV3" + std::to_string(_currentStep) + ".vtk";
    // // std::string vtkFile = "../../data/VTK_VonKarman/dataTestBullShit" +
    // // std::to_string(_currentStep) + ".vtk";
    // feExporterVTK writer(vtkFile, _mesh, _sol, _metaNumber, spaces);
  }

  return FE_STATUS_OK;
}

DC2FSolver::DC2FSolver(feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem,
                       feSolution *sol, std::vector<feNorm *> &norms, feMesh *mesh,
                       feExportData exportData, double t0, double tEnd, int nTimeSteps)
  : TimeIntegrator(tol, metaNumber, linearSystem, sol, norms, mesh, exportData, t0, tEnd,
                   nTimeSteps)
{
  // Initialize the solution container
  int nSol = 3;
  _solutionContainerBDF1 = new feSolutionBDF1(nSol, _sol->getCurrentTime(), _metaNumber);
  _solutionContainer = new feSolutionDC2F(nSol, _sol->getCurrentTime(), _metaNumber);
  _solutionContainerBDF1->initialize(_sol, _mesh, _metaNumber);
  _solutionContainer->initialize(_sol, _mesh, _metaNumber);
  fePstClc(_sol, _linearSystem, _solutionContainerBDF1);

  _normL2.resize(norms.size());
  for(auto &n : _normL2) n.resize(_nTimeSteps, 0.);

  printf("Initializing DC2F solver with BDF1 : integrating from t0 = %f to tEnd = %f in %d steps\n",
         _t0, _tEnd, _nTimeSteps);
}

feStatus DC2FSolver::makeSteps(int nSteps)
{
  printf("DC2 : Advancing %d steps from t = %f to t = %f\n", nSteps, _tCurrent,
         _tCurrent + nSteps * _dt);

  if(_currentStep == 0) {
    _linearSystem->setRecomputeStatus(true);
    _solutionContainerBDF1->initialize(_sol, _mesh, _metaNumber);
    _sol->setSolFromContainer(_solutionContainerBDF1);
    // _solutionContainerBDF1->rotate(_dt);

    _solutionContainer->initialize(_sol, _mesh, _metaNumber);
    _sol->setSolFromContainer(_solutionContainer);

    for(size_t i = 0; i < _norms.size(); ++i) {
      _norms[i]->computeL2Norm(_metaNumber, _sol, _mesh);
      _normL2[i][0] = _norms[i]->getNorm();
    }
    std::cout << "========dt vaux ==========" << _dt << std::endl;
    printf("\n");
    printf("Current step = %d/%d : t = %f\n", _currentStep, nSteps, _tCurrent);
  }

  std::vector<double> tK1K2(nSteps + 1, 0.0);
  for(int i = 0; i < nSteps + 1; ++i) tK1K2[i] = _t0 + i * _dt;
  for(int i = 0; i < nSteps; i += 2)
    if(i + 2 < nSteps + 1) tK1K2[i + 1] = tK1K2[i] + _f * (tK1K2[i + 2] - tK1K2[i]);

  // Integration with BDF1 and DC2F
  for(int i = 0; i < nSteps; ++i) {
    _solutionContainerBDF1->rotate(_dt);
    _solutionContainer->rotate(_dt);
    initializeBDF1(_sol, _metaNumber, _mesh,
                   dynamic_cast<feSolutionBDF1 *>(_solutionContainerBDF1));
    printf("\n");
    printf("Étape 1 - recomputeMatrix = %s : Solution BDF1 - t = %6.6e\n",
           _linearSystem->getRecomputeStatus() ? "true" : "false", _sol->getCurrentTime());
    solveQNBDF(_solutionContainerBDF1, _tol, _metaNumber, _linearSystem, _sol, _mesh);
    fePstClc(_sol, _linearSystem, _solutionContainerBDF1);
    _sol->setSolFromContainer(_solutionContainerBDF1);
    // Compute L2 norm of BDF1 solution
    _norms[0]->computeL2Norm(_metaNumber, _sol, _mesh);
    _normL2[0][_currentStep] = _norms[0]->getNorm();
    if(_norms.size() > 2) {
      _norms[2]->computeL2Norm(_metaNumber, _sol, _mesh);
      _normL2[2][_currentStep] = _norms[2]->getNorm();
      _norms[3]->computeL2Norm(_metaNumber, _sol, _mesh);
      _normL2[3][_currentStep] = _norms[3]->getNorm();
    }
    initializeDC2F(_sol, _metaNumber, _mesh, dynamic_cast<feSolutionBDF1 *>(_solutionContainerBDF1),
                   dynamic_cast<feSolutionDC2F *>(_solutionContainer));
    printf("\n");
    printf("Étape 2 - recomputeMatrix = %s : Solution DC2F - t = %6.6e\n",
           _linearSystem->getRecomputeStatus() ? "true" : "false", _sol->getCurrentTime());
    solveQNBDF(_solutionContainer, _tol, _metaNumber, _linearSystem, _sol, _mesh);
    fePstClc(_sol, _linearSystem, _solutionContainer);
    // Compute L2 norm of DC2F solution
    _norms[1]->computeL2Norm(_metaNumber, _sol, _mesh);
    _normL2[1][_currentStep] = _norms[1]->getNorm();
    if(_norms.size() > 3) {
      _norms[4]->computeL2Norm(_metaNumber, _sol, _mesh);
      _normL2[4][_currentStep] = _norms[4]->getNorm();
      _norms[5]->computeL2Norm(_metaNumber, _sol, _mesh);
      _normL2[5][_currentStep] = _norms[5]->getNorm();
    }
    _sol->setSolFromContainer(_solutionContainer);
    _tCurrent += _dt;
    ++_currentStep;
    if(K1K2) _dt = tK1K2[i + 1] - tK1K2[i];
    std::cout << "========dt vaux ==========" << _dt << std::endl;
    printf("\n");
    printf("Current step = %d : t = %f\n", _currentStep, _tCurrent);

    // std::string vtkFile = "../../data/cylindreAdapt" + std::to_string(_currentStep) + ".vtk";
    // feExporterVTK writer(vtkFile, _mesh, _sol, _metaNumber, spaces);
  }

  return FE_STATUS_OK;
}

DC3FSolver::DC3FSolver(feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem,
                       feSolution *sol, std::vector<feNorm *> &norms, feMesh *mesh,
                       feExportData exportData, double t0, double tEnd, int nTimeSteps)
  : TimeIntegrator(tol, metaNumber, linearSystem, sol, norms, mesh, exportData, t0, tEnd,
                   nTimeSteps)
{
  // Initialize the solution container
  int nSol = 6;
  _solutionContainerBDF1 = new feSolutionBDF1(nSol, _sol->getCurrentTime(), _metaNumber);
  _solutionContainerDC2F = new feSolutionDC2F(nSol, _sol->getCurrentTime(), _metaNumber);
  _solutionContainer = new feSolutionDC2F(nSol, _sol->getCurrentTime(), _metaNumber);
  _solutionContainerBDF1->initialize(_sol, _mesh, _metaNumber);
  _solutionContainerDC2F->initialize(_sol, _mesh, _metaNumber);
  _solutionContainer->initialize(_sol, _mesh, _metaNumber);
  fePstClc(_sol, _linearSystem, _solutionContainerBDF1);
  fePstClc(_sol, _linearSystem, _solutionContainerDC2F);

  _normL2.resize(norms.size());
  for(auto &n : _normL2) n.resize(_nTimeSteps, 0.);

  printf("Initializing DC3F solver with BDF1 and DC2F : integrating from t0 = %f to tEnd = %f in "
         "%d steps\n",
         _t0, _tEnd, _nTimeSteps);
}

feStatus DC3FSolver::makeSteps(int nSteps)
{
  printf("DC2 : Advancing %d steps from t = %f to t = %f\n", nSteps, _tCurrent,
         _tCurrent + nSteps * _dt);

  if(_currentStep == 0) {
    _linearSystem->setRecomputeStatus(true);
    _solutionContainerBDF1->initialize(_sol, _mesh, _metaNumber);
    _solutionContainerBDF1->rotate(_dt);
    _sol->setSolFromContainer(_solutionContainerBDF1);
    _solutionContainerDC2F->initialize(_sol, _mesh, _metaNumber);
    _solutionContainerDC2F->rotate(_dt);
    _sol->setSolFromContainer(_solutionContainerDC2F);
    _solutionContainer->initialize(_sol, _mesh, _metaNumber);
    _solutionContainer->rotate(_dt);
    _sol->setSolFromContainer(_solutionContainer);

    for(size_t i = 0; i < _norms.size(); ++i) {
      // _norms[i]->computeL2Norm0D( _sol);
      _norms[i]->computeL2Norm(_metaNumber, _sol, _mesh);
      _normL2[i][0] = _norms[i]->getNorm();
    }
    printf("-----------------------------------------------------------------");
    printf("\n");
    printf("Current step = %d/%d : t = %f\n", _currentStep, nSteps, _tCurrent);
  }

  std::vector<double> tK1K2(nSteps + 1, 0.0);
  for(int i = 0; i < nSteps + 1; ++i) tK1K2[i] = _t0 + i * _dt;
  for(int i = 0; i < nSteps; i += 2)
    if(i + 2 < nSteps + 1) tK1K2[i + 1] = tK1K2[i] + _f * (tK1K2[i + 2] - tK1K2[i]);

  // Continue with DC2F
  for(int i = 0; i < 1; ++i) {
    initializeBDF1(_sol, _metaNumber, _mesh,
                   dynamic_cast<feSolutionBDF1 *>(_solutionContainerBDF1));
    printf("\n");
    printf("Étape 1 - recomputeMatrix = %s : Solution BDF1 - t = %6.6e\n",
           _linearSystem->getRecomputeStatus() ? "true" : "false", _sol->getCurrentTime());
    solveQNBDF(_solutionContainerBDF1, _tol, _metaNumber, _linearSystem, _sol, _mesh);
    fePstClc(_sol, _linearSystem, _solutionContainerBDF1);
    // Compute L2 norm of BDF1 solution
    for(int k = 0; k < _norms.size() - 2; k += 3) {
      if(_norms[k]->getNbFields() > 1) {
        if(_norms[k]->getTypeNorm() == "NormH1") {
          _norms[k]->computeH1NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "NormSemiH1") {
          _norms[k]->computeSemiH1NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "NormL2") {
          _norms[k]->computeL2NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else {
          printf("%s\n", "Norm doesn't exist");
        }

      } else {
        _norms[k]->computeL2Norm(_metaNumber, _sol, _mesh);
        _normL2[k][_currentStep] = _norms[k]->getNorm();
      }
    }
    initializeDC2F(_sol, _metaNumber, _mesh, dynamic_cast<feSolutionBDF1 *>(_solutionContainerBDF1),
                   dynamic_cast<feSolutionDC2F *>(_solutionContainerDC2F));
    printf("\n");
    printf("Étape 2 - recomputeMatrix = %s : Solution DC2F - t = %6.6e\n",
           _linearSystem->getRecomputeStatus() ? "true" : "false", _sol->getCurrentTime());
    solveQNBDF(_solutionContainerDC2F, _tol, _metaNumber, _linearSystem, _sol, _mesh);
    fePstClc(_sol, _linearSystem, _solutionContainerDC2F);
    // Compute L2 norm of DC2F solution
    for(int k = 1; k < _norms.size() - 1; k += 3) {
      if(_norms[k]->getNbFields() > 1) {
        if(_norms[k]->getTypeNorm() == "NormH1") {
          _norms[k]->computeH1NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "NormSemiH1") {
          _norms[k]->computeSemiH1NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "NormL2") {
          _norms[k]->computeL2NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
          std::cout << "Norm DC2F in ini = " << _normL2[k][_currentStep] << std::endl;
        } else {
          printf("%s\n", "Norm doesn't exist");
        }

      } else {
        _norms[k]->computeL2Norm(_metaNumber, _sol, _mesh);
        _normL2[k][_currentStep] = _norms[k]->getNorm();
      }
    }

    _sol->setSolFromContainer(_solutionContainerBDF1);
    _sol->setSolFromContainer(_solutionContainerDC2F);
    _tCurrent += _dt;
    ++_currentStep;
    printf("\n");
    printf("Current step = %d : t = %f\n", _currentStep, _tCurrent);
    // std::string vtkFile = "../../data/cylindreAdapt" + std::to_string(_currentStep) + ".vtk";
    // feExporterVTK writer(vtkFile, _mesh, _sol, _metaNumber, spaces);
  }

  // Continue with DC3F
  for(int i = 1; i < nSteps; ++i) {
    _solutionContainerBDF1->rotate(_dt);
    _solutionContainerDC2F->rotate(_dt);
    initializeBDF1(_sol, _metaNumber, _mesh,
                   dynamic_cast<feSolutionBDF1 *>(_solutionContainerBDF1));
    printf("\n");
    printf("Étape 1 - recomputeMatrix = %s : Solution BDF1 - t = %6.6e\n",
           _linearSystem->getRecomputeStatus() ? "true" : "false", _sol->getCurrentTime());
    solveQNBDF(_solutionContainerBDF1, _tol, _metaNumber, _linearSystem, _sol, _mesh);
    fePstClc(_sol, _linearSystem, _solutionContainerBDF1);
    // Calcul of norms BDF1
    for(int k = 0; k < _norms.size() - 2; k += 3) {
      if(_norms[k]->getNbFields() > 1) {
        if(_norms[k]->getTypeNorm() == "NormH1") {
          _norms[k]->computeH1NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "NormSemiH1") {
          _norms[k]->computeSemiH1NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "NormL2") {
          _norms[k]->computeL2NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else {
          printf("%s\n", "Norm doesn't exist");
        }

      } else {
        _norms[k]->computeL2Norm(_metaNumber, _sol, _mesh);
        _normL2[k][_currentStep] = _norms[k]->getNorm();
      }
    }
    initializeDC2F(_sol, _metaNumber, _mesh, dynamic_cast<feSolutionBDF1 *>(_solutionContainerBDF1),
                   dynamic_cast<feSolutionDC2F *>(_solutionContainerDC2F));
    printf("\n");
    printf("Étape 2 - recomputeMatrix = %s : Solution DC2F - t = %6.6e\n",
           _linearSystem->getRecomputeStatus() ? "true" : "false", _sol->getCurrentTime());
    solveQNBDF(_solutionContainerDC2F, _tol, _metaNumber, _linearSystem, _sol, _mesh);
    fePstClc(_sol, _linearSystem, _solutionContainerDC2F);

    // Calcul of norms DC2F
    for(int k = 1; k < _norms.size() - 1; k += 3) {
      if(_norms[k]->getNbFields() > 1) {
        if(_norms[k]->getTypeNorm() == "NormH1") {
          _norms[k]->computeH1NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "NormSemiH1") {
          _norms[k]->computeSemiH1NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "NormL2") {
          _norms[k]->computeL2NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
          std::cout << "Norm DC2F in ini = " << _normL2[k][_currentStep] << std::endl;
        } else {
          printf("%s\n", "Norm doesn't exist");
        }

      } else {
        std::cout << "on est la pour " << k << std::endl;
        _norms[k]->computeL2Norm(_metaNumber, _sol, _mesh);
        _normL2[k][_currentStep] = _norms[k]->getNorm();
      }
    }
    // at this moment we fixe the DC3F à the first time step

    if(i == 1) {
      initializeDC3F_centered(
        _sol, _metaNumber, _mesh, dynamic_cast<feSolutionDC2F *>(_solutionContainerDC2F),
        dynamic_cast<feSolutionDC2F *>(
          _solutionContainer)); // soit on creer un fonction soit on passe un paramêtre a
                                // la fonction car seul tn change en réalité
      printf("\n");
      printf("Retour a l'iteration precedante pour calculer t1 - recomputeMatrix = %s : Solution "
             "DC3F - t = %6.6e\n",
             _linearSystem->getRecomputeStatus() ? "true" : "false", _sol->getCurrentTime());
      solveQNBDF(_solutionContainer, _tol, _metaNumber, _linearSystem, _sol, _mesh);
      fePstClc(_sol, _linearSystem, _solutionContainer);
      // Calcul of norms DC3F
      for(int k = 2; k < _norms.size(); k += 3) {
        if(_norms[k]->getNbFields() > 1) {
          if(_norms[k]->getTypeNorm() == "NormH1") {
            _norms[k]->computeH1NormVec(_metaNumber, _sol, _mesh);
            _normL2[k][_currentStep - 1] = _norms[k]->getNorm();
          } else if(_norms[k]->getTypeNorm() == "NormSemiH1") {
            _norms[k]->computeSemiH1NormVec(_metaNumber, _sol, _mesh);
            _normL2[k][_currentStep - 1] = _norms[k]->getNorm();
          } else if(_norms[k]->getTypeNorm() == "NormL2") {
            _norms[k]->computeL2NormVec(_metaNumber, _sol, _mesh);
            _normL2[k][_currentStep - 1] = _norms[k]->getNorm();
          } else {
            printf("%s\n", "Norm doesn't exist");
          }

        } else {
          _norms[k]->computeL2Norm(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep - 1] = _norms[k]->getNorm();
        }
      }
    }
    _solutionContainer->rotate(_dt);
    initializeDC3F(_sol, _metaNumber, _mesh, dynamic_cast<feSolutionDC2F *>(_solutionContainerDC2F),
                   dynamic_cast<feSolutionDC2F *>(_solutionContainer));
    printf("\n");
    printf("Étape 3 - recomputeMatrix = %s : Solution DC3F - t = %6.6e\n",
           _linearSystem->getRecomputeStatus() ? "true" : "false", _sol->getCurrentTime());
    solveQNBDF(_solutionContainer, _tol, _metaNumber, _linearSystem, _sol, _mesh);
    fePstClc(_sol, _linearSystem, _solutionContainer);
    // Calcul of norms
    for(int k = 2; k < _norms.size(); k += 3) {
      if(_norms[k]->getNbFields() > 1) {
        if(_norms[k]->getTypeNorm() == "NormH1") {
          _norms[k]->computeH1NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "NormSemiH1") {
          _norms[k]->computeSemiH1NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "NormL2") {
          _norms[k]->computeL2NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else {
          printf("%s\n", "Norm doesn't exist");
        }

      } else {
        std::cout << "on est la pour " << k << std::endl;
        _norms[k]->computeL2Norm(_metaNumber, _sol, _mesh);
        _normL2[k][_currentStep] = _norms[k]->getNorm();
      }
    }
    // Compute L2 norm of the solution
    _sol->setSolFromContainer(_solutionContainerBDF1);
    _sol->setSolFromContainer(_solutionContainerDC2F);
    _sol->setSolFromContainer(_solutionContainer);
    _tCurrent += _dt;
    ++_currentStep;
    if(K1K2) _dt = tK1K2[i + 1] - tK1K2[i];
    std::cout << "========dt vaux ==========" << _dt << std::endl;
    printf("\n");
    printf("Current step = %d : t = %f\n", _currentStep, _tCurrent);

    // std::string vtkFile = "../../data/cylindreAdapt" + std::to_string(_currentStep) + ".vtk";
    // feExporterVTK writer(vtkFile, _mesh, _sol, _metaNumber, spaces);
  }
  return FE_STATUS_OK;
}

DC3FSolver_centered::DC3FSolver_centered(feTolerances tol, feMetaNumber *metaNumber,
                                         feLinearSystem *linearSystem, feSolution *sol,
                                         std::vector<feNorm *> &norms, feMesh *mesh,
                                         feExportData exportData, double t0, double tEnd,
                                         int nTimeSteps)
  : TimeIntegrator(tol, metaNumber, linearSystem, sol, norms, mesh, exportData, t0, tEnd,
                   nTimeSteps)
{
  // Initialize the solution container
  int nSol = 5;
  _solutionContainerBDF1 = new feSolutionBDF1(nSol, _sol->getCurrentTime(), _metaNumber);
  _solutionContainerDC2F = new feSolutionDC2F(nSol, _sol->getCurrentTime(), _metaNumber);
  _solutionContainer = new feSolutionDC2F(nSol, _sol->getCurrentTime(), _metaNumber);
  _solutionContainerBDF1->initialize(_sol, _mesh, _metaNumber);
  _solutionContainerDC2F->initialize(_sol, _mesh, _metaNumber);
  _solutionContainer->initialize(_sol, _mesh, _metaNumber);
  fePstClc(_sol, _linearSystem, _solutionContainerBDF1);
  fePstClc(_sol, _linearSystem, _solutionContainerDC2F);

  _normL2.resize(norms.size());
  for(auto &n : _normL2) n.resize(_nTimeSteps, 0.);

  printf("Initializing DC3F solver with BDF1 and DC2F : integrating from t0 = %f to tEnd = %f in "
         "%d steps\n",
         _t0, _tEnd, _nTimeSteps);
}

feStatus DC3FSolver_centered::makeSteps(int nSteps)
{
  printf("DC2 : Advancing %d steps from t = %f to t = %f\n", nSteps, _tCurrent,
         _tCurrent + nSteps * _dt);

  if(_currentStep == 0) {
    _linearSystem->setRecomputeStatus(true);
    _solutionContainerBDF1->initialize(_sol, _mesh, _metaNumber);
    _solutionContainerBDF1->rotate(_dt);
    _sol->setSolFromContainer(_solutionContainerBDF1);
    _solutionContainerDC2F->initialize(_sol, _mesh, _metaNumber);
    _solutionContainerDC2F->rotate(_dt);
    _sol->setSolFromContainer(_solutionContainerDC2F);
    _solutionContainer->initialize(_sol, _mesh, _metaNumber);
    _solutionContainer->rotate(_dt);
    _sol->setSolFromContainer(_solutionContainer);
    // ++_currentStep;

    for(size_t i = 0; i < _norms.size(); ++i) {
      _norms[i]->computeL2Norm(_metaNumber, _sol, _mesh);
      _normL2[i][0] = _norms[i]->getNorm();
    }
    printf("-----------------------------------------------------------------");
    printf("\n");
    printf("Current step = %d/%d : t = %f\n", _currentStep, nSteps, _tCurrent);
  }

  std::vector<double> tK1K2(nSteps + 1, 0.0);
  for(int i = 0; i < nSteps + 1; ++i) tK1K2[i] = _t0 + i * _dt;
  for(int i = 0; i < nSteps; i += 2)
    if(i + 2 < nSteps + 1) tK1K2[i + 1] = tK1K2[i] + _f * (tK1K2[i + 2] - tK1K2[i]);

  // Continue with DC2F
  for(int i = 0; i < 1; ++i) {
    initializeBDF1(_sol, _metaNumber, _mesh,
                   dynamic_cast<feSolutionBDF1 *>(_solutionContainerBDF1));
    printf("\n");
    printf("Étape 1 - recomputeMatrix = %s : Solution BDF1 - t = %6.6e\n",
           _linearSystem->getRecomputeStatus() ? "true" : "false", _sol->getCurrentTime());
    solveQNBDF(_solutionContainerBDF1, _tol, _metaNumber, _linearSystem, _sol, _mesh);
    fePstClc(_sol, _linearSystem, _solutionContainerBDF1);
    // Compute L2 norm of BDF1 solution
    _norms[0]->computeL2Norm(_metaNumber, _sol, _mesh);
    _normL2[0][_currentStep] = _norms[0]->getNorm();
    initializeDC2F(_sol, _metaNumber, _mesh, dynamic_cast<feSolutionBDF1 *>(_solutionContainerBDF1),
                   dynamic_cast<feSolutionDC2F *>(_solutionContainerDC2F));
    printf("\n");
    printf("Étape 2 - recomputeMatrix = %s : Solution DC2F - t = %6.6e\n",
           _linearSystem->getRecomputeStatus() ? "true" : "false", _sol->getCurrentTime());
    solveQNBDF(_solutionContainerDC2F, _tol, _metaNumber, _linearSystem, _sol, _mesh);
    fePstClc(_sol, _linearSystem, _solutionContainerDC2F);

    // Compute L2 norm of DC2F solution
    _norms[1]->computeL2Norm(_metaNumber, _sol, _mesh);
    _normL2[1][_currentStep] = _norms[1]->getNorm();

    _sol->setSolFromContainer(_solutionContainerBDF1);
    _sol->setSolFromContainer(_solutionContainerDC2F);
    _tCurrent += _dt;
    ++_currentStep;
    printf("\n");
    printf("Current step = %d : t = %f\n", _currentStep, _tCurrent);
    // std::string vtkFile = "../../data/cylindreAdapt" + std::to_string(_currentStep) + ".vtk";
    // feExporterVTK writer(vtkFile, _mesh, _sol, _metaNumber, spaces);
  }

  // Continue with DC3F
  for(int i = 1; i < nSteps; ++i) {
    _solutionContainerBDF1->rotate(_dt);
    _solutionContainerDC2F->rotate(_dt);
    initializeBDF1(_sol, _metaNumber, _mesh,
                   dynamic_cast<feSolutionBDF1 *>(_solutionContainerBDF1));
    printf("\n");
    printf("Étape 1 - recomputeMatrix = %s : Solution BDF1 - t = %6.6e\n",
           _linearSystem->getRecomputeStatus() ? "true" : "false", _sol->getCurrentTime());
    solveQNBDF(_solutionContainerBDF1, _tol, _metaNumber, _linearSystem, _sol, _mesh);
    fePstClc(_sol, _linearSystem, _solutionContainerBDF1);
    // Compute L2 norm of BDF1 solution
    _norms[0]->computeL2Norm(_metaNumber, _sol, _mesh);
    _normL2[0][_currentStep] = _norms[0]->getNorm();
    if(_norms.size() > 3) {
      _norms[3]->computeL2Norm(_metaNumber, _sol, _mesh);
      _normL2[3][_currentStep] = _norms[3]->getNorm();
      _norms[4]->computeL2Norm(_metaNumber, _sol, _mesh);
      _normL2[4][_currentStep] = _norms[4]->getNorm();
    }
    initializeDC2F(_sol, _metaNumber, _mesh, dynamic_cast<feSolutionBDF1 *>(_solutionContainerBDF1),
                   dynamic_cast<feSolutionDC2F *>(_solutionContainerDC2F));
    printf("\n");
    printf("Étape 2 - recomputeMatrix = %s : Solution DC2F - t = %6.6e\n",
           _linearSystem->getRecomputeStatus() ? "true" : "false", _sol->getCurrentTime());
    solveQNBDF(_solutionContainerDC2F, _tol, _metaNumber, _linearSystem, _sol, _mesh);
    fePstClc(_sol, _linearSystem, _solutionContainerDC2F);

    // Compute L2 norm of DC2F solution
    _norms[1]->computeL2Norm(_metaNumber, _sol, _mesh);
    _normL2[1][_currentStep] = _norms[1]->getNorm();
    if(_norms.size() > 5) {
      _norms[5]->computeL2Norm(_metaNumber, _sol, _mesh);
      _normL2[5][_currentStep] = _norms[5]->getNorm();
      _norms[6]->computeL2Norm(_metaNumber, _sol, _mesh);
      _normL2[6][_currentStep] = _norms[6]->getNorm();
    }
    // at this moment we fixe the DC3F à the first time step

    if(i == 1) {
      initializeDC3F_centered(
        _sol, _metaNumber, _mesh, dynamic_cast<feSolutionDC2F *>(_solutionContainerDC2F),
        dynamic_cast<feSolutionDC2F *>(
          _solutionContainer)); // soit on creer un fonction soit on passe un paramêtre a
                                // la fonction car seul tn change en réalité
      printf("\n");
      printf("Retour a l'iteration precedante pour calculer t1 - recomputeMatrix = %s : Solution "
             "DC3F - t = %6.6e\n",
             _linearSystem->getRecomputeStatus() ? "true" : "false", _sol->getCurrentTime());
      solveQNBDF(_solutionContainer, _tol, _metaNumber, _linearSystem, _sol, _mesh);
      fePstClc(_sol, _linearSystem, _solutionContainer);
      // Compute L2 norm of DC3F solution
      std::cout << "la norme vaut " << _norms[2]->getNorm() << std::endl;
      _norms[2]->computeL2Norm(_metaNumber, _sol, _mesh);
      _normL2[2][_currentStep - 1] = _norms[2]->getNorm();
      if(_norms.size() > 7) {
        _norms[7]->computeL2Norm(_metaNumber, _sol, _mesh);
        _normL2[7][_currentStep] = _norms[7]->getNorm();
        _norms[8]->computeL2Norm(_metaNumber, _sol, _mesh);
        _normL2[8][_currentStep] = _norms[8]->getNorm();
      }
    }
    if(i > 1 && i < nSteps) {
      _solutionContainer->rotate(_dt);
      initializeDC3F_centered(_sol, _metaNumber, _mesh,
                              dynamic_cast<feSolutionDC2F *>(_solutionContainerDC2F),
                              dynamic_cast<feSolutionDC2F *>(_solutionContainer));
      printf("\n");
      printf("Étape 3 - recomputeMatrix = %s : Solution DC3F - t = %6.6e\n",
             _linearSystem->getRecomputeStatus() ? "true" : "false", _sol->getCurrentTime());
      solveQNBDF(_solutionContainer, _tol, _metaNumber, _linearSystem, _sol, _mesh);
      fePstClc(_sol, _linearSystem, _solutionContainer);
      // Compute L2 norm of DC3F solution
      _norms[2]->computeL2Norm(_metaNumber, _sol, _mesh);
      _normL2[2][_currentStep] = _norms[2]->getNorm();
      if(_norms.size() > 7) {
        _norms[7]->computeL2Norm(_metaNumber, _sol, _mesh);
        _normL2[7][_currentStep] = _norms[7]->getNorm();
        _norms[8]->computeL2Norm(_metaNumber, _sol, _mesh);
        _normL2[8][_currentStep] = _norms[8]->getNorm();
      }
    }
    // Compute L2 norm of the solution
    _sol->setSolFromContainer(_solutionContainerBDF1);
    _sol->setSolFromContainer(_solutionContainerDC2F);
    _sol->setSolFromContainer(_solutionContainer);
    _tCurrent += _dt;
    ++_currentStep;
    if(K1K2) _dt = tK1K2[i + 1] - tK1K2[i];
    std::cout << "========dt vaux ==========" << _dt << std::endl;
    printf("\n");
    printf("Current step = %d : t = %f\n", _currentStep, _tCurrent);

    // std::string vtkFile = "../../data/cylindreAdapt" + std::to_string(_currentStep) + ".vtk";
    // feExporterVTK writer(vtkFile, _mesh, _sol, _metaNumber, spaces);
  }
  // std::cout<<"Solution finale du DC2"<<std::endl;
  // _sol->printSol();

  return FE_STATUS_OK;
}

DC3Solver::DC3Solver(feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem,
                     feSolution *sol, std::vector<feNorm *> &norms, feMesh *mesh,
                     feExportData exportData, double t0, double tEnd, int nTimeSteps,
                     std::string CodeIni)
  : TimeIntegrator(tol, metaNumber, linearSystem, sol, norms, mesh, exportData, t0, tEnd,
                   nTimeSteps, CodeIni)
{
  // Initialize the solution container
  int nSol = 5;
  //=====================================================================//
  // Define the refinement of the initialization
  int Nb_pas_de_temps = 2;
  //=====================================================================//
  _solutionContainer = new feSolutionDCF(nSol, _sol->getCurrentTime(), _metaNumber);
  _solutionContainerBDF2 = new feSolutionBDF2(nSol, _sol->getCurrentTime(), _metaNumber);
  _solutionContainerDC3F =
    new feSolutionDCF(Nb_pas_de_temps + 1, _sol->getCurrentTime(), _metaNumber);
  _solutionContainerBDF2->initialize(_sol, _mesh, _metaNumber);
  _solutionContainerDC3F->initialize(_sol, _mesh, _metaNumber);
  _solutionContainer->initialize(_sol, _mesh, _metaNumber);
  if(_CodeIni != "BDF1/DCF") {
    fePstClc(_sol, _linearSystem, _solutionContainerBDF2); // To initialize residual if no BDF1/DCF
  }
  fePstClc(_sol, _linearSystem, _solutionContainer);

  _normL2.resize(norms.size());
  for(auto &n : _normL2) n.resize(_nTimeSteps, 0.);

  printf("Initializing DC3 solver : integrating from t0 = %f to tEnd = %f in %d steps\n", _t0,
         _tEnd, _nTimeSteps);
}

feStatus DC3Solver::makeSteps(int nSteps)
{
  printf("BDF2 : Advancing %d steps from t = %f to t = %f\n", nSteps, _tCurrent,
         _tCurrent + nSteps * _dt);
  if(_currentStep == 0) {
    if(_CodeIni == "BDF1/DCF") {
      printf("\n");
      printf(" ----------------------------- \n");
      printf("Using DC3F to initialize DC3 \n");
      printf(" ----------------------------- \n");
      printf("\n");
      // _solutionContainer->rotate(_dt);
      // _solutionContainerBDF2->rotate(_dt);
      double true_dt = _dt;
      int Nb_pas_de_temps = 2;
      _dt = _dt / (Nb_pas_de_temps - 1);
      double _t_ini = _t0 + (Nb_pas_de_temps)*_dt;
      std::cout << "t ini = " << _t_ini << " et _dt = " << _dt << std::endl;
      // feNorm *normBDF1 = new feNorm(&U_M1D, mesh, nQuad, funSol);
      // feNorm *normDC2F = new feNorm(&U_M1D, mesh, nQuad, funSol);
      // feNorm *normDC3F = new feNorm(&U_M1D, mesh, nQuad, funSol);
      std::vector<feNorm *> norms2 = {_norms[1], _norms[1],
                                      _norms[1]}; // Il faut un vecteur norme de taille 3 sinon pb
      DC3FSolver solver(_tol, _metaNumber, _linearSystem, _sol, norms2, _mesh, _exportData, _t0,
                        _t_ini, Nb_pas_de_temps);
      solver.makeSteps(Nb_pas_de_temps);
      _solutionContainerDC3F = solver.getSolutionContainer();
      std::cout << "taille du container = " << _solutionContainerDC3F->getNbSol() << std::endl;
      _solutionContainerBDF2->setSol(0, _solutionContainerDC3F->getSolution(1));
      _solutionContainerBDF2->setSol(1, _solutionContainerDC3F->getSolution(Nb_pas_de_temps));
      _solutionContainerBDF2->setSol(2, _solutionContainerDC3F->getSolution(Nb_pas_de_temps));
      _solutionContainerBDF2->setResidual(0, _solutionContainerDC3F->getResidual(1));
      _solutionContainerBDF2->setResidual(1, _solutionContainerDC3F->getResidual(Nb_pas_de_temps));
      _solutionContainerBDF2->setResidual(2, _solutionContainerDC3F->getResidual(Nb_pas_de_temps));
      // std::cout<<"residu vaut "<<  _solutionContainerDC3F->getRes(1, 0) << std::endl;
      _solutionContainer->setSol(0, _solutionContainerDC3F->getSolution(1));
      _solutionContainer->setSol(1, _solutionContainerDC3F->getSolution(Nb_pas_de_temps));
      _solutionContainer->setSol(2, _solutionContainerDC3F->getSolution(Nb_pas_de_temps));
      // _norms[0]->computeL2Norm(_metaNumber, _sol, _mesh);
      // _norms[0] = norms2[1];
      // _norms[1] = norms2[2];
      // _normL2[0][_currentStep] = _norms[0]->getNorm();
      // _norms[1]->computeL2Norm(_metaNumber, _sol, _mesh);
      // _normL2[1][_currentStep] =_norms[1]->getNorm();
      // std::vector<double> &testBDF1 = solver.getNorm(0);
      std::vector<double> &DC2F = solver.getNorm(1);
      std::vector<double> &DC3F = solver.getNorm(2);
      _normL2[0][_currentStep] = DC2F[0]; // value at t1
      _normL2[1][_currentStep] = DC3F[0]; // value ay t1
      // for(int j=0; j<testBDF1.size() ; j++)std::cout<<"Norme BDF1  "<<testBDF1[j]<<std::endl;
      for(int j = 0; j < DC2F.size(); j++) std::cout << "Norme DC2F  " << DC2F[j] << std::endl;
      for(int j = 0; j < DC3F.size(); j++) std::cout << "DC3F   " << DC3F[j] << std::endl;
      // for(int j=0; j<2 ; j++)std::cout<<"Norme DC2F  "<<_norms[0]->getNorm()<<std::endl;
      // for(int j=0; j<2 ; j++)std::cout<<"Norme DC3F  "<<_norms[1]->getNorm()<<std::endl;

      _dt = true_dt;
      _solutionContainer->rotate(_dt);
      _solutionContainerBDF2->rotate(_dt);
      _sol->setSolFromContainer(_solutionContainer);

    } else {
      _solutionContainer->rotate(_dt);
      // _solutionContainer->rotate(_dt);
      _solutionContainerBDF2->rotate(_dt);
      _sol->setSolFromContainer(_solutionContainer);
      _solutionContainerBDF2->initialize(_sol, _mesh, _metaNumber);
      _solutionContainer->initialize(_sol, _mesh, _metaNumber);
      fePstClc(_sol, _linearSystem, _solutionContainerBDF2);
    }

    _linearSystem->setRecomputeStatus(true);
    _sol->setSolFromContainer(_solutionContainerBDF2);

    _tCurrent += _dt;
    ++_currentStep;
    printf("Current step = %d/%d : t = %f\n", _currentStep, nSteps, _tCurrent);

    --nSteps; // To advance the same number of steps than if currentStep != 0
  }

  std::vector<double> tK1K2(nSteps + 1, 0.0);
  for(int i = 0; i < nSteps + 1; ++i) tK1K2[i] = _t0 + i * _dt;
  for(int i = 0; i < nSteps; i += 2)
    if(i + 2 < nSteps + 1) tK1K2[i + 1] = tK1K2[i] + _f * (tK1K2[i + 2] - tK1K2[i]);

  // Integration with BDF2 and DC3F
  for(int i = 0; i < nSteps; ++i) {
    if(i == 0) {
      _linearSystem->setRecomputeStatus(true);
      if(_CodeIni == "BDF1/DCF") {
        _solutionContainerBDF2->setSol(1, _solutionContainerBDF2->getSolution(2));
        _solutionContainerBDF2->setResidual(1, _solutionContainerBDF2->getResidual(2));
        _solutionContainer->setSol(1, _solutionContainer->getSolution(2));
      }
    } else
      _linearSystem->setRecomputeStatus(false);

    _solutionContainerBDF2->rotate(_dt);
    _solutionContainer->rotate(_dt);
    initializeBDF2(_sol, _metaNumber, _mesh,
                   dynamic_cast<feSolutionBDF2 *>(_solutionContainerBDF2));

    printf("\n");
    printf("Étape 1 - recomputeMatrix = %s : Solution BDF2 - t = %6.6e\n",
           _linearSystem->getRecomputeStatus() ? "true" : "false", _sol->getCurrentTime());
    solveQNBDF(_solutionContainerBDF2, _tol, _metaNumber, _linearSystem, _sol, _mesh);
    fePstClc(_sol, _linearSystem, _solutionContainerBDF2);
    _sol->setSolFromContainer(_solutionContainerBDF2);
    // Compute L2 norm of BDF1 solution
    for(int k = 0; k < _norms.size() - 1; k += 2) {
      if(_norms[k]->getNbFields() > 1) {
        if(_norms[k]->getTypeNorm() == "NormH1") {
          _norms[k]->computeH1NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "NormSemiH1") {
          _norms[k]->computeSemiH1NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "NormL2") {
          _norms[k]->computeL2NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else {
          printf("%s\n", "Norm doesn't exist");
        }

      } else {
        if(_norms[k]->getTypeNorm() == "IntFlux") {
          _norms[k]->computeNormLambda(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "IntForce") {
          _norms[k]->computeNormLambdaNS(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "IntLambda") {
          _norms[k]->computeIntegralNum(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else {
          _norms[k]->computeL2Norm(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        }
      }
    }
    initializeDC3(_sol, _metaNumber, _mesh, dynamic_cast<feSolutionBDF2 *>(_solutionContainerBDF2),
                  dynamic_cast<feSolutionDCF *>(_solutionContainer));
    printf("\n");
    printf("Étape 2 - recomputeMatrix = %s : Solution DC3 - t = %6.6e\n",
           _linearSystem->getRecomputeStatus() ? "true" : "false", _sol->getCurrentTime());
    solveQNBDF(_solutionContainer, _tol, _metaNumber, _linearSystem, _sol, _mesh);
    fePstClc(_sol, _linearSystem, _solutionContainer);

    for(int k = 1; k < _norms.size(); k += 2) {
      if(_norms[k]->getNbFields() > 1) {
        if(_norms[k]->getTypeNorm() == "NormH1") {
          _norms[k]->computeH1NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "NormSemiH1") {
          _norms[k]->computeSemiH1NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "NormL2") {
          _norms[k]->computeL2NormVec(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else {
          printf("%s\n", "Norm doesn't exist");
        }

      } else {
        if(_norms[k]->getTypeNorm() == "IntFlux") {
          _norms[k]->computeNormLambda(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else if(_norms[k]->getTypeNorm() == "IntForce") {
          _norms[k]->computeNormLambdaNS(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
          std::cout << "Norm L-DC3" << _normL2[k][_currentStep] << std::endl;
        } else if(_norms[k]->getTypeNorm() == "IntLambda") {
          _norms[k]->computeIntegralNum(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        } else {
          _norms[k]->computeL2Norm(_metaNumber, _sol, _mesh);
          _normL2[k][_currentStep] = _norms[k]->getNorm();
        }
      }
    }

    _sol->setSolFromContainer(_solutionContainer);

    _tCurrent += _dt;
    ++_currentStep;
    if(K1K2) _dt = tK1K2[i + 1] - tK1K2[i];
    printf("\n");
    printf("Current step = %d : t = %f\n", _currentStep, _tCurrent);

    // std::string vtkFile = "../../data/TestNS_DC3/NS_Square" + std::to_string(_currentStep) +
    // ".vtk"; feExporterVTK writer(vtkFile, _mesh, _sol, _metaNumber, spaces);
  }

  return FE_STATUS_OK;
}
