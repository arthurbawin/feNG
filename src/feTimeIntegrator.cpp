#include "feTimeIntegrator.h"
#include "feExporter.h"
#include "feNG.h"


// STATIONARY, BDF1, BDF2, DC2F, DC3F, DC4F, DC5F, DC3, DC4, DC5, DC6
feStatus createTimeIntegratorV2(feTimeIntegrator *&timeIntegrator, feMesh *mesh, feSolverV2 *solver, feLinearSystem *linearSystem, feMetaNumber *metaNumber, feSolution *currentSolution, timeIntegratorSchemeV2 scheme, std::string correcType, double t0, double tEnd, int nTimeStep, std::vector<feNormV2 *> &norms, feExportData exporter, std::string codeDt)
{
  switch(scheme) {
    case STATIONARY_V2:
      timeIntegrator = new feTimeIntegratorStationary(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, t0, norms, exporter);
      break;
    
    case BDF1_V2:
      timeIntegrator = new feTimeIntegratorBDF1(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, t0, tEnd, nTimeStep, norms, exporter, codeDt);
      break;
    
    case BDF2_V2:
      timeIntegrator = new feTimeIntegratorBDF2(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, t0, tEnd, nTimeStep, norms, exporter, codeDt);
      break;
    
    case DC2F_V2:
      timeIntegrator = new feTimeIntegratorDC2F(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, correcType, t0, tEnd, nTimeStep, norms, exporter, codeDt);  
      break;
    
  	case DC3F_V2:
		  timeIntegrator = new feTimeIntegratorDC3F(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, correcType, t0, tEnd, nTimeStep, norms, exporter, codeDt);
		  break;

		case DC4F_V2:
	      timeIntegrator = new feTimeIntegratorDC4F(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, correcType, t0, tEnd, nTimeStep, norms, exporter, codeDt);
	      break;

		case DC5F_V2:
	      timeIntegrator = new feTimeIntegratorDC5F(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, correcType,  t0, tEnd, nTimeStep, norms, exporter, codeDt);
	      break;

    case DC3_V2:
      timeIntegrator = new feTimeIntegratorDC3(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, correcType, t0, tEnd, nTimeStep, norms, exporter, codeDt);
      break;
    
    case DC4_V2:
      timeIntegrator = new feTimeIntegratorDC4(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, correcType, t0, tEnd, nTimeStep, norms, exporter, codeDt);
      break;
    
    case DC5_V2:
      timeIntegrator = new feTimeIntegratorDC5(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, correcType, t0, tEnd, nTimeStep, norms, exporter, codeDt);
      break;
    
    case DC6_V2:
      timeIntegrator = new feTimeIntegratorDC6(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, correcType, t0, tEnd, nTimeStep, norms, exporter, codeDt);
      break;

    default:
	      return feErrorMsg(FE_STATUS_ERROR, "Unsupported time integration scheme.");
  }
  return FE_STATUS_OK;
}



feTimeIntegrator::feTimeIntegrator(feMesh *mesh, feSolverV2 *solver, feLinearSystem *linearSystem, feMetaNumber *metaNumber, feSolution *currentSolution, timeIntegratorSchemeV2 scheme, double t0, double tEnd, int nTimeStep, std::vector<feNormV2 *> &norms, feExportData exporter, std::string codeDt): 
		_mesh(mesh), _solver(solver), _linearSystem(linearSystem), _metaNumber(metaNumber), _currentSolution(currentSolution), _scheme(scheme), _norms(norms), _exportData(exporter), _t0(t0), _tEnd(tEnd), _nTimeStep(nTimeStep), _codeDt(codeDt)
{
	_nbDDL = _metaNumber->getNbDOFs();
	_nbInc = _metaNumber->getNbUnknowns();

	this->initializeTime();

}


void feTimeIntegrator::pstClc(feSolutionContainerV2 *solutionContainer)
{
  _currentSolution->setSolDotToZero();
  _linearSystem->setResidualToZero();
  _linearSystem->assembleResiduals(_currentSolution);
  _linearSystem->assignResidualToDCResidualV2(solutionContainer);
}





// ===============================================================
// TimeIntegrator Initialization 
// ===============================================================

void feTimeIntegrator::initialize()
{
  _t[0]=_tC;
  _currentSolution->setCurrentTime(_tC);
  _currentSolution->initializeUnknowns(_mesh, _metaNumber);
  _currentSolution->initializeEssentialBC(_mesh, _metaNumber);
}


void feTimeIntegrator::initializeUDOT()
{
	_solver->solveUDOT(_currentSolution);
	// _currentSolution->displaySolDot();
}



void feTimeIntegratorStationary::initializeStat()
{
	this->initialize();
	_linearSystem->setRecomputeStatus(true);
	
	_solutionContainerStat->copyCurrentSolution(_currentSolution);
	
	_solutionContainerStat->copyCurrentSolutionDot(_currentSolution);
}

void feTimeIntegratorBDF1::initializeBDF1()
{
	this->initialize();
	_linearSystem->setRecomputeStatus(true);
	
	_solutionContainerBDF1->copyCurrentSolution(_currentSolution);
	
	_solutionContainerBDF1->copyCurrentSolutionDot(_currentSolution);
}


void feTimeIntegratorBDF2::initializeBDF2()
{
	this->initialize();
	_linearSystem->setRecomputeStatus(true);
		
	_solutionContainerBDF2->copyCurrentSolution(_currentSolution);
	
	_solutionContainerBDF2->copyCurrentSolutionDot(_currentSolution);
}

void feTimeIntegratorDC2F::initializeDC2F()
{
	_linearSystem->setRecomputeStatus(true);

	this->initialize();

	if (_correctionType =="SOLDOT"){
		this->initializeUDOT();
	}

	_solutionContainerBDF1->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC2F->copyCurrentSolution(_currentSolution);
 	
 	_solutionContainerBDF1->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC2F->copyCurrentSolutionDot(_currentSolution);

 	if (_correctionType=="RESIDUAL"){
 		this->pstClc(_solutionContainerBDF1);
 	}
}

void feTimeIntegratorDC3F::initializeDC3F()
{
	_linearSystem->setRecomputeStatus(true);
	
	this->initialize();
	
	// if (_correctionType =="SOLDOT"){
	// 	this->initializeUDOT();
	// }
	
	_solutionContainerBDF1->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC2F->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC3F->copyCurrentSolution(_currentSolution);
 	
 	_solutionContainerBDF1->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC2F->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC3F->copyCurrentSolutionDot(_currentSolution);

 	if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF1);
			this->pstClc(_solutionContainerDC2F);
	}
}


void feTimeIntegratorDC4F::initializeDC4F()
{
	_linearSystem->setRecomputeStatus(true);
	
	this->initialize();
	
	if (_correctionType =="SOLDOT"){
		this->initializeUDOT();
	}
	
	_solutionContainerBDF1->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC2F->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC3F->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC4F->copyCurrentSolution(_currentSolution);
 	
 	_solutionContainerBDF1->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC2F->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC3F->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC4F->copyCurrentSolutionDot(_currentSolution);

 	if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF1);
			this->pstClc(_solutionContainerDC2F);
			this->pstClc(_solutionContainerDC3F);
	}
}


void feTimeIntegratorDC5F::initializeDC5F()
{
	_linearSystem->setRecomputeStatus(true);
	
	this->initialize();
	
	// if (_correctionType =="SOLDOT"){
	// 	this->initializeUDOT();
	// }
	
	_solutionContainerBDF1->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC2F->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC3F->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC4F->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC5F->copyCurrentSolution(_currentSolution);
 	
 	_solutionContainerBDF1->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC2F->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC3F->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC4F->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC5F->copyCurrentSolutionDot(_currentSolution);

 	if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF1);
			this->pstClc(_solutionContainerDC2F);
			this->pstClc(_solutionContainerDC3F);
			this->pstClc(_solutionContainerDC4F);
	}
}


void feTimeIntegratorDC3::initializeDC3()
{
	_linearSystem->setRecomputeStatus(true);
	
	this->initialize();
	
	// if (_correctionType =="SOLDOT"){
	// 	this->initializeUDOT();
	// }

	_solutionContainerBDF2->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC3->copyCurrentSolution(_currentSolution);

 	_solutionContainerBDF2->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC3->copyCurrentSolutionDot(_currentSolution);

 	if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF2);
	}
}


void feTimeIntegratorDC4::initializeDC4()
{
	_linearSystem->setRecomputeStatus(true);
	
	this->initialize();
	
	// if (_correctionType =="SOLDOT"){
	// 	this->initializeUDOT();
	// }

	_solutionContainerBDF2->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC3->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC4->copyCurrentSolution(_currentSolution);

 	_solutionContainerBDF2->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC3->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC4->copyCurrentSolutionDot(_currentSolution);

 	if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF2);
			this->pstClc(_solutionContainerDC3);
	}
}


void feTimeIntegratorDC5::initializeDC5()
{
	_linearSystem->setRecomputeStatus(true);
	
	this->initialize();
	
	// if (_correctionType =="SOLDOT"){
	// 	this->initializeUDOT();
	// }

	_solutionContainerBDF2->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC3->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC4->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC5->copyCurrentSolution(_currentSolution);

 	_solutionContainerBDF2->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC3->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC4->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC5->copyCurrentSolutionDot(_currentSolution);

 	if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF2);
			this->pstClc(_solutionContainerDC3);
			this->pstClc(_solutionContainerDC4);
	}
}


void feTimeIntegratorDC6::initializeDC6()
{
	_linearSystem->setRecomputeStatus(true);
	
	this->initialize();
	
	// if (_correctionType =="SOLDOT"){
	// 	this->initializeUDOT();
	// }

	_solutionContainerBDF2->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC3->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC4->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC5->copyCurrentSolution(_currentSolution);
 	_solutionContainerDC6->copyCurrentSolution(_currentSolution);

 	_solutionContainerBDF2->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC3->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC4->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC5->copyCurrentSolutionDot(_currentSolution);
 	_solutionContainerDC6->copyCurrentSolutionDot(_currentSolution);

 	if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF2);
			this->pstClc(_solutionContainerDC3);
			this->pstClc(_solutionContainerDC4);
			this->pstClc(_solutionContainerDC5);
	}
}

// ===============================================================
// TimeIntegrator MakeStep 
// ===============================================================

feStatus feTimeIntegratorStationary::makeSteps(std::string startCode, std::string codeIniDC)
{
	_solutionContainerStat->computeCoeffBDF(_k);
	_solutionContainerStat->computeCorrection();

	_solver->solve(_currentSolution, _solutionContainerStat);

	for (size_t i=0; i<_norms.size(); i++){
			_resultNormStat[i][0]=_norms[i]->compute(_tC, _solutionContainerStat);
		}
	
	if(_exportData.exporter != nullptr) {
        std::string fileName = _exportData.fileNameRoot + ".vtk";
        _exportData.exporter->writeStep(fileName);
  }
}


feStatus feTimeIntegratorBDF1::makeSteps(std::string startCode, std::string codeIniDC)
{
	for(int n=0; n<_nTimeStep; n++){
		printf("=========================================\n");
		this->updateTime(n);

		_solutionContainerBDF1->rotate();

		_solutionContainerBDF1->computeCoeffBDF(_k);

		_solutionContainerBDF1->initializeBDFSolution(_mesh, _metaNumber, _currentSolution);


		feInfo("Solving BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);
		_solutionContainerBDF1->computeCorrection();
		_solver->solve(_currentSolution, _solutionContainerBDF1);

		for (size_t i=0; i<_norms.size(); i++){
			_resultNormBDF1[i][n]=_norms[i]->compute(_tC, _solutionContainerBDF1);
			// feInfo("%10.12f",_resultNormBDF1[i][n]);
		}

		if(_exportData.exporter != nullptr && ((n + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_BDF1_" + std::to_string(n+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
  	}

	}

	//maj exporter
}


feStatus feTimeIntegratorBDF2::makeSteps(std::string startCode, std::string codeIniDC)
{
	this->startBDF2(startCode);
	for(int n=1; n<_nTimeStep; n++){
		printf("=========================================\n");
		this->updateTime(n);

		_solutionContainerBDF2->rotate();

		_solutionContainerBDF2->computeCoeffBDF(_k);

		_solutionContainerBDF2->initializeBDFSolution(_mesh, _metaNumber, _currentSolution);


		feInfo("Solving BDF2 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);
		_solutionContainerBDF2->computeCorrection();
		_solver->solve(_currentSolution, _solutionContainerBDF2);

		for (size_t i=0; i<_norms.size(); i++){
			_resultNormBDF2[i][n]=_norms[i]->compute(_tC, _solutionContainerBDF2);
		}

		if(_exportData.exporter != nullptr && ((n + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_BDF2_" + std::to_string(n+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
  	}

	}
}



feStatus feTimeIntegratorDC2F::makeSteps(std::string startCode, std::string codeIniDC)
{
	for(int n=0; n<_nTimeStep; n++){
		printf("=========================================\n");
		this->updateTime(n);

		_solutionContainerBDF1->rotate();
		_solutionContainerDC2F->rotate();

		_solutionContainerBDF1->computeCoeffBDF(_k);
		_solutionContainerDC2F->computeCoeffBDF(_k);
		
		feInfo("Solving BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerBDF1->initializeBDFSolution(_mesh, _metaNumber, _currentSolution);
		
		_solver->solve(_currentSolution, _solutionContainerBDF1);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF1);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC2/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);
		_solutionContainerDC2F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC2F->computeCorrection(_solutionContainerBDF1, _k);

		_solver->solve(_currentSolution, _solutionContainerDC2F);


		for (size_t i=0; i<_norms.size(); i++){
			_resultNormBDF1[i][n]=_norms[i]->compute(_tC, _solutionContainerBDF1, _solutionContainerDC2F);
			_resultNormDC2F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC2F, _solutionContainerDC2F);
		}

		if(_exportData.exporter != nullptr && ((n + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC2F_" + std::to_string(n+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
  	}

	}

	//maj exporter
}



feStatus feTimeIntegratorDC3F::makeSteps(std::string startCode, std::string codeIniDC)
{
	this->startDC3F(startCode, codeIniDC);
	for(int n=1; n<_nTimeStep; n++){
		printf("=========================================\n");
		this->updateTime(n);

		_solutionContainerBDF1->rotate();
		_solutionContainerDC2F->rotate();
		_solutionContainerDC3F->rotate();

		_solutionContainerBDF1->computeCoeffBDF(_k);
		_solutionContainerDC2F->computeCoeffBDF(_k);
		_solutionContainerDC3F->computeCoeffBDF(_k);

		// _solutionContainerDC3F->printSolution(0);
		// printf("\n");
		// printf("\n");
		// _solutionContainerDC3F->printSolution(1);
		// printf("\n");
		// printf("\n");
		// _solutionContainerDC3F->printSolution(2);

		feInfo("Solving BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerBDF1->initializeBDFSolution(_mesh, _metaNumber, _currentSolution);
		
		_solutionContainerBDF1->computeCorrection();

		_solver->solve(_currentSolution, _solutionContainerBDF1);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF1);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC2/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerDC2F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC2F->computeCorrection(_solutionContainerBDF1, _k);
		
		_solver->solve(_currentSolution, _solutionContainerDC2F);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerDC2F);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC3/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerDC3F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC3F->computeCorrection(_solutionContainerDC2F, _k);
		
		_solver->solve(_currentSolution, _solutionContainerDC3F);


		for (size_t i=0; i<_norms.size(); i++){
			_resultNormBDF1[i][n]=_norms[i]->compute(_tC, _solutionContainerBDF1, _solutionContainerDC3F);
			_resultNormDC2F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC2F, _solutionContainerDC3F);
			_resultNormDC3F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC3F, _solutionContainerDC3F);
		}

		if(_exportData.exporter != nullptr && ((n + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC3F_" + std::to_string(n+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
  	}

	}

	//maj exporter
}




feStatus feTimeIntegratorDC4F::makeSteps(std::string startCode, std::string codeIniDC)
{
	this->startDC4F(startCode, codeIniDC);
	for(int n=2; n<_nTimeStep; n++){
		printf("=========================================\n");
		this->updateTime(n);

		_solutionContainerBDF1->rotate();
		_solutionContainerDC2F->rotate();
		_solutionContainerDC3F->rotate();
		_solutionContainerDC4F->rotate();

		_solutionContainerBDF1->computeCoeffBDF(_k);
		_solutionContainerDC2F->computeCoeffBDF(_k);
		_solutionContainerDC3F->computeCoeffBDF(_k);
		_solutionContainerDC4F->computeCoeffBDF(_k);

		feInfo("Solving BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerBDF1->initializeBDFSolution(_mesh, _metaNumber, _currentSolution);
		
		_solutionContainerBDF1->computeCorrection();

		_solver->solve(_currentSolution, _solutionContainerBDF1);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF1);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC2/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerDC2F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC2F->computeCorrection(_solutionContainerBDF1, _k);

		_solver->solve(_currentSolution, _solutionContainerDC2F);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerDC2F);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC3/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerDC3F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC3F->computeCorrection(_solutionContainerDC2F, _k);
		
		_solver->solve(_currentSolution, _solutionContainerDC3F);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerDC3F);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC4/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerDC4F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC4F->computeCorrection(_solutionContainerDC3F, _k);
		
		_solver->solve(_currentSolution, _solutionContainerDC4F);


		for (size_t i=0; i<_norms.size(); i++){
			_resultNormBDF1[i][n]=_norms[i]->compute(_tC, _solutionContainerBDF1, _solutionContainerDC4F);
			_resultNormDC2F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC2F, _solutionContainerDC4F);
			_resultNormDC3F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC3F, _solutionContainerDC4F);
			_resultNormDC4F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC4F, _solutionContainerDC4F);
		}

		if(_exportData.exporter != nullptr && ((n + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC4F_" + std::to_string(n+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  }

	}

}


feStatus feTimeIntegratorDC5F::makeSteps(std::string startCode, std::string codeIniDC)
{
	this->startDC5F(startCode, codeIniDC);
	for(int n=3; n<_nTimeStep; n++){
		printf("=========================================\n");
		this->updateTime(n);

		_solutionContainerBDF1->rotate();
		_solutionContainerDC2F->rotate();
		_solutionContainerDC3F->rotate();
		_solutionContainerDC4F->rotate();
		_solutionContainerDC5F->rotate();

		_solutionContainerBDF1->computeCoeffBDF(_k);
		_solutionContainerDC2F->computeCoeffBDF(_k);
		_solutionContainerDC3F->computeCoeffBDF(_k);
		_solutionContainerDC4F->computeCoeffBDF(_k);
		_solutionContainerDC5F->computeCoeffBDF(_k);
		

		feInfo("Solving BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerBDF1->initializeBDFSolution(_mesh, _metaNumber, _currentSolution);
		
		_solutionContainerBDF1->computeCorrection();

		_solver->solve(_currentSolution, _solutionContainerBDF1);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF1);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC2/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerDC2F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC2F->computeCorrection(_solutionContainerBDF1, _k);
		
		_solver->solve(_currentSolution, _solutionContainerDC2F);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerDC2F);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC3/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerDC3F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC3F->computeCorrection(_solutionContainerDC2F, _k);
		
		_solver->solve(_currentSolution, _solutionContainerDC3F);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerDC3F);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC4/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerDC4F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC4F->computeCorrection(_solutionContainerDC3F, _k);
		
		_solver->solve(_currentSolution, _solutionContainerDC4F);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerDC4F);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC5/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerDC5F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC5F->computeCorrection(_solutionContainerDC4F, _k);
		
		_solver->solve(_currentSolution, _solutionContainerDC5F);


		for (size_t i=0; i<_norms.size(); i++){
			_resultNormBDF1[i][n]=_norms[i]->compute(_tC, _solutionContainerBDF1, _solutionContainerDC5F);
			_resultNormDC2F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC2F, _solutionContainerDC5F);
			_resultNormDC3F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC3F, _solutionContainerDC5F);
			_resultNormDC4F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC4F, _solutionContainerDC5F);
			_resultNormDC5F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC5F, _solutionContainerDC5F);
		}

		if(_exportData.exporter != nullptr && ((n + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC5F_" + std::to_string(n+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  }

	}
}


feStatus feTimeIntegratorDC3::makeSteps(std::string startCode, std::string codeIniDC)
{
	this->startDC3(startCode);
	for(int n=1; n<_nTimeStep; n++){
		printf("=========================================\n");
		this->updateTime(n);

		_solutionContainerBDF2->rotate();
		_solutionContainerDC3->rotate();

		_solutionContainerBDF2->computeCoeffBDF(_k);
		_solutionContainerDC3->computeCoeffBDF(_k);

		_solutionContainerDC3->computeCorrectionCoeff(_k);
		

		feInfo("Solving BDF2 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerBDF2->initializeBDFSolution(_mesh, _metaNumber, _currentSolution);
		
		_solutionContainerBDF2->computeCorrection();

		_solver->solve(_currentSolution, _solutionContainerBDF2);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF2);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC3/BDF2 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerDC3->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC3->computeCorrection(_solutionContainerBDF2, _k);
		
		_solver->solve(_currentSolution, _solutionContainerDC3);


		for (size_t i=0; i<_norms.size(); i++){
			_resultNormBDF2[i][n]= _norms[i]->compute(_tC, _solutionContainerBDF2, _solutionContainerDC3);
			_resultNormDC3[i][n] = _norms[i]->compute(_tC, _solutionContainerDC3, _solutionContainerDC3);
		}

		if(_exportData.exporter != nullptr && ((n + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC3_" + std::to_string(n+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  }

	}

	//maj exporter
}



feStatus feTimeIntegratorDC4::makeSteps(std::string startCode, std::string codeIniDC)
{
	this->startDC4(startCode);
	for(int n=2; n<_nTimeStep; n++){
		printf("=========================================\n");
		this->updateTime(n);

		_solutionContainerBDF2->rotate();
		_solutionContainerDC3->rotate();
		_solutionContainerDC4->rotate();

		_solutionContainerBDF2->computeCoeffBDF(_k);
		_solutionContainerDC3->computeCoeffBDF(_k);
		_solutionContainerDC4->computeCoeffBDF(_k);

		_solutionContainerDC3->computeCorrectionCoeff(_k);
		_solutionContainerDC4->computeCorrectionCoeff(_k);
		

		feInfo("Solving BDF2 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerBDF2->initializeBDFSolution(_mesh, _metaNumber, _currentSolution);
		
		_solutionContainerBDF2->computeCorrection();

		_solver->solve(_currentSolution, _solutionContainerBDF2);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF2);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC3/BDF2 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerDC3->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC3->computeCorrection(_solutionContainerBDF2, _k);
		
		_solver->solve(_currentSolution, _solutionContainerDC3);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerDC3);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC4/BDF2 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerDC4->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC4->computeCorrection(_solutionContainerDC3, _k);
		
		_solver->solve(_currentSolution, _solutionContainerDC4);


		for (size_t i=0; i<_norms.size(); i++){
			_resultNormBDF2[i][n]= _norms[i]->compute(_tC, _solutionContainerBDF2, _solutionContainerDC4);
			_resultNormDC3[i][n] = _norms[i]->compute(_tC, _solutionContainerDC3, _solutionContainerDC4);
			_resultNormDC4[i][n] = _norms[i]->compute(_tC, _solutionContainerDC4, _solutionContainerDC4);
		}

		if(_exportData.exporter != nullptr && ((n + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC4_" + std::to_string(n+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  }

	}
}



feStatus feTimeIntegratorDC5::makeSteps(std::string startCode, std::string codeIniDC)
{
	this->startDC5(startCode);
	for(int n=3; n<_nTimeStep; n++){
		printf("=========================================\n");
		this->updateTime(n);

		_solutionContainerBDF2->rotate();
		_solutionContainerDC3->rotate();
		_solutionContainerDC4->rotate();
		_solutionContainerDC5->rotate();

		_solutionContainerBDF2->computeCoeffBDF(_k);
		_solutionContainerDC3->computeCoeffBDF(_k);
		_solutionContainerDC4->computeCoeffBDF(_k);
		_solutionContainerDC5->computeCoeffBDF(_k);

		_solutionContainerDC3->computeCorrectionCoeff(_k);
		_solutionContainerDC4->computeCorrectionCoeff(_k);
		_solutionContainerDC5->computeCorrectionCoeff(_k);
		

		feInfo("Solving BDF2 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerBDF2->initializeBDFSolution(_mesh, _metaNumber, _currentSolution);
		
		_solutionContainerBDF2->computeCorrection();

		_solver->solve(_currentSolution, _solutionContainerBDF2);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF2);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC3/BDF2 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerDC3->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC3->computeCorrection(_solutionContainerBDF2, _k);
		
		_solver->solve(_currentSolution, _solutionContainerDC3);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerDC3);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC4/BDF2 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerDC4->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC4->computeCorrection(_solutionContainerDC3, _k);
		
		_solver->solve(_currentSolution, _solutionContainerDC4);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerDC4);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC5/BDF2 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerDC5->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC5->computeCorrection(_solutionContainerDC4, _k);
		
		_solver->solve(_currentSolution, _solutionContainerDC5);

		for (size_t i=0; i<_norms.size(); i++){
			_resultNormBDF2[i][n]= _norms[i]->compute(_tC, _solutionContainerBDF2, _solutionContainerDC5);
			_resultNormDC3[i][n] = _norms[i]->compute(_tC, _solutionContainerDC3, _solutionContainerDC5);
			_resultNormDC4[i][n] = _norms[i]->compute(_tC, _solutionContainerDC4, _solutionContainerDC5);
			_resultNormDC5[i][n] = _norms[i]->compute(_tC, _solutionContainerDC5, _solutionContainerDC5);
		}

		if(_exportData.exporter != nullptr && ((n + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC5_" + std::to_string(n+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  }

	}
}


feStatus feTimeIntegratorDC6::makeSteps(std::string startCode, std::string codeIniDC)
{
	this->startDC6(startCode);
	for(int n=4; n<_nTimeStep; n++){
		printf("=========================================\n");
		this->updateTime(n);

		_solutionContainerBDF2->rotate();
		_solutionContainerDC3->rotate();
		_solutionContainerDC4->rotate();
		_solutionContainerDC5->rotate();
		_solutionContainerDC6->rotate();

		_solutionContainerBDF2->computeCoeffBDF(_k);
		_solutionContainerDC3->computeCoeffBDF(_k);
		_solutionContainerDC4->computeCoeffBDF(_k);
		_solutionContainerDC5->computeCoeffBDF(_k);
		_solutionContainerDC6->computeCoeffBDF(_k);

		_solutionContainerDC3->computeCorrectionCoeff(_k);
		_solutionContainerDC4->computeCorrectionCoeff(_k);
		_solutionContainerDC5->computeCorrectionCoeff(_k);
		_solutionContainerDC6->computeCorrectionCoeff(_k);
		

		feInfo("Solving BDF2 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerBDF2->initializeBDFSolution(_mesh, _metaNumber, _currentSolution);
		
		_solutionContainerBDF2->computeCorrection();

		_solver->solve(_currentSolution, _solutionContainerBDF2);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF2);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC3/BDF2 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerDC3->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC3->computeCorrection(_solutionContainerBDF2, _k);
		
		_solver->solve(_currentSolution, _solutionContainerDC3);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerDC3);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC4/BDF2 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerDC4->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC4->computeCorrection(_solutionContainerDC3, _k);
		
		_solver->solve(_currentSolution, _solutionContainerDC4);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerDC4);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC5/BDF2 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerDC5->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC5->computeCorrection(_solutionContainerDC4, _k);
		
		_solver->solve(_currentSolution, _solutionContainerDC5);

		if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerDC5);
		}

		printf("----------------------------------------\n");
		feInfo("Solving DC6/BDF2 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

		_solutionContainerDC6->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC6->computeCorrection(_solutionContainerDC5, _k);
		
		_solver->solve(_currentSolution, _solutionContainerDC6);

		for (size_t i=0; i<_norms.size(); i++){
			_resultNormBDF2[i][n]= _norms[i]->compute(_tC, _solutionContainerBDF2, _solutionContainerDC6);
			_resultNormDC3[i][n] = _norms[i]->compute(_tC, _solutionContainerDC3, _solutionContainerDC6);
			_resultNormDC4[i][n] = _norms[i]->compute(_tC, _solutionContainerDC4, _solutionContainerDC6);
			_resultNormDC5[i][n] = _norms[i]->compute(_tC, _solutionContainerDC5, _solutionContainerDC6);
			_resultNormDC6[i][n] = _norms[i]->compute(_tC, _solutionContainerDC6, _solutionContainerDC6);
		}

		if(_exportData.exporter != nullptr && ((n + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC6_" + std::to_string(n+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  }

	}
}


//=======================================================================================
// Start
//=======================================================================================
void feTimeIntegratorBDF2::startBDF2(std::string startCode)
{
	printf("=========================================\n");
	printf("             Starting BDF2               \n");
	printf("=========================================\n");
	
	if (startCode=="Analytique"){
		this->updateTime(0);

		_solutionContainerBDF2->rotate();

		this->initializeBDF2();

		if(_exportData.exporter != nullptr && (1 % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_BDF2_" + std::to_string(1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  }
	
	} else if (startCode=="Consistant"){
		feInfo("TODO : demarrage consistant DC3F");
	
	} else {
		feErrorMsg(FE_STATUS_ERROR, "Unknown startCode for DC3F");
	}

	printf("=========================================\n");
	printf("          End Starting DC3F              \n");
	printf("=========================================\n");
	
}




void feTimeIntegratorDC3F::startDC3F(std::string startCode, std::string codeIniDC)
{
	printf("=========================================\n");
	printf("             Starting DC3F               \n");
	printf("=========================================\n");
	
	if (startCode=="Analytique"){
		this->updateTime(0);

		_solutionContainerBDF1->rotate();
		_solutionContainerDC2F->rotate();
		_solutionContainerDC3F->rotate();

		this->initializeDC3F();

		if(_exportData.exporter != nullptr && (1 % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC3F_" + std::to_string(1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  }
	
	} 

	else if (startCode=="Consistant"){
		feInfo("coucou");

		for (int n=0; n<2; n++){
			printf("=========================================\n");
			this->updateTime(n);

			_solutionContainerBDF1->rotate();
			_solutionContainerDC2F->rotate();

			_solutionContainerBDF1->computeCoeffBDF(_k);
			_solutionContainerDC2F->computeCoeffBDF(_k);

			feInfo("Solving BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);
			_solutionContainerBDF1->initializeBDFSolution(_mesh, _metaNumber, _currentSolution);

			_solutionContainerBDF1->computeCorrection();

			_solver->solve(_currentSolution, _solutionContainerBDF1);

			if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF1);
			}

			printf("----------------------------------------\n");
			feInfo("Solving DC2/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

			_solutionContainerDC2F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

			_solutionContainerDC2F->computeCorrection(_solutionContainerBDF1, _k);
			
			_solver->solve(_currentSolution, _solutionContainerDC2F);

			if (_correctionType=="RESIDUAL"){
				this->pstClc(_solutionContainerDC2F);
			}

			for (size_t i=0; i<_norms.size(); i++){
				_resultNormBDF1[i][n]=_norms[i]->compute(_tC, _solutionContainerBDF1, _solutionContainerDC3F);
				_resultNormDC2F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC2F, _solutionContainerDC3F);
			}

		}

		_solutionContainerDC3F->computeStartCorrection(_solutionContainerDC2F, _k);
		this->rebootTime();

		this->updateTime(0);

		_solutionContainerDC3F->rotate();
		_solutionContainerDC3F->computeCoeffBDF(_k);

		printf("----------------------------------------\n");
		feInfo("Solving DC3/BDF1 - Step = %d/%d  (%f s)", 0+1, _nTimeStep, _tC);

		_solutionContainerDC3F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);
		
		_solver->solve(_currentSolution, _solutionContainerDC3F);


		for (size_t i=0; i<_norms.size(); i++){
				_resultNormDC3F[i][0]=_norms[i]->compute(_tC, _solutionContainerDC3F, _solutionContainerDC3F);
				// feInfo("normDC3 : %10.10f",_resultNormDC3F[i][0]);
			}

		if(_exportData.exporter != nullptr && ((0 + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC3F_" + std::to_string(0+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
  	}


  	_solutionContainerBDF1->invRotate();
		_solutionContainerDC2F->invRotate();

	
	} else {
		feErrorMsg(FE_STATUS_ERROR, "Unknown startCode for DC3F");
	}

	printf("=========================================\n");
	printf("          End Starting DC3F              \n");
	printf("=========================================\n");
	
}



void feTimeIntegratorDC4F::startDC4F(std::string startCode, std::string codeIniDC)
{
	printf("=========================================\n");
	printf("             Starting DC4F               \n");
	printf("=========================================\n");
	
	if (startCode=="Analytique"){
		for(int n=0; n<2; n++){
			this->updateTime(n);

			_solutionContainerBDF1->rotate();
			_solutionContainerDC2F->rotate();
			_solutionContainerDC3F->rotate();
			_solutionContainerDC4F->rotate();

			this->initializeDC4F();

			if(_exportData.exporter != nullptr && ((n+1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC4F_" + std::to_string(n+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  	}
		}
	
	} 

	else if (startCode=="Consistant"){
		printf("=========================================\n");
		printf("       Consistant Starting DC3F          \n");
		printf("=========================================\n");
		for (int n=0; n<2; n++){
			printf("=========================================\n");
			this->updateTime(n);

			_solutionContainerBDF1->rotate();
			_solutionContainerDC2F->rotate();

			_solutionContainerBDF1->computeCoeffBDF(_k);
			_solutionContainerDC2F->computeCoeffBDF(_k);

			feInfo("Solving BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);
			_solutionContainerBDF1->initializeBDFSolution(_mesh, _metaNumber, _currentSolution);

			_solutionContainerBDF1->computeCorrection();

			_solver->solve(_currentSolution, _solutionContainerBDF1);

			if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF1);
			}

			printf("----------------------------------------\n");
			feInfo("Solving DC2/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

			_solutionContainerDC2F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

			_solutionContainerDC2F->computeCorrection(_solutionContainerBDF1, _k);
			
			_solver->solve(_currentSolution, _solutionContainerDC2F);

			if (_correctionType=="RESIDUAL"){
				this->pstClc(_solutionContainerDC2F);
			}

			for (size_t i=0; i<_norms.size(); i++){
				_resultNormBDF1[i][n]=_norms[i]->compute(_tC, _solutionContainerBDF1, _solutionContainerDC4F);
				_resultNormDC2F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC2F, _solutionContainerDC4F);
			}

		}

		_solutionContainerDC3F->computeStartCorrection(_solutionContainerDC2F, _k);
		this->rebootTime();

		this->updateTime(0);

		_solutionContainerDC3F->rotate();
		_solutionContainerDC3F->computeCoeffBDF(_k);

		printf("----------------------------------------\n");
		feInfo("Solving DC3/BDF1 - Step = %d/%d  (%f s)", 0+1, _nTimeStep, _tC);

		_solutionContainerDC3F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);
		
		_solver->solve(_currentSolution, _solutionContainerDC3F);


		for (size_t i=0; i<_norms.size(); i++){
				_resultNormDC3F[i][0]=_norms[i]->compute(_tC, _solutionContainerDC3F, _solutionContainerDC3F);
				// feInfo("normDC3 : %10.10f",_resultNormDC3F[i][0]);
			}

		if(_exportData.exporter != nullptr && ((0 + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC3F_" + std::to_string(0+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
  	}


  	_solutionContainerBDF1->invRotate();
		_solutionContainerDC2F->invRotate();	
		
		printf("=========================================\n");
		printf("     End Consistant Starting DC3F        \n");
		printf("=========================================\n");

		for (int n=1; n<3; n++){
			printf("=========================================\n");
			this->updateTime(n);

			_solutionContainerBDF1->rotate();
			_solutionContainerDC2F->rotate();
			_solutionContainerDC3F->rotate();

			_solutionContainerBDF1->computeCoeffBDF(_k);
			_solutionContainerDC2F->computeCoeffBDF(_k);
			_solutionContainerDC3F->computeCoeffBDF(_k);

			feInfo("Solving BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);
			_solutionContainerBDF1->initializeBDFSolution(_mesh, _metaNumber, _currentSolution);

			_solutionContainerBDF1->computeCorrection();

			_solver->solve(_currentSolution, _solutionContainerBDF1);

			if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF1);
			}

			printf("----------------------------------------\n");
			feInfo("Solving DC2/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

			_solutionContainerDC2F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

			_solutionContainerDC2F->computeCorrection(_solutionContainerBDF1, _k);
			
			_solver->solve(_currentSolution, _solutionContainerDC2F);

			if (_correctionType=="RESIDUAL"){
				this->pstClc(_solutionContainerDC2F);
			}

			printf("----------------------------------------\n");
			feInfo("Solving DC3/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

			_solutionContainerDC3F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

			_solutionContainerDC3F->computeCorrection(_solutionContainerDC2F, _k);
			
			_solver->solve(_currentSolution, _solutionContainerDC3F);

			if (_correctionType=="RESIDUAL"){
				this->pstClc(_solutionContainerDC3F);
			}


			for (size_t i=0; i<_norms.size(); i++){
				_resultNormBDF1[i][n]=_norms[i]->compute(_tC, _solutionContainerBDF1, _solutionContainerDC4F);
				_resultNormDC2F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC2F, _solutionContainerDC4F);
				_resultNormDC3F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC3F, _solutionContainerDC4F);
			}

		}

		_solutionContainerDC4F->computeStartCorrection(_solutionContainerDC3F, _k);
		this->rebootTime();

		printf("----------------------------------------\n");
		feInfo("Solving DC4/BDF1 - Step = %d/%d  (%f s)", 0+1, _nTimeStep, _tC);
		this->updateTime(0);

		_solutionContainerDC4F->rotate();

		_solutionContainerDC4F->computeCoeffBDF(_k);

		_solutionContainerDC4F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC4F->setStartCorrection(0);

		_solver->solve(_currentSolution, _solutionContainerDC4F);

		for (size_t i=0; i<_norms.size(); i++){
			_resultNormDC4F[i][0]=_norms[i]->compute(_tC, _solutionContainerDC4F, _solutionContainerDC4F);
		}

		if(_exportData.exporter != nullptr && ((0 + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC4F_" + std::to_string(0+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  }

	  printf("----------------------------------------\n");
		feInfo("Solving DC4/BDF1 - Step = %d/%d  (%f s)", 1+1, _nTimeStep, _tC);
		this->updateTime(1);

		_solutionContainerDC4F->rotate();

		_solutionContainerDC4F->computeCoeffBDF(_k);

		_solutionContainerDC4F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC4F->setStartCorrection(1);

		_solver->solve(_currentSolution, _solutionContainerDC4F);

		for (size_t i=0; i<_norms.size(); i++){
			_resultNormDC4F[i][1]=_norms[i]->compute(_tC, _solutionContainerDC4F, _solutionContainerDC4F);
		}

		if(_exportData.exporter != nullptr && ((1 + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC4F_" + std::to_string(1+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  }


	  _solutionContainerBDF1->invRotate();
		_solutionContainerDC2F->invRotate();
		_solutionContainerDC3F->invRotate();

		// exit(-1);

	} 

	else {
		feErrorMsg(FE_STATUS_ERROR, "Unknown startCode for DC3F");
	}

	printf("=========================================\n");
	printf("          End Starting DC4F              \n");
	printf("=========================================\n");
}


void feTimeIntegratorDC5F::startDC5F(std::string startCode, std::string codeIniDC)
{
	printf("=========================================\n");
	printf("             Starting DC5F               \n");
	printf("=========================================\n");
	
	if (startCode=="Analytique"){
		for(int n=0; n<3; n++){
			this->updateTime(n);

			_solutionContainerBDF1->rotate();
			_solutionContainerDC2F->rotate();
			_solutionContainerDC3F->rotate();
			_solutionContainerDC4F->rotate();
			_solutionContainerDC5F->rotate();

			this->initializeDC5F();

			if(_exportData.exporter != nullptr && ((n+1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC5F_" + std::to_string(n+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  	}

		}
	
	} else if (startCode=="Consistant"){
		printf("=========================================\n");
		printf("       Consistant Starting DC4F          \n");
		printf("=========================================\n");

		printf("=========================================\n");
		printf("       Consistant Starting DC3F          \n");
		printf("=========================================\n");
		for (int n=0; n<2; n++){
			printf("=========================================\n");
			this->updateTime(n);

			_solutionContainerBDF1->rotate();
			_solutionContainerDC2F->rotate();

			_solutionContainerBDF1->computeCoeffBDF(_k);
			_solutionContainerDC2F->computeCoeffBDF(_k);

			feInfo("Solving BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);
			_solutionContainerBDF1->initializeBDFSolution(_mesh, _metaNumber, _currentSolution);

			_solutionContainerBDF1->computeCorrection();

			_solver->solve(_currentSolution, _solutionContainerBDF1);

			if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF1);
			}

			printf("----------------------------------------\n");
			feInfo("Solving DC2/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

			_solutionContainerDC2F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

			_solutionContainerDC2F->computeCorrection(_solutionContainerBDF1, _k);
			
			_solver->solve(_currentSolution, _solutionContainerDC2F);

			if (_correctionType=="RESIDUAL"){
				this->pstClc(_solutionContainerDC2F);
			}

			for (size_t i=0; i<_norms.size(); i++){
				_resultNormBDF1[i][n]=_norms[i]->compute(_tC, _solutionContainerBDF1, _solutionContainerDC4F);
				_resultNormDC2F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC2F, _solutionContainerDC4F);
			}

		}

		_solutionContainerDC3F->computeStartCorrection(_solutionContainerDC2F, _k);
		this->rebootTime();

		this->updateTime(0);

		_solutionContainerDC3F->rotate();
		_solutionContainerDC3F->computeCoeffBDF(_k);

		printf("----------------------------------------\n");
		feInfo("Solving DC3/BDF1 - Step = %d/%d  (%f s)", 0+1, _nTimeStep, _tC);

		_solutionContainerDC3F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);
		
		_solver->solve(_currentSolution, _solutionContainerDC3F);


		for (size_t i=0; i<_norms.size(); i++){
				_resultNormDC3F[i][0]=_norms[i]->compute(_tC, _solutionContainerDC3F, _solutionContainerDC3F);
				// feInfo("normDC3 : %10.10f",_resultNormDC3F[i][0]);
			}

		if(_exportData.exporter != nullptr && ((0 + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC3F_" + std::to_string(0+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
  	}


  	_solutionContainerBDF1->invRotate();
		_solutionContainerDC2F->invRotate();	
		
		printf("=========================================\n");
		printf("     End Consistant Starting DC3F        \n");
		printf("=========================================\n");

		for (int n=1; n<3; n++){
			printf("=========================================\n");
			this->updateTime(n);

			_solutionContainerBDF1->rotate();
			_solutionContainerDC2F->rotate();
			_solutionContainerDC3F->rotate();

			_solutionContainerBDF1->computeCoeffBDF(_k);
			_solutionContainerDC2F->computeCoeffBDF(_k);
			_solutionContainerDC3F->computeCoeffBDF(_k);

			feInfo("Solving BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);
			_solutionContainerBDF1->initializeBDFSolution(_mesh, _metaNumber, _currentSolution);

			_solutionContainerBDF1->computeCorrection();

			_solver->solve(_currentSolution, _solutionContainerBDF1);

			if (_correctionType=="RESIDUAL"){
			this->pstClc(_solutionContainerBDF1);
			}

			printf("----------------------------------------\n");
			feInfo("Solving DC2/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

			_solutionContainerDC2F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

			_solutionContainerDC2F->computeCorrection(_solutionContainerBDF1, _k);
			
			_solver->solve(_currentSolution, _solutionContainerDC2F);

			if (_correctionType=="RESIDUAL"){
				this->pstClc(_solutionContainerDC2F);
			}

			printf("----------------------------------------\n");
			feInfo("Solving DC3/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

			_solutionContainerDC3F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

			_solutionContainerDC3F->computeCorrection(_solutionContainerDC2F, _k);
			
			_solver->solve(_currentSolution, _solutionContainerDC3F);

			if (_correctionType=="RESIDUAL"){
				this->pstClc(_solutionContainerDC3F);
			}


			for (size_t i=0; i<_norms.size(); i++){
				_resultNormBDF1[i][n]=_norms[i]->compute(_tC, _solutionContainerBDF1, _solutionContainerDC4F);
				_resultNormDC2F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC2F, _solutionContainerDC4F);
				_resultNormDC3F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC3F, _solutionContainerDC4F);
			}

		}

		_solutionContainerDC4F->computeStartCorrection(_solutionContainerDC3F, _k);
		this->rebootTime();

		printf("----------------------------------------\n");
		feInfo("Solving DC4/BDF1 - Step = %d/%d  (%f s)", 0+1, _nTimeStep, _tC);
		this->updateTime(0);

		_solutionContainerDC4F->rotate();

		_solutionContainerDC4F->computeCoeffBDF(_k);

		_solutionContainerDC4F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC4F->setStartCorrection(0);

		_solver->solve(_currentSolution, _solutionContainerDC4F);

		for (size_t i=0; i<_norms.size(); i++){
			_resultNormDC4F[i][0]=_norms[i]->compute(_tC, _solutionContainerDC4F, _solutionContainerDC4F);
		}

		if(_exportData.exporter != nullptr && ((0 + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC4F_" + std::to_string(0+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  }

	  printf("----------------------------------------\n");
		feInfo("Solving DC4/BDF1 - Step = %d/%d  (%f s)", 1+1, _nTimeStep, _tC);
		this->updateTime(1);

		_solutionContainerDC4F->rotate();

		_solutionContainerDC4F->computeCoeffBDF(_k);

		_solutionContainerDC4F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC4F->setStartCorrection(1);

		_solver->solve(_currentSolution, _solutionContainerDC4F);

		for (size_t i=0; i<_norms.size(); i++){
			_resultNormDC4F[i][1]=_norms[i]->compute(_tC, _solutionContainerDC4F, _solutionContainerDC4F);
		}

		if(_exportData.exporter != nullptr && ((1 + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC4F_" + std::to_string(1+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  }


	  _solutionContainerBDF1->invRotate();
		_solutionContainerDC2F->invRotate();
		_solutionContainerDC3F->invRotate();

		printf("=========================================\n");
		printf("     End Consistant Start DC4F           \n");
		printf("=========================================\n");





		for (int n=2; n<4; n++){
			printf("=========================================\n");
			this->updateTime(n);

			_solutionContainerBDF1->rotate();
			_solutionContainerDC2F->rotate();
			_solutionContainerDC3F->rotate();
			_solutionContainerDC4F->rotate();
			_solutionContainerDC5F->rotate();

			_solutionContainerBDF1->computeCoeffBDF(_k);
			_solutionContainerDC2F->computeCoeffBDF(_k);
			_solutionContainerDC3F->computeCoeffBDF(_k);
			_solutionContainerDC4F->computeCoeffBDF(_k);
			_solutionContainerDC5F->computeCoeffBDF(_k);
			

			feInfo("Solving BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

			_solutionContainerBDF1->initializeBDFSolution(_mesh, _metaNumber, _currentSolution);
			
			_solutionContainerBDF1->computeCorrection();

			_solver->solve(_currentSolution, _solutionContainerBDF1);

			if (_correctionType=="RESIDUAL"){
				this->pstClc(_solutionContainerBDF1);
			}

			printf("----------------------------------------\n");
			feInfo("Solving DC2/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

			_solutionContainerDC2F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

			_solutionContainerDC2F->computeCorrection(_solutionContainerBDF1, _k);
			
			_solver->solve(_currentSolution, _solutionContainerDC2F);

			if (_correctionType=="RESIDUAL"){
				this->pstClc(_solutionContainerDC2F);
			}

			printf("----------------------------------------\n");
			feInfo("Solving DC3/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

			_solutionContainerDC3F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

			_solutionContainerDC3F->computeCorrection(_solutionContainerDC2F, _k);
			
			_solver->solve(_currentSolution, _solutionContainerDC3F);

			if (_correctionType=="RESIDUAL"){
				this->pstClc(_solutionContainerDC3F);
			}

			printf("----------------------------------------\n");
			feInfo("Solving DC4/BDF1 - Step = %d/%d  (%f s)", n+1, _nTimeStep, _tC);

			_solutionContainerDC4F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

			_solutionContainerDC4F->computeCorrection(_solutionContainerDC3F, _k);
			
			_solver->solve(_currentSolution, _solutionContainerDC4F);

			if (_correctionType=="RESIDUAL"){
				this->pstClc(_solutionContainerDC4F);
			}


			for (size_t i=0; i<_norms.size(); i++){
			_resultNormBDF1[i][n]=_norms[i]->compute(_tC, _solutionContainerBDF1, _solutionContainerDC5F);
			_resultNormDC2F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC2F, _solutionContainerDC5F);
			_resultNormDC3F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC3F, _solutionContainerDC5F);
			_resultNormDC4F[i][n]=_norms[i]->compute(_tC, _solutionContainerDC4F, _solutionContainerDC5F);
			}
		}

		_solutionContainerDC5F->computeStartCorrection(_solutionContainerDC4F, _k);
		this->rebootTime();

		printf("----------------------------------------\n");
		feInfo("Solving DC5/BDF1 - Step = %d/%d  (%f s)", 0+1, _nTimeStep, _tC);
		this->updateTime(0);

		_solutionContainerDC5F->rotate();

		_solutionContainerDC5F->computeCoeffBDF(_k);

		_solutionContainerDC5F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC5F->setStartCorrection(0);

		_solver->solve(_currentSolution, _solutionContainerDC5F);

		for (size_t i=0; i<_norms.size(); i++){
			_resultNormDC5F[i][0]=_norms[i]->compute(_tC, _solutionContainerDC5F, _solutionContainerDC5F);
		}

		if(_exportData.exporter != nullptr && ((0 + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC5F_" + std::to_string(0+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  }


	  printf("----------------------------------------\n");
		feInfo("Solving DC5/BDF1 - Step = %d/%d  (%f s)", 1+1, _nTimeStep, _tC);
		this->updateTime(1);

		_solutionContainerDC5F->rotate();

		_solutionContainerDC5F->computeCoeffBDF(_k);

		_solutionContainerDC5F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC5F->setStartCorrection(1);

		_solver->solve(_currentSolution, _solutionContainerDC5F);

		for (size_t i=0; i<_norms.size(); i++){
			_resultNormDC5F[i][1]=_norms[i]->compute(_tC, _solutionContainerDC5F, _solutionContainerDC5F);
		}

		if(_exportData.exporter != nullptr && ((1 + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC5F_" + std::to_string(1+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  }

	  printf("----------------------------------------\n");
		feInfo("Solving DC5/BDF1 - Step = %d/%d  (%f s)", 2+1, _nTimeStep, _tC);
		this->updateTime(2);

		_solutionContainerDC5F->rotate();

		_solutionContainerDC5F->computeCoeffBDF(_k);

		_solutionContainerDC5F->initializeDCSolution(_mesh, _metaNumber, _currentSolution, codeIniDC);

		_solutionContainerDC5F->setStartCorrection(2);

		_solver->solve(_currentSolution, _solutionContainerDC5F);

		for (size_t i=0; i<_norms.size(); i++){
			_resultNormDC5F[i][2]=_norms[i]->compute(_tC, _solutionContainerDC5F, _solutionContainerDC5F);
		}

		if(_exportData.exporter != nullptr && ((2 + 1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC5F_" + std::to_string(2+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  }

	  _solutionContainerBDF1->invRotate();
		_solutionContainerDC2F->invRotate();
		_solutionContainerDC3F->invRotate();
		_solutionContainerDC4F->invRotate();

		// exit(-1);
		
	} else {
		feErrorMsg(FE_STATUS_ERROR, "Unknown startCode for DC5F");
	}

	printf("=========================================\n");
	printf("          End Starting DC5F              \n");
	printf("=========================================\n");
	
}


void feTimeIntegratorDC3::startDC3(std::string startCode)
{
	printf("=========================================\n");
	printf("             Starting DC3                \n");
	printf("=========================================\n");
	
	if (startCode=="Analytique"){
		this->updateTime(0);

		_solutionContainerBDF2->rotate();
		_solutionContainerDC3->rotate();

		this->initializeDC3();

		if(_exportData.exporter != nullptr && ((1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC3_" + std::to_string(1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  	}
	
	} else if (startCode=="Consistant"){
		feInfo("TODO : demarrage consistant DC3");
	
	} else {
		feErrorMsg(FE_STATUS_ERROR, "Unknown startCode for DC3");
	}

	printf("=========================================\n");
	printf("          End Starting DC3               \n");
	printf("=========================================\n");
	
}


void feTimeIntegratorDC4::startDC4(std::string startCode)
{
	printf("=========================================\n");
	printf("             Starting DC4                \n");
	printf("=========================================\n");
	
	if (startCode=="Analytique"){
		for (int n=0; n<2; n++){
			this->updateTime(n);

			_solutionContainerBDF2->rotate();
			_solutionContainerDC3->rotate();
			_solutionContainerDC4->rotate();

			this->initializeDC4();

			if(_exportData.exporter != nullptr && ((n+1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC4_" + std::to_string(n+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  	}
		}
		
	
	} else if (startCode=="Consistant"){
		feInfo("TODO : demarrage consistant DC4");
	
	} else {
		feErrorMsg(FE_STATUS_ERROR, "Unknown startCode for DC4");
	}

	printf("=========================================\n");
	printf("          End Starting DC4               \n");
	printf("=========================================\n");
	
}

void feTimeIntegratorDC5::startDC5(std::string startCode)
{
	printf("=========================================\n");
	printf("             Starting DC5                \n");
	printf("=========================================\n");
	
	if (startCode=="Analytique"){
		for (int n=0; n<3; n++){
			this->updateTime(n);

			_solutionContainerBDF2->rotate();
			_solutionContainerDC3->rotate();
			_solutionContainerDC4->rotate();
			_solutionContainerDC5->rotate();

			this->initializeDC5();

			if(_exportData.exporter != nullptr && ((n+1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC5_" + std::to_string(n+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  	}
		}
		
	
	} else if (startCode=="Consistant"){
		feInfo("TODO : demarrage consistant DC5");
	
	} else {
		feErrorMsg(FE_STATUS_ERROR, "Unknown startCode for DC5");
	}

	printf("=========================================\n");
	printf("          End Starting DC5               \n");
	printf("=========================================\n");
	
}


void feTimeIntegratorDC6::startDC6(std::string startCode)
{
	printf("=========================================\n");
	printf("             Starting DC6                \n");
	printf("=========================================\n");
	
	if (startCode=="Analytique"){
		for (int n=0; n<4; n++){
			this->updateTime(n);

			_solutionContainerBDF2->rotate();
			_solutionContainerDC3->rotate();
			_solutionContainerDC4->rotate();
			_solutionContainerDC5->rotate();
			_solutionContainerDC6->rotate();

			this->initializeDC6();

			if(_exportData.exporter != nullptr && ((n+1) % _exportData.exportEveryNSteps) == 0) {
        std::string fileName = _exportData.fileNameRoot + "_DC6_" + std::to_string(n+1) + ".vtk";
        _exportData.exporter->writeStep(fileName);
	  	}
		}
		
	
	} else if (startCode=="Consistant"){
		feInfo("TODO : demarrage consistant DC6");
	
	} else {
		feErrorMsg(FE_STATUS_ERROR, "Unknown startCode for DC6");
	}

	printf("=========================================\n");
	printf("          End Starting DC6               \n");
	printf("=========================================\n");
	
}

//=======================================================================================
// Time rotation and update
//=======================================================================================

void feTimeIntegrator::initializeTime()
{
	_tC = _t0;

	if (_codeDt == "Constant") {
		_dt.resize(1);
		_dt[0] = (_tEnd - _t0)/_nTimeStep;
	
	}else if (_codeDt =="Alt") {
		if (_nTimeStep%2 != 0) {
			feWarning("Pour pas de temps alternée necessité d'avoir un nombre paire de nTimeStep ==> nTimeSteps +1");
			_nTimeStep += 1;
		}
		int ratio = 4;
		double dt2 = (_tEnd-_t0)/_nTimeStep *2/(ratio +1);
		double dt1 = ratio*dt2;
		_dt.resize(2);
		_dt = {dt1, dt2};

	}else if (_codeDt =="AltV2") {
		int r = _nTimeStep%3;
		if ( r != 0) {
			feWarning("nTimeStep not 3 multiple ==> +/- 1 Step");
			if (r==1){
				_nTimeStep -= 1;
			} else if (r==2){
				_nTimeStep += 1;
			}
			
		}
		double k1k2Ratio = 4;
		double k2k3Ratio = 0.4;
		double dt2 = (_tEnd-_t0)/_nTimeStep *3/(k1k2Ratio + 1 + 1/k2k3Ratio);
		double dt1 = k1k2Ratio*dt2;
		double dt3 = dt2/k2k3Ratio;
		_dt.resize(3);
		_dt = {dt1, dt2, dt3};

	} else if(_codeDt == "Increase"){
		_gamma = 2.;
		_r = pow(_gamma,1./(_nTimeStep - 1.));
		_k0 = _tEnd*(_r - 1.)/(pow(_r,_nTimeStep) - 1.);
		_dt.resize(1);
		_dt[0] = _k0;


	}else if (_codeDt == "InDecrease"){
		if (_nTimeStep%2 != 0){
			feInfo("nTimeStep not even => adding 1 step");
			_nTimeStep += 1;
		}
		double tM = (_tEnd + _t0)/2.;
		_nTimeStepM = _nTimeStep/2.;

		_gamma = 2.;
		_r = pow(_gamma, 1./(_nTimeStepM - 1.));
		_k0 = tM*(_r - 1.)/(pow(_r,_nTimeStepM) - 1.);
		_dt.resize(1);
		_dt[0] = _k0;


	}else if(_codeDt == "InDecreaseV2"){
		_nbInterv = 10;
		_nTimeStepPerInterv = _nTimeStep;

		if (_nTimeStepPerInterv%2 !=0){
			feWarning("nTimeStep not even => adding 1 step");
			_nTimeStepPerInterv += 1;
		}

		_nTimeStep = _nTimeStepPerInterv * _nbInterv;
		
		double I = (_tEnd - _t0)/_nbInterv;
		_nTimeStepM = _nTimeStepPerInterv/2;
		double tF1 = _t0 + I;
		double tM = (tF1+_t0)/2.;

		_gamma = 2;
		_r = pow(_gamma, 1./(_nTimeStepM - 1.));
		_k0 = tM*(_r - 1.)/(pow(_r,_nTimeStepM) - 1.);
		_dt.resize(1);
		_dt[0] = _k0;


	}else{
		feWarning("codeDt non reconnu ==> par defaut dt=Constant");
		_dt.resize(1);
		_dt[0] = (_tEnd - _t0)/_nTimeStep;
	}

	_tG.resize(_nTimeStep);
}



void feTimeIntegrator::updateTime(int n)
{
	for (int i=_nStorage-1; i>0; i--){
		_t[i]=_t[i-1];
		_k[i]=_k[i-1];
	}

	if (_codeDt=="Alt"){
		double dtS = _dt[0];
		_dt[0] = _dt[1];
		_dt[1] = dtS;

	
	} else if (_codeDt =="AltV2"){
		double dtS = _dt[0];
		_dt[0] = _dt[1];
		_dt[1] = _dt[2];
		_dt[2] = dtS;		


	} else if (_codeDt == "Increase"){
		_dt[0] = pow(_r,n)*_k0;
	
	} else if (_codeDt =="InDecrease"){
		int exp=n;
		if (n>_nTimeStepM-1){
			exp = _nTimeStep -1 - n;
		}

		_dt[0] = pow(_r,exp)*_k0;
	
	} else if (_codeDt == "InDecreaseV2"){
		int nM = n - n/_nTimeStepPerInterv * _nTimeStepPerInterv;

		int exp= nM;
		if (nM>_nTimeStepM-1){
			exp = _nTimeStepPerInterv -1 - nM;
		}

		_dt[0] = pow(_r,exp)*_k0;
	}

	_tC += _dt[0];
	_t[0]= _tC;
	_k[0]=_t[0]-_t[1];
	_tG[n] = _tC;	

	_currentSolution->setCurrentTime(_tC);
}


void feTimeIntegrator::rebootTime()
{
	// this->initializeTime();
	_tC = _t0;

	for (size_t i=0; i<_t.size(); i++){
		_t[i] =0.0;
	}

	for (size_t i=0; i<_k.size(); i++){
		_k[i] =0.0;
	}

	_t[0] = _t0;
}