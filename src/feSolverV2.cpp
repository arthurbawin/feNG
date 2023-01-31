#include "feSolverV2.h"
// #include "feSolutionContainer.h"
#include "feExporter.h"
#include "feNG.h"

extern int FE_VERBOSE;

feStatus createSolver(feSolverV2 *&solver, feLinearSystem *linearSystem, feMetaNumber *metaNumber, solverType codeSolver, feTolerancesV2 tol, feMatrixComputeCriteria *matrixComputeCriteria)
{ 
    switch(codeSolver) {
        case Newton :
            solver = new feNewtonSolver(linearSystem, metaNumber, codeSolver, tol, matrixComputeCriteria); 
            break;

        default:
            return feErrorMsg(FE_STATUS_ERROR, "Unsupported solver type");
    }
    return FE_STATUS_OK;
}





void feNewtonSolver::solve(feSolution *currentSolution, feSolutionContainerV2 *solutionContainer)
{
    int iter = 0, linearSystemIter;
    double normDx, normResidual, normAxb;
    
    bool newton=true;
    while(newton) {
        _linearSystem->setToZero();

        solutionContainer->computeSolDot(_linearSystem);

        currentSolution->setSolDotFromContainer(solutionContainer, 0);
        
        _linearSystem->assemble(currentSolution);
        _linearSystem->solve(&normDx, &normResidual, &normAxb, &linearSystemIter);
        
        _linearSystem->correctSolution(currentSolution);

        solutionContainer->copyCurrentSolution(currentSolution);

        std::string status = _linearSystem->getRecomputeStatus() ? "true" : "false";

        feInfoCond(FE_VERBOSE > 0, "    iter %2d : ||A*dx-res|| = %10.10e (%4d iter.) \t ||dx|| = %10.10e \t ||res|| = %10.10e  (%s)", ++iter, normAxb, linearSystemIter, normDx, normResidual, status.c_str());
        
        newton = !((normDx <= _tol.tolDx && normResidual <= _tol.tolResidual) || iter > _tol.maxIter);

        if (_matrixComputeCriteria != nullptr)
            this->setMatrixComputeStatus(newton, iter, normDx, normResidual);

    }

    if(iter > _tol.maxIter) {
        feWarning("Not converged at iter %2d : ||A*dx-res|| = %10.10e  (%4d iter.) \t ||dx|| = "
                  "%10.10e \t ||res|| = %10.10e",
                  iter, normAxb, linearSystemIter, normDx, normResidual);
        exit(0);

    } else {
        feInfoCond(FE_VERBOSE > 0,
                   "Converged in %2d Newton iterations : ||dx|| = %10.10e \t ||res|| = %10.10e", iter,
                   normDx, normResidual);
    }
}






void feSolverV2::setMatrixComputeStatus(bool newton, int iter, double normDx, double normResidual)
{
    if (newton == false){
        if (iter<_matrixComputeCriteria->convergedBeforeNIter){
            // feInfo("coucou5");
            _linearSystem->setRecomputeStatus(false);
        
        } else {
            // feInfo("coucou6");
            _linearSystem->setRecomputeStatus(true);
        }
    
    } else {
        if ((normDx < _matrixComputeCriteria->criteriaNormU) && (normResidual < _matrixComputeCriteria->criteriaNormRes)){
            if (iter == 1){
                // feInfo("coucou1");
                _linearSystem->setRecomputeStatus(false);
                _matrixComputeCriteria->cntBeforRecompute += 1;
                _matrixComputeCriteria->perviousNormU = normDx;
                _matrixComputeCriteria->perviousNormRes = normResidual;
            
            } else if ((normDx/_matrixComputeCriteria->perviousNormU) < _matrixComputeCriteria->ratioNormU && (normResidual/_matrixComputeCriteria->perviousNormRes)<_matrixComputeCriteria->ratioNormRes && (_matrixComputeCriteria->cntBeforRecompute < _matrixComputeCriteria->nbIterBeforeRecompute)){
                // feInfo("coucou2");
                _linearSystem->setRecomputeStatus(false);
                _matrixComputeCriteria->cntBeforRecompute += 1;
                _matrixComputeCriteria->perviousNormU   = normDx;
                _matrixComputeCriteria->perviousNormRes = normResidual;
            
            } else {
                // feInfo("coucou3");
                _linearSystem->setRecomputeStatus(true);
                _matrixComputeCriteria->cntBeforRecompute = 0;
            }
        
        } else {
            // feInfo("coucou4");
            _matrixComputeCriteria->cntBeforRecompute = 0;
            _linearSystem->setRecomputeStatus(true);
        }
    }

}



void feNewtonSolver::solveUDOT(feSolution *currentSolution)
{
    int iter = 0, linearSystemIter;
    double normDx, normResidual, normAxb;

    currentSolution->setC0(1.);
    
    bool newton=true;
    while(newton) {
        _linearSystem->setToZero();

        // solutionContainer->computeSolDot(_linearSystem);

        // currentSolution->setSolDotFromContainer(solutionContainer, 0);
        
        _linearSystem->assembleIni(currentSolution);
        _linearSystem->solve(&normDx, &normResidual, &normAxb, &linearSystemIter);
        
        _linearSystem->correctSolutionDot(currentSolution);

        std::string status = _linearSystem->getRecomputeStatus() ? "true" : "false";

        feInfoCond(FE_VERBOSE > 0, "    iter %2d : ||A*dx-res|| = %10.10e (%4d iter.) \t ||dx|| = %10.10e \t ||res|| = %10.10e  (%s)", ++iter, normAxb, linearSystemIter, normDx, normResidual, status.c_str());
        
        newton = !((normDx <= _tol.tolDx && normResidual <= _tol.tolResidual) || iter > _tol.maxIter);
        // _linearSystem->setRecomputeStatus(false);
        // if (_matrixComputeCriteria != nullptr)
        //     this->setMatrixComputeStatus(newton, iter, normDx, normResidual);

    }

    if(iter > _tol.maxIter) {
        feWarning("Not converged at iter %2d : ||A*dx-res|| = %10.10e  (%4d iter.) \t ||dx|| = "
                  "%10.10e \t ||res|| = %10.10e",
                  iter, normAxb, linearSystemIter, normDx, normResidual);
        exit(0);

    } else {
        feInfoCond(FE_VERBOSE > 0,
                   "Converged in %2d Newton iterations : ||dx|| = %10.10e \t ||res|| = %10.10e", iter,
                   normDx, normResidual);
    }
}