#ifndef _FESOLVERV2_
#define _FESOLVERV2_

#include "feMessage.h"
#include "feMesh.h"
#include "feLinearSystem.h"
#include "feBilinearForm.h"
#include "feNumber.h"
// #include "feNorm.h"
#include "feSpace.h"
#include "feExporter.h"
// #include "feComputer.h"
#include "feSolutionContainerV2.h"
#include "feNG.h"


typedef enum {Newton} solverType;

typedef struct feTolV2 {
  double tolDx;
  double tolResidual;
  double maxIter;
} feTolerancesV2;


typedef struct feMatComputeCriteria{
	double criteriaNormU;
	double criteriaNormRes;
	double ratioNormU;
	double ratioNormRes;
	int nbIterBeforeRecompute;
	int convergedBeforeNIter;
	int cntBeforRecompute  = 0;
	double perviousNormU   = 0.;
	double perviousNormRes = 0.;
} feMatrixComputeCriteria;


class feSolverV2 {
protected:
	feLinearSystem *_linearSystem;
	feTolerancesV2 _tol;
	feMatrixComputeCriteria *_matrixComputeCriteria;
	solverType _solverType;

	int _nbInc;
	

public:
	feSolverV2(feLinearSystem *linearSystem, feMetaNumber *metaNumber, solverType codeSolver, feTolerancesV2 tol, feMatrixComputeCriteria *matrixComputeCriteria=nullptr):  //, 
		_linearSystem(linearSystem), _tol(tol), _solverType(codeSolver), _nbInc(metaNumber->getNbUnknowns()), _matrixComputeCriteria(matrixComputeCriteria) {} //, 

	virtual void solve(feSolution *currentSolution, feSolutionContainerV2 *solutionContainer){};
	virtual void solveUDOT(feSolution *currentSolution){};

	void setMatrixComputeStatus(bool Newton, int iter, double normDx, double normResidual);

	virtual ~feSolverV2() {}

};



class feNewtonSolver : public feSolverV2 {
protected:

public:
	feNewtonSolver(feLinearSystem *linearSystem, feMetaNumber *metaNumber,  solverType codeSolver, feTolerancesV2 tol, feMatrixComputeCriteria *matrixComputeCriteria)
		:feSolverV2(linearSystem, metaNumber, codeSolver, tol, matrixComputeCriteria){}

	virtual void solve(feSolution *currentSolution, feSolutionContainerV2 *solutionContainer);
	virtual void solveUDOT(feSolution *currentSolution);

	virtual ~feNewtonSolver() {}
};



feStatus createSolver(feSolverV2 *&solver, feLinearSystem *linearSystem, feMetaNumber *metaNumber, solverType codeSolver, feTolerancesV2 tol, feMatrixComputeCriteria *matrixComputeCriteria);

#endif