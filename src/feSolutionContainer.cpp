#include "feSolutionContainer.h"

feSolutionContainer::feSolutionContainer(int nSol, double tn, feMetaNumber *metaNumber)
  : _nSol(nSol), _nDofs(metaNumber->getNbDOFs()){
	_t.resize(_nSol);
	_t[0] = tn;
	_sol.resize(_nSol);
	_fResidual.resize(_nSol);
	_cn.resize(3);
	_d.resize(_nDofs);
}

void feSolutionContainer::initialize(feSolution *sol, feMesh *mesh, feMetaNumber *metaNumber){
	sol->setCurrentTime(_t[0]);
	sol->initializeUnknowns(mesh, metaNumber);
	sol->initializeEssentialBC(mesh, metaNumber);
	_sol.push_back(sol->getSolutionCopy());
}

void feStationarySolution::computeSolTimeDerivative(feSolution *sol){
  for(int i = 0; i < _nDofs; ++i)
    sol->setSolDotAtDOF(i, 0.0);
}

void feSolutionBDF2::computeSolTimeDerivative(feSolution *sol){
  for(int i = 0; i < _nDofs; ++i)
    sol->setSolDotAtDOF(i, 0.0);
}

void feSolutionDCF::computeSolTimeDerivative(feSolution *sol){
  for(int i = 0; i < _nDofs; ++i)
    sol->setSolDotAtDOF(i, 0.0);
}