#include "feTimeIntegrator.h"

void feTimeIntegrator::initialize(feSolution *sol, feMesh *mesh, feMetaNumber *metaNumber){
	sol->setCurrentTime(_t[0]);
	sol->initializeUnknowns(mesh, metaNumber);
	sol->initializeEssentialBC(mesh, metaNumber);
	for(int i = 0; i < _nDofs; ++i)
		_sol[i] = sol->getSolAtDOF(i);
}