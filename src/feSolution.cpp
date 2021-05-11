#include "feSolution.h"

void feSolution::initializeUnknowns(feMesh *mesh, feMetaNumber *metaNumber){
	for(feSpace* const &fS : _space){
		// std::cout<<"Space "<<fS->getFieldID()<<" - "<<fS->getCncGeoID()<<std::endl;
		int nElm = mesh->getNbElm(fS->getCncGeoID()); // On pourrait donner un _nElm au feSpace
		std::vector<double> coor = fS->getLcoor();
		// std::cout<<"coor"<<std::endl;
		// for(auto const &val : coor)
		// 	std::cout<<val<<std::endl;
		feSpace* getGeometricSpace = mesh->getGeometricSpace(fS->getCncGeoID());
		for(int iElm = 0; iElm < nElm; ++iElm){
			// Call initializeAddressingVector() on each element of each feSpace
			fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), iElm);

			std::vector<double> coocnc = mesh->getCoord(fS->getCncGeoID(), iElm);
			// std::cout<<"coocnc"<<std::endl;
			// for(auto const &val : coocnc)
			// 	std::cout<<val<<std::endl;
			// std::cout<<"size"<<coocnc.size()<<std::endl;

			

			for(int j = 0; j < fS->getNbFunctions(); ++j){
				double r[3] = {coor[j],0.,0.};
				// std::cout<<r[0]<<" "<<r[1]<<" "<<r[2]<<	" "<<std::endl;
				std::vector<double> x(3,0.0);
				x[0] = getGeometricSpace->interpolateField(coocnc, r);
				// std::cout<<x[0]<<std::endl;
				double val = fS->evalFun(_tn, x);
				_sol[fS->getAddressingVectorAt(j)] = val;
			}
		}
	}
}

void feSolution::initializeEssentialBC(feMesh *mesh, feMetaNumber *metaNumber){
	for(feSpace* const &fS : _essBC){
		int nElm = mesh->getNbElm(fS->getCncGeoID()); // On pourrait donner un _nElm au feSpace
		std::vector<double> coor = fS->getLcoor();
		feSpace* getGeometricSpace = mesh->getGeometricSpace(fS->getCncGeoID());
		for(int iElm = 0; iElm < nElm; ++iElm){
			fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), iElm);
			std::vector<double> coocnc = mesh->getCoord(fS->getCncGeoID(), iElm);
			for(int j = 0; j < fS->getNbFunctions(); ++j){
				double r[3] = {coor[j],0.,0.};
				std::vector<double> x(3,0.0);
				x[0] = getGeometricSpace->interpolateField(coocnc, r);
				_sol[fS->getAddressingVectorAt(j)] = fS->evalFun(_tn, x);;
			}
		}
	}
}

void feSolution::printSol(){
	for(auto const &val : _sol)
		std::cout<<val<<std::endl;
}