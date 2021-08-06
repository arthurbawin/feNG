#include "feSolution.h"

void feSolution::initializeTemporalSolution(double t0, double t1, int nTimeSteps) {
  _t0 = t0;
  _t1 = t1;
  _nTimeSteps = nTimeSteps;
  _dt = (t1 - t0) / (double)nTimeSteps;
  _tn = t0;
}

void feSolution::initializeUnknowns(feMesh *mesh, feMetaNumber *metaNumber) {
  for(feSpace *const &fS : _space) {
    // std::cout<<"Space "<<fS->getFieldID()<<" - "<<fS->getCncGeoID()<<std::endl;
    int nElm = mesh->getNbElm(fS->getCncGeoID()); // On pourrait donner un _nElm au feSpace
    std::vector<double> coor = fS->getLcoor();
    // std::cout<<"coor"<<std::endl;
    // for(auto const &val : coor)
    // 	std::cout<<val<<std::endl;
    feSpace *geoSpace = mesh->getGeometricSpace(fS->getCncGeoID());
    // std::cout<<geoSpace->getNbFunctions()<<"?"<<std::endl;
    for(int iElm = 0; iElm < nElm; ++iElm) {
      // Call initializeAddressingVector() on each element of each feSpace
      fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), iElm);

      std::vector<double> localCoord = mesh->getCoord(fS->getCncGeoID(), iElm);
      // std::cout<<"size"<<localCoord.size()<<std::endl;

      for(int j = 0; j < fS->getNbFunctions(); ++j) {
        double r[3] = {coor[3 * j], coor[3 * j + 1], coor[3 * j + 2]};
        std::vector<double> x(3, 0.0);
        // x[0] = geoSpace->interpolateField(localCoord, r);
        geoSpace->interpolateVectorField(localCoord, r, x);
        // std::cout<<x[0]<<std::endl;
        double val = fS->evalFun(_tn, x);
        _sol[fS->getAddressingVectorAt(j)] = val;
      }
    }
  }
}

void feSolution::initializeEssentialBC(feMesh *mesh, feMetaNumber *metaNumber,
                                       feSolutionContainer *solContainer) {
  for(feSpace *const &fS : _essBC) {
    int nElm = mesh->getNbElm(fS->getCncGeoID());
    std::vector<double> coor = fS->getLcoor();
    feSpace *geoSpace = mesh->getGeometricSpace(fS->getCncGeoID());
    for(int iElm = 0; iElm < nElm; ++iElm) {
      fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), iElm);
      std::vector<double> localCoord = mesh->getCoord(fS->getCncGeoID(), iElm);
      for(int j = 0; j < fS->getNbFunctions(); ++j) {
        double r[3] = {coor[3 * j], coor[3 * j + 1], coor[3 * j + 2]};
        std::vector<double> x(3, 0.0);
        geoSpace->interpolateVectorField(localCoord, r, x);
        _sol[fS->getAddressingVectorAt(j)] = fS->evalFun(_tn, x);
        if(solContainer != nullptr)
          solContainer->_sol[0][fS->getAddressingVectorAt(j)] = fS->evalFun(_tn, x);
        ;
      }
    }
  }
}

void feSolution::setSolFromContainer(feSolutionContainer *solContainer, int iSol) {
  size_t nDOFs = _sol.size();
  for(size_t i = 0; i < nDOFs; ++i) _sol[i] = solContainer->_sol[iSol][i];
}

void feSolution::setSolDotToZero() {
  size_t nDOFs = _dsoldt.size();
  for(size_t i = 0; i < nDOFs; ++i) _dsoldt[i] = 0.0;
}

void feSolution::printSol() {
  for(auto const &val : _sol) printf("%12.16f\n", val);
}
