#include "feSpaceTriangle.h"
#include "feMesh.h"
#include "feNumber.h"

std::vector<double> feSpaceTriP1::L(double r[3]){ 
  return { 1.0 - r[0] - r[1], r[0], r[1] }; 
};

std::vector<double> feSpaceTriP1::dLdr(double r[3]){ return { -1.0, 1.0, 0.0 }; };
std::vector<double> feSpaceTriP1::dLds(double r[3]){ return { -1.0, 0.0, 1.0 }; };
std::vector<double> feSpaceTriP1::dLdt(double r[3]){ return {   0.,  0.,  0. }; };

void feSpaceTriP1::initializeNumberingUnknowns(feNumber *number){
	for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i){
    number->defDDLSommet(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet(_mesh, _cncGeoID, i, 1);
    number->defDDLSommet(_mesh, _cncGeoID, i, 2);
  }
}

void feSpaceTriP1::initializeNumberingEssential(feNumber *number){
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i){
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 1);
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 2);
  }
}

void feSpaceTriP1::initializeAddressingVector(feNumber *number, int numElem){
  _adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  _adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
  _adr[2] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 2);
}