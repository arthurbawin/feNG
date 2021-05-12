#include "feNorm.h"
#include <cmath>
#include "feQuadrature.h"

feNorm::feNorm(feSpace *intSpace, feMesh *mesh, int nQuad)
  : _intSpace(intSpace), cncGeoTag(intSpace->getCncGeoTag()), geoSpace(mesh->getGeometricSpace(cncGeoTag)),
  nElmGeo(mesh->getNbElm(cncGeoTag)), dim(mesh->getDim()), nNodePerElem(intSpace->getNbNodePerElem()),
  _nQuad(nQuad)
{
  feQuadrature2 *rule = new feQuadrature2(_nQuad);
  w = rule->getWeights();
  x = rule->getPoints();
  intSpace->setQuadratureRule(rule);
  geoSpace->setQuadratureRule(rule);
  delete rule;
}

void feNorm::computeL2Norm(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh){
  double normL2 = 0.0, solInt, solRef, xInt, J, t = sol->getCurrentTime();
  int nElm = _intSpace->getNbElm();

  for(int iElm = 0; iElm < nElm; ++iElm){
    _intSpace->initializeAddressingVector(metaNumber->getNumbering(_intSpace->getFieldID()), iElm);
    _intSpace->initializeSolution(sol);
    geoCoord = mesh->getCoord(cncGeoTag, iElm);
    for(int k = 0; k < _nQuad; ++k){
      solInt = _intSpace->interpolateSolutionAtQuadNode(k);
      // x-positions to interpolate
      std::vector<double> xGeoCoord(nNodePerElem,0.0); 
      for(int i = 0; i < nNodePerElem; ++i)
        xGeoCoord[i] = geoCoord[dim*i];

      xInt = geoSpace->interpolateFieldAtQuadNode(xGeoCoord, k);
      J = geoSpace->interpolateFieldAtQuadNode_rDerivative(xGeoCoord, k);

      solRef = _intSpace->evalFun(t, {xInt, 0., 0.});
      normL2 += (solInt - solRef)*(solInt - solRef) * J * w[k];
    }
  }
  norm = sqrt(normL2);
}