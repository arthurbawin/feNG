#include "feNorm.h"
#include <cmath>

feNorm::feNorm(feSpace *intSpace, feMesh *mesh, int nQuad)
  : _intSpace(intSpace), cncGeoTag(intSpace->getCncGeoTag()), geoSpace(mesh->getGeometricSpace(cncGeoTag)),
  nElmGeo(mesh->getNbElm(cncGeoTag)), dim(mesh->getDim()), nNodePerElem(intSpace->getNbNodePerElem()),
  _nQuad(nQuad)
{
  feQuadrature *rule = new feQuadrature(nQuad, mesh->getDim());
  w = rule->getWeights();
  // x = rule->getPoints();
  // Attention : ça changle la règle de quadrature pour les interpolants avant la résolution !
  // Il faudrait choisir idépendamment les deux sans conséquences
  intSpace->setQuadratureRule(rule);
  geoSpace->setQuadratureRule(rule);
  delete rule;
}

void feNorm::computeL2Norm(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh){
  double normL2 = 0.0, solInt, solRef, xInt, yInt, zInt, J, t = sol->getCurrentTime();
  int nElm = _intSpace->getNbElm();

  for(int iElm = 0; iElm < nElm; ++iElm){
    _intSpace->initializeAddressingVector(metaNumber->getNumbering(_intSpace->getFieldID()), iElm);
    _intSpace->initializeSolution(sol);
    geoCoord = mesh->getCoord(cncGeoTag, iElm);
    for(int k = 0; k < _nQuad; ++k){
      solInt = _intSpace->interpolateSolutionAtQuadNode(k);
      // // x-positions to interpolate
      // std::vector<double> xGeoCoord(nNodePerElem,0.0);
      // for(int i = 0; i < nNodePerElem; ++i)
      //   xGeoCoord[i] = geoCoord[dim*i];

      // xInt = geoSpace->interpolateFieldAtQuadNode(xGeoCoord, k);
      std::vector<double> x(3,0.0);
      geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      // J = geoSpace->interpolateFieldAtQuadNode_rDerivative(xGeoCoord, k);
      if(dim == 1){
        std::vector<double> j(3, 0.0);
        geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, j); // TODO : complete for dim > 1
        J = j[0];
      } else if(dim == 2){
        std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
        std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
        geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
        geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
        J = dxdr[0]*dxds[1] - dxdr[1]*dxds[0];
      }

      solRef = _intSpace->evalFun(t, x);
      normL2 += (solInt - solRef)*(solInt - solRef) * J * w[k];
    }
  }
  norm = sqrt(normL2);
}

void feNorm::computeArea(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh){
  double area = 0.0, solInt, solRef, xInt, yInt, zInt, J, t = sol->getCurrentTime();
  int nElm = _intSpace->getNbElm();

  for(int iElm = 0; iElm < nElm; ++iElm){
    _intSpace->initializeAddressingVector(metaNumber->getNumbering(_intSpace->getFieldID()), iElm);
    _intSpace->initializeSolution(sol);
    geoCoord = mesh->getCoord(cncGeoTag, iElm);
    for(int k = 0; k < _nQuad; ++k){
      solInt = _intSpace->interpolateSolutionAtQuadNode(k);
      // // x-positions to interpolate
      // std::vector<double> xGeoCoord(nNodePerElem,0.0);
      // for(int i = 0; i < nNodePerElem; ++i)
      //   xGeoCoord[i] = geoCoord[dim*i];

      // xInt = geoSpace->interpolateFieldAtQuadNode(xGeoCoord, k);
      std::vector<double> x(3,0.0);
      geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      // J = geoSpace->interpolateFieldAtQuadNode_rDerivative(xGeoCoord, k);
      // std::vector<double> j(3, 0.0);
      // geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, j); // TODO : complete for dim > 1
      // J = j[0];
      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
      J = dxdr[0]*dxds[1] - dxdr[1]*dxds[0];

      area += J * w[k];
    }
  }
  norm = area;
}