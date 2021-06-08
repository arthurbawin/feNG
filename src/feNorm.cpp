#include "feNorm.h"
#include <cmath>
#include "feQuadrature.h"


// feNorm::feNorm(feSpace *intSpace, feMesh *mesh, int degQuad)
//   : _intSpace(intSpace), cncGeoTag(intSpace->getCncGeoTag()), geoSpace(mesh->getGeometricSpace(cncGeoTag)),
//   nElmGeo(mesh->getNbElm(cncGeoTag)), dim(mesh->getDim()), nNodePerElem(intSpace->getNbNodePerElem()),
//   _degQuad(degQuad)

feNorm::feNorm(feSpace *intSpace, feMesh *mesh, int nQuad, feFunction *solRef)
  : _intSpace(intSpace), cncGeoTag(intSpace->getCncGeoTag()), geoSpace(mesh->getGeometricSpace(cncGeoTag)),
  nElmGeo(mesh->getNbElm(cncGeoTag)), dim(mesh->getDim()), nNodePerElem(intSpace->getNbNodePerElem()),
  _nQuad(nQuad), _solRef(solRef)

{

  feQuadrature *rule = new feQuadrature(_degQuad,_dimQuad, form);
  w = rule->getWeights();
  x = rule->getXPoints();
  _nQuad=rule->getNQuad();

  // x = rule->getPoints();
  y = rule->getYPoints();
  z = rule->getZPoints();
  // Attention : ça changle la règle de quadrature pour les interpolants avant la résolution !
  // Il faudrait choisir idépendamment les deux sans conséquences

  intSpace->setQuadratureRule(rule);
  geoSpace->setQuadratureRule(rule);
  delete rule;
}

void feNorm::computeL2Norm(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh){
  double normL2 = 0.0, solInt, solRef, xInt, yInt, zInt, J, t = sol->getCurrentTime();
  int nElm = _intSpace->getNbElm();

  // std::cout<<"Computing norm of "<<_intSpace->getFieldID()<<" on "<<nElm<<" elements"<<std::endl;
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
      // solRef = _intSpace->evalFun(t, x);
      solRef = _solRef->eval(t, x);
      // printf("Solution at (%f,%f,%f) = %10.10f - ref = %10.10f\n", x[0], x[1], x[2], solInt, solRef);
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

      // std::cout<<"J = "<<J<<std::endl;

      area += J * w[k];
    }
  }
  norm = area;
}

void feNorm::computeIntegral(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh, feFunction *fun){
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

      area += fun->eval(t,x) * J * w[k];
    }
  }
  norm = area;
}