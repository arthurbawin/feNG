#include "feNorm.h"
#include <cmath>
#include "feQuadrature.h"

feNorm::feNorm(feSpace *intSpace, feMesh *mesh, int degQuad, feFunction *solRef, feVectorFunction *gradRef, std::string TypeNorm)
  : _intSpace(intSpace), cncGeoTag(intSpace->getCncGeoTag()),
    geoSpace(mesh->getGeometricSpace(cncGeoTag)), nElmGeo(mesh->getNbElm(cncGeoTag)),
    dim(intSpace->getDim()), nNodePerElem(intSpace->getNbNodePerElem()), _solRef(solRef),
    _degQuad(degQuad), _TypeNorm(TypeNorm) {
  feQuadrature *rule =
    new feQuadrature(_degQuad, intSpace->getDim(), intSpace->getCncGeo()->getForme());
  w = rule->getWeights();
  x = rule->getXPoints();
  _nQuad = rule->getNQuad();
  y = rule->getYPoints();
  z = rule->getZPoints();
  if(gradRef != nullptr){
    _gradRef = gradRef;
  }
  // Attention : ça changle la règle de quadrature pour les interpolants avant la résolution !
  // Il faudrait choisir idépendamment les deux sans conséquences
  intSpace->setQuadratureRule(rule);
  geoSpace->setQuadratureRule(rule);
  delete rule;
}

void feNorm::computeL2Norm0D(feSolution *sol) {
  double normL2 = 0.0, solInt, solRef, t = sol->getCurrentTime();

  // solRef = _intSpace->evalFun(t, x);
  solRef = (_solRef != nullptr) ? _solRef->eval(t, x) : 0.0;
  std::cout << "la solution de ref est " << solRef << " au temps " << t << std::endl;
  solInt = sol->getSolAtDOF(0);
  std::cout << "la solution de numérique est " << solInt << " au temps " << t << std::endl;
  // printf("Solution at (%f,%f,%f) = %10.10f - ref = %10.10f\n", x[0], x[1], x[2], solInt,
  // solRef);
  normL2 += (abs(solInt) - abs(solRef)) * (abs(solInt) - abs(solRef));

  norm = sqrt(normL2);
}

void feNorm::computeL2Norm(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh) {
  double normL2 = 0.0, solInt, solRef, J, t = sol->getCurrentTime();
  int nElm = _intSpace->getNbElm();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _intSpace->initializeAddressingVector(metaNumber->getNumbering(_intSpace->getFieldID()), iElm);
    _intSpace->initializeSolution(sol);
    geoCoord = mesh->getCoord(cncGeoTag, iElm);
    for(int k = 0; k < _nQuad; ++k) {
      solInt = _intSpace->interpolateSolutionAtQuadNode(k);

      std::vector<double> x(3, 0.0);
      geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      if(dim == 0) {
        J = 1.;
      } 
      else if(dim == 1) {
          std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
          geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
          J = sqrt(dxdr[0] * dxdr[0]+dxdr[1] * dxdr[1]);
        } 
      else if(dim == 2) {
        std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
        std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
        geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
        geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
        J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];
      } else {
      }
      solRef = (_solRef != nullptr) ? _solRef->eval(t, x) : 0.0;  
      normL2 += (solInt - solRef) * (solInt - solRef) * J * w[k];
    }
  }
  norm = sqrt(normL2);
}

void feNorm::computeNormLambda(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh) {
  double normL2 = 0.0, solInt, solRef, J, t = sol->getCurrentTime();
  int nElm = _intSpace->getNbElm();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _intSpace->initializeAddressingVector(metaNumber->getNumbering(_intSpace->getFieldID()), iElm);
    _intSpace->initializeSolution(sol);
    geoCoord = mesh->getCoord(cncGeoTag, iElm);
    for(int k = 0; k < _nQuad; ++k) {
      solInt = _intSpace->interpolateSolutionAtQuadNode(k);

      std::vector<double> x(3, 0.0);
      geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      if(dim == 0) {
        J = 1.;
      } 
      else if(dim == 1) {
          std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
          geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
          J = sqrt(dxdr[0] * dxdr[0]+dxdr[1] * dxdr[1]);
        } 
      else {
      }

      solRef = (_solRef != nullptr) ? _solRef->eval(t, x) : 0.0;
      double Kd;
      std::vector<double> gradf(3, 0);
      if (_gradRef != nullptr) {
        _gradRef->eval(t, x, gradf);
        Kd = _gradRef->getParam()[0];

      }
      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      double Nx = dxdr[1] ;
      double Ny = -dxdr[0];
      double N= sqrt(Nx*Nx + Ny*Ny);
      Nx/=N;
      Ny/=N;

      if(_gradRef != nullptr) {
        if(_TypeNorm == "IntFlux"){
          normL2 += (solInt + Kd*(gradf[0]*Nx +  gradf[1]*Ny))* J* w[k];
        }
        else if(_TypeNorm == "NormL2"){
          normL2 += (solInt + Kd*(gradf[0]*Nx +  gradf[1]*Ny )) * (solInt + Kd*(gradf[0]*Nx +  gradf[1]*Ny)) *J * w[k];
        }
      }
    }
  }
  if(_TypeNorm == "NormL2"){
    norm = sqrt(normL2);
  }
  else if(_TypeNorm == "IntFlux"){
    norm = normL2;
  }

}

void feNorm::computeArea(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh) {
  double area = 0.0, solInt, J;
  int nElm = _intSpace->getNbElm();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _intSpace->initializeAddressingVector(metaNumber->getNumbering(_intSpace->getFieldID()), iElm);
    _intSpace->initializeSolution(sol);
    geoCoord = mesh->getCoord(cncGeoTag, iElm);
    for(int k = 0; k < _nQuad; ++k) {
      solInt = _intSpace->interpolateSolutionAtQuadNode(k);

      std::vector<double> x(3, 0.0);
      geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      // J = geoSpace->interpolateFieldAtQuadNode_rDerivative(xGeoCoord, k);
      // std::vector<double> j(3, 0.0);
      // geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, j); // TODO : complete
      // for dim > 1 J = j[0];
      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
      J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];

      area += J * w[k];
    }
  }
  norm = area;
}

void feNorm::computeIntegral(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh,
                             feFunction *fun) {
  double area = 0.0, solInt, J, t = sol->getCurrentTime();
  int nElm = _intSpace->getNbElm();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _intSpace->initializeAddressingVector(metaNumber->getNumbering(_intSpace->getFieldID()), iElm);
    _intSpace->initializeSolution(sol);
    geoCoord = mesh->getCoord(cncGeoTag, iElm);
    for(int k = 0; k < _nQuad; ++k) {
      solInt = _intSpace->interpolateSolutionAtQuadNode(k);

      std::vector<double> x(3, 0.0);
      geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      // J = geoSpace->interpolateFieldAtQuadNode_rDerivative(xGeoCoord, k);
      // std::vector<double> j(3, 0.0);
      // geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, j); // TODO : complete
      // for dim > 1 J = j[0];
      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
      J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];

      area += fun->eval(t, x) * J * w[k];
    }
  }
  norm = area;
}