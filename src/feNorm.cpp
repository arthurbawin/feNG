#include "feNorm.h"
#include <cmath>
#include "feQuadrature.h"

feNorm::feNorm(feSpace *intSpace, feMesh *mesh, int degQuad, feFunction *solRef)
  : _intSpace(intSpace), cncGeoTag(intSpace->getCncGeoTag()),
    geoSpace(mesh->getGeometricSpace(cncGeoTag)), nElmGeo(mesh->getNbElm(cncGeoTag)),
    dim(mesh->getDim()), nNodePerElem(intSpace->getNbNodePerElem()), _solRef(solRef),
    _degQuad(degQuad)
{
  feQuadrature *rule =
    new feQuadrature(_degQuad, intSpace->getDim(), intSpace->getCncGeo()->getForme());
  w = rule->getWeights();
  x = rule->getXPoints();
  _nQuad = rule->getNQuad();
  y = rule->getYPoints();
  z = rule->getZPoints();
  // Attention : ça changle la règle de quadrature pour les interpolants avant la résolution !
  // Il faudrait choisir idépendamment les deux sans conséquences
  geoSpace->setQuadratureRule(rule);
  intSpace->setQuadratureRule(rule);
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
  normL2 += (solInt - solRef) * (solInt - solRef);

  norm = sqrt(normL2);
  std::cout << "norme L2   " << norm << std::endl;
}

void feNorm::computeL2Norm(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh)
{
  double normL2 = 0.0, solInt, solRef, J, t = sol->getCurrentTime();
  int nElm = _intSpace->getNbElm();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _intSpace->initializeAddressingVector(metaNumber->getNumbering(_intSpace->getFieldID()), iElm);
    _intSpace->initializeSolution(sol);
    geoCoord = mesh->getCoord(cncGeoTag, iElm);

    for(int k = 0; k < _nQuad; ++k) {
      if(_intSpace->useGlobalFunctions()){
        solInt = _intSpace->interpolateSolutionAtQuadNode(iElm,k);
        // solInt = _intSpace->interpolateSolutionAtQuadNode_xDerivative(iElm,k);
      } else{
        solInt = _intSpace->interpolateSolutionAtQuadNode(k);  
      }

      std::vector<double> x(3, 0.0);
      geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      if(dim == 0) {
        J = 1.;
      } else if(dim == 1) {
        std::vector<double> j(3, 0.0);
        geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k,
                                                               j); // TODO : complete for dim > 1
        J = j[0];
      } else if(dim == 2) {
        std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
        std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
        geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
        geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
        J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];
      }
      // solRef = _intSpace->evalFun(t, x);
      solRef = (_solRef != nullptr) ? _solRef->eval(t, x) : 0.0;

      if(J < 0){
        printf("In feCncGeo::computeJacobians : Error - Element jacobian = %+-12.12e\n", J);
        exit(-1);
      }

      // printf("Solution at (%f,%f,%f) (ref = %f,%f,%f) = %10.10f - ref = %10.10f\n", x[0], x[1], x[2], _rQuad[k], _sQuad[k], 0.0, solInt, solRef);
      normL2 += (solInt - solRef) * (solInt - solRef) * J * w[k];
    }
  }
  norm = sqrt(normL2);
}

void feNorm::computeArea(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh) {
  double area = 0.0, solInt, J;
  int nElm = _intSpace->getNbElm();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _intSpace->initializeAddressingVector(metaNumber->getNumbering(_intSpace->getFieldID()), iElm);
    _intSpace->initializeSolution(sol);
    geoCoord = mesh->getCoord(cncGeoTag, iElm);
    for(int k = 0; k < _nQuad; ++k) {
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
                             feFunction *fun)
{
  double area = 0.0, solInt, J, t = sol->getCurrentTime();
  int nElm = _intSpace->getNbElm();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _intSpace->initializeAddressingVector(metaNumber->getNumbering(_intSpace->getFieldID()), iElm);
    _intSpace->initializeSolution(sol);
    geoCoord = mesh->getCoord(cncGeoTag, iElm);
    for(int k = 0; k < _nQuad; ++k) {
      solInt = _intSpace->interpolateSolutionAtQuadNode(k);
      // solInt = _intSpace->interpolateSolutionAtQuadNode(iElm,k);

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

      // area += fun->eval(t, x) * J * w[k];
      area += solInt * J * w[k];
    }
  }
  norm = area;
}

/* Estimates the L2 norm of the error taking an external solution as the reference solution. */
void feNorm::computeErrorNormFromExternalSolution(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh, 
  feMetaNumber *refMN, feSolution *refSol, feMesh *refMesh, const std::vector<feSpace *> refSpaces)
{
  double normL2 = 0.0, solInt, solRef;
  int nElm = _intSpace->getNbElm();

  if(fabs(sol->getCurrentTime() - refSol->getCurrentTime()) > 1e-3){
    printf("In feNorm::computeErrorNormFromExternalSolution : Warning - Solutions time differ by more than 1e-3 : tRef = %f - tSol = %f.\n", refSol->getCurrentTime(), sol->getCurrentTime());
    printf("In feNorm::computeErrorNormFromExternalSolution : Warning - Maybe the solution file does not match the current solution.\n");
  }

  std::vector<double> x(3, 0.0);
  std::vector<double> r(3, 0.0);
  std::vector<double> J = mesh->getCncGeoByTag(cncGeoTag)->getJacobians();

  // Identify the feSpace from among the set of spaces from the reference solution (we should use a map)
  // Aussi : la numérotation doit être recréée pour assigner correctement les ddl aux espaces.
  // Donc soit elle est la même que pour le problème courant, ce qui suppose que les fespaces sont les memes, ce qui a du sens, soit il faut redéfinir
  // un nouveau set de feSpace à partir d'une autre simulation, ce qui semble hasardeux. 
  bool matchingSpace = false;
  for(feSpace *fS : refSpaces){
    if(fS->getFieldID() == _intSpace->getFieldID() && fS->getCncGeoID() == _intSpace->getCncGeoID()){
      matchingSpace = true;

      feNumber *n = refMN->getNumbering(fS->getFieldID());

      for(int iElm = 0; iElm < nElm; ++iElm) {
        _intSpace->initializeAddressingVector(metaNumber->getNumbering(_intSpace->getFieldID()), iElm);
        _intSpace->initializeSolution(sol);
        geoCoord = mesh->getCoord(cncGeoTag, iElm);
        for(int k = 0; k < _nQuad; ++k) {
          geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
          // Current solution
          solInt = _intSpace->interpolateSolutionAtQuadNode(k);
          // Reference solution evaluated at quad node
          int elm = -1;
          bool isFound = static_cast<feMesh2DP1*>(refMesh)->locateVertex(x,elm,r);
          if(!isFound){
            printf("In feNorm::computeErrorNormFromExternalSolution : Warning - Point (%f, %f, %f) was not found in the mesh.\n", x[0], x[1], x[2]);
            solRef = solInt; // Points outside the mesh do not contribute
          } else{
            fS->initializeAddressingVector(n, elm);
            fS->initializeSolution(refSol);
            solRef = fS->interpolateSolution(r.data());
            normL2 += (solInt - solRef) * (solInt - solRef) * J[_nQuad * iElm + k] * w[k];
          }

          // normL2 += (solInt - solRef) * (solInt - solRef) * J[_nQuad * iElm + k] * w[k];
        }
      }

      norm = sqrt(normL2);

      break;
    }
  }

  if(!matchingSpace)
    printf("In feNorm::computeErrorNormFromExternalSolution : Error - No finite element space in the reference set matches the target space on this connectivity.\n");
}