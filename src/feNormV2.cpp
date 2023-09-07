#include "feNormV2.h"
#include <cmath>


// Option available ;

//    - L2Norm_1Field
//    - L2Norm_MultipleFields


double feNormV2::compute(double tC, feSolutionContainerV2 *solutionContainer, feSolutionContainerV2 *referenceSolution)
{
  if(_normType == "L2Norm_1Field") {
     return computeL2ErrorNorm(tC, solutionContainer);
  
  } else if(_normType == "L2Norm_MultipleFields") {
     return computeL2ErrorNormVec(tC, solutionContainer);
  
  } else if(_normType == "fromMaxDC_1Field") {
      return computeL2ErrorNormFromMaxDC(tC, solutionContainer, referenceSolution);

  } else if(_normType == "fromMaxDC_MultipleFields") {
      return computeL2ErrorNormVecFromMaxDC(tC, solutionContainer, referenceSolution);

  }  else if(_normType == "IntMultipli") {
      return integrationLagrangeMupltiplicator(tC, solutionContainer);
  
  } else {
    feWarning("normType does't recognized for compute feNormV2");
  }
}

/* Compute the L2 norm of the error function e = uh - u at current solution time,
where uh is the discrete solution stored in sol and u is the solution given in referenceSolution. */
double feNormV2::computeL2ErrorNorm(double tC, feSolutionContainerV2 *solutionContainer)
{
  double L2Error = 0.0, solInt, solRef;
  std::vector<double> geoCoord;
  
  std::vector<double> x;
  
  std::vector<feInt> adr0;
  std::vector<double> sol0;
  std::vector<double> &solVec = solutionContainer->getSol(0);

  if(_referenceSolutions[0] == nullptr) feWarning("Reference solution is NULL.");

#if defined(HAVE_OMP)
#pragma omp parallel for private(geoCoord, adr0, sol0, solInt, solRef, x) reduction(+ : L2Error)schedule(dynamic)
#endif
  for(int iElm = 0; iElm < _nElm; ++iElm) {
    adr0.resize(_intSpaces[0]->getNbFunctions());
    _intSpaces[0]->initializeAddressingVector(_metaNumber->getNumbering(_intSpaces[0]->getFieldID()), iElm, adr0);

    sol0.resize(_intSpaces[0]->getNbFunctions());
    for(size_t i = 0; i<sol0.size(); ++i){
      sol0[i]=solVec[adr0[i]];
    }
    
    geoCoord.resize(3*_nNodePerElem);
    _mesh->getCoord(_intSpaces[0]->getCncGeoTag(), iElm, geoCoord);

    x.resize(3);
    for(int k = 0; k < _nQuad; ++k) {
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      // for (int i=0; i<3; i++){
      //   feInfo("%10.10f", x[i]);
      // }

      solInt = _intSpaces[0]->interpolateFieldAtQuadNode(sol0, k);
      solRef = (_referenceSolutions[0] != nullptr) ? _referenceSolutions[0]->eval(tC, x) : 0.0;

      L2Error +=(solInt - solRef) * (solInt - solRef) * _cnc->getJacobians()[_nQuad * iElm + k] * _wQuad[k];
      // feInfo("%10.10f", solInt - solRef);
    }
  }
  // feInfo("%10.10f",sqrt(L2Error));
  return sqrt(L2Error);
}



double feNormV2::computeL2ErrorNormVec(double tC, feSolutionContainerV2 *solutionContainer)
{
  double normL2 = 0.0, solIntU, solIntV, solRefU, solRefV;

  std::vector<double> geoCoord;
  std::vector<double> x;

  std::vector<feInt> adr0;
  std::vector<double> sol0;
  std::vector<feInt> adr1;
  std::vector<double> sol1;
  
  std::vector<double> &solVec = solutionContainer->getSol(0);

#if defined(HAVE_OMP)
#pragma omp parallel for private(geoCoord, adr0, adr1, sol0, sol1, solIntU, solIntV, solRefU, solRefV, x) reduction(+ : normL2) schedule(dynamic)
#endif
  for(int iElm = 0; iElm < _nElm; ++iElm) {
    adr0.resize(_intSpaces[0]->getNbFunctions());
    _intSpaces[0]->initializeAddressingVector(_metaNumber->getNumbering(_intSpaces[0]->getFieldID()), iElm, adr0);
    sol0.resize(_intSpaces[0]->getNbFunctions());
    for(size_t i = 0; i<sol0.size(); ++i){
      sol0[i]=solVec[adr0[i]];
    }

    adr1.resize(_intSpaces[1]->getNbFunctions()); 
    _intSpaces[1]->initializeAddressingVector(_metaNumber->getNumbering(_intSpaces[1]->getFieldID()), iElm, adr1);
    sol1.resize(_intSpaces[1]->getNbFunctions());
    for(size_t i = 0; i<sol1.size(); ++i){
      sol1[i]=solVec[adr1[i]];
    }

    geoCoord.resize(3*_nNodePerElem);
    _mesh->getCoord(_intSpaces[0]->getCncGeoTag(), iElm, geoCoord);
    
    x.resize(3);
    for(int k = 0; k < _nQuad; ++k) {
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      solIntU = _intSpaces[0]->interpolateFieldAtQuadNode(sol0, k);
      solIntV = _intSpaces[1]->interpolateFieldAtQuadNode(sol1, k);  
      
    
      solRefU = (_referenceSolutions[0] != nullptr) ? _referenceSolutions[0]->eval(tC, x) : 0.0;
      solRefV = (_referenceSolutions[0] != nullptr) ? _referenceSolutions[1]->eval(tC, x) : 0.0;

      // printf("%10.12f\n",solIntU - solRefU);
     
      normL2 += ((solIntU - solRefU) * (solIntU - solRefU) + (solIntV - solRefV) * (solIntV - solRefV)) * _cnc->getJacobians()[_nQuad * iElm + k] * _wQuad[k];
    }
  }
  // exit(-1);
  return sqrt(normL2);
}



double feNormV2::computeL2ErrorNormFromMaxDC(double tC, feSolutionContainerV2 *solutionContainer, feSolutionContainerV2 *referenceSolution)
{
  double L2Error = 0.0;
  double interpolatedCrtSol;
  double interpolatedRefSol;

  std::vector<double> geoCoord;
  
  std::vector<double> x;
  
  std::vector<feInt> adr;

  std::vector<double> crtSol;
  std::vector<double> refSol;

  std::vector<double> &crtSolVec = solutionContainer->getSol(0);
  std::vector<double> &refSolVec = referenceSolution->getSol(0);


#if defined(HAVE_OMP)
#pragma omp parallel for private(geoCoord, adr, crtSol, refSol, interpolatedCrtSol, interpolatedRefSol, x) reduction(+ : L2Error)schedule(dynamic)
#endif
  for(int iElm = 0; iElm < _nElm; ++iElm) 
  {
    adr.resize(_intSpaces[0]->getNbFunctions());

    _intSpaces[0]->initializeAddressingVector(_metaNumber->getNumbering(_intSpaces[0]->getFieldID()), iElm, adr);

    crtSol.resize(_intSpaces[0]->getNbFunctions());
    refSol.resize(_intSpaces[0]->getNbFunctions());
    for(size_t i = 0; i<crtSol.size(); ++i){
      crtSol[i] = crtSolVec[adr[i]];
      refSol[i] = refSolVec[adr[i]];
    }
    
    geoCoord.resize(3*_nNodePerElem);
    _mesh->getCoord(_intSpaces[0]->getCncGeoTag(), iElm, geoCoord);

    x.resize(3);
    for(int k = 0; k < _nQuad; ++k) {
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      interpolatedCrtSol = _intSpaces[0]->interpolateFieldAtQuadNode(crtSol, k);
      interpolatedRefSol = _intSpaces[0]->interpolateFieldAtQuadNode(refSol, k);

      L2Error +=(interpolatedCrtSol - interpolatedRefSol) * (interpolatedCrtSol - interpolatedRefSol) * _cnc->getJacobians()[_nQuad * iElm + k] * _wQuad[k];
    }
  }

  return sqrt(L2Error);
}



double feNormV2::computeL2ErrorNormVecFromMaxDC(double tC, feSolutionContainerV2 *solutionContainer, feSolutionContainerV2 *referenceSolution)
{
  double L2Error = 0.0;
  double interpolatedCrtSolU;
  double interpolatedRefSolU;

  double interpolatedCrtSolV;
  double interpolatedRefSolV;

  std::vector<double> geoCoord;
  std::vector<double> x;
  
  std::vector<feInt> adrU;
  std::vector<feInt> adrV;

  std::vector<double> crtSolU;
  std::vector<double> refSolU;
  std::vector<double> crtSolV;
  std::vector<double> refSolV;

  std::vector<double> &crtSolVec = solutionContainer->getSol(0);
  std::vector<double> &refSolVec = referenceSolution->getSol(0);


#if defined(HAVE_OMP)
#pragma omp parallel for private(geoCoord, adrU, adrV, crtSolU, refSolU, crtSolV, refSolV, interpolatedCrtSolU, interpolatedRefSolU, interpolatedCrtSolV, interpolatedRefSolV, x) reduction(+ : L2Error)schedule(dynamic)
#endif
  for(int iElm = 0; iElm < _nElm; ++iElm) 
  {
    adrU.resize(_intSpaces[0]->getNbFunctions());
    adrV.resize(_intSpaces[1]->getNbFunctions());

    _intSpaces[0]->initializeAddressingVector(_metaNumber->getNumbering(_intSpaces[0]->getFieldID()), iElm, adrU);
    _intSpaces[1]->initializeAddressingVector(_metaNumber->getNumbering(_intSpaces[1]->getFieldID()), iElm, adrV);

    crtSolU.resize(_intSpaces[0]->getNbFunctions());
    refSolU.resize(_intSpaces[0]->getNbFunctions());
    for(size_t i = 0; i<crtSolU.size(); ++i){
      crtSolU[i] = crtSolVec[adrU[i]];
      refSolU[i] = refSolVec[adrU[i]];
    }

    crtSolV.resize(_intSpaces[1]->getNbFunctions());
    refSolV.resize(_intSpaces[1]->getNbFunctions());
    for(size_t i = 0; i<crtSolV.size(); ++i){
      crtSolV[i] = crtSolVec[adrV[i]];
      refSolV[i] = refSolVec[adrV[i]];
    }
    

    geoCoord.resize(3*_nNodePerElem);
    _mesh->getCoord(_intSpaces[0]->getCncGeoTag(), iElm, geoCoord);

    x.resize(3);
    for(int k = 0; k < _nQuad; ++k) {
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      interpolatedCrtSolU = _intSpaces[0]->interpolateFieldAtQuadNode(crtSolU, k);
      interpolatedRefSolU = _intSpaces[0]->interpolateFieldAtQuadNode(refSolU, k);

      interpolatedCrtSolV = _intSpaces[1]->interpolateFieldAtQuadNode(crtSolV, k);
      interpolatedRefSolV = _intSpaces[1]->interpolateFieldAtQuadNode(refSolV, k);

      L2Error += (pow(interpolatedCrtSolU - interpolatedRefSolU, 2) + pow(interpolatedCrtSolV - interpolatedRefSolV, 2))* _cnc->getJacobians()[_nQuad * iElm + k] * _wQuad[k];
                 
    }
  }

  return sqrt(L2Error);
}



double feNormV2::integrationLagrangeMupltiplicator(double tC, feSolutionContainerV2 *solutionContainer)
{
  double intLamda = 0.0;
  double interpolatedCrtSol;

  std::vector<double> geoCoord;
  std::vector<double> x;
  
  std::vector<feInt> adr;

  std::vector<double> crtSol;

  std::vector<double> &crtSolVec = solutionContainer->getSol(0);

// #if defined(HAVE_OMP)
// #pragma omp parallel for private(geoCoord, adrU, adrV, crtSolU, refSolU, crtSolV, refSolV, interpolatedCrtSolU, interpolatedRefSolU, interpolatedCrtSolV, interpolatedRefSolV, x) reduction(+ : L2Error)schedule(dynamic)
// #endif
  for(int iElm = 0; iElm < _nElm; ++iElm) 
  {
    adr.resize(_intSpaces[0]->getNbFunctions());

    _intSpaces[0]->initializeAddressingVector(_metaNumber->getNumbering(_intSpaces[0]->getFieldID()), iElm, adr);

    crtSol.resize(_intSpaces[0]->getNbFunctions());
    for(size_t i = 0; i<crtSol.size(); ++i){
      crtSol[i] = crtSolVec[adr[i]];
    }
    
    geoCoord.resize(3*_nNodePerElem);
    _mesh->getCoord(_intSpaces[0]->getCncGeoTag(), iElm, geoCoord);

    x.resize(3);
    for(int k = 0; k < _nQuad; ++k) {
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      interpolatedCrtSol = _intSpaces[0]->interpolateFieldAtQuadNode(crtSol, k);

      intLamda += interpolatedCrtSol* _cnc->getJacobians()[_nQuad * iElm + k] * _wQuad[k];
      feInfo("Ccouco");
                 
    }
  }

  return intLamda;
}