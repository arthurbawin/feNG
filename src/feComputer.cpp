#include "feComputer.h"
#include <cmath>
#include "feQuadrature.h"

// Option available ;

//    - L2Norm_1Field
//    - L2Norm_MultipleFields
//    - H1Norm_MultipleFields
//    - SemiH1Norm_MultipleFields
//    - Int
//    - IntNormalx
//    - IntNormaly
//    - NormLambda
//    - NormLambdaNS

void feComputer::Compute(feSolution *sol, int i)
{
  if(_computerType == "L2Norm_1Field") {
    _Result[i] = computeL2ErrorNorm(sol);
  } else if(_computerType == "L2Norm_MultipleFields") {
    _Result[i] = computeL2ErrorNormVec(sol);
  } else if(_computerType == "H1Norm_MultipleFields") {
    _Result[i] = computeH1ErrorNormVec(sol);
  } else if(_computerType == "SemiH1Norm_MultipleFields") {
    _Result[i] = computeSemiH1ErrorNormVec(sol);
  } else if(_computerType == "NormL2_LagrangeMultiplier") {
    _Result[i] = computeNormL2Lambda(sol);
  } else if(_computerType == "Int_LagrangeMultiplier") {
    _Result[i] = computeIntLagrangeMultiplier(sol);  
  }else if(_computerType == "Int_MultiplierNS") {
    _Result[i] = computeIntMultiplierNS(sol);  
  } else {
  }
}

/* Compute the L2 norm of the error function e = uh - u at current solution time,
where uh is the discrete solution stored in sol and u is the solution given in referenceSolution. */
double feComputer::computeL2ErrorNorm(feSolution *sol)
{
  double L2Error = 0.0, ERROR = 0.0, solInt, solRef, t = sol->getCurrentTime();
  std::vector<double> geoCoord;
  std::vector<double> &w = _intSpace->getQuadratureWeights();
  int nElm = _intSpace->getNbElm();
  int nQuad = _geoSpace->getNbQuadPoints();
  std::vector<double> x;
  
  std::vector<feInt> adr0;
  std::vector<double> sol0;
  std::vector<double> &solVec= sol->getSolutionReference();

  if(_referenceSolution == nullptr) feWarning("Reference solution is NULL.");

#if defined(HAVE_OMP)
#pragma omp parallel for private(geoCoord, adr0, sol0, solInt, solRef, x) reduction(+ : ERROR) reduction(+ : L2Error)schedule(dynamic)
#endif
  for(int iElm = 0; iElm < nElm; ++iElm) {
    adr0.resize(_intSpace->getNbFunctions());
    _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()), iElm, adr0);
    sol0.resize(_intSpace->getNbFunctions());
    for(size_t i = 0; i<sol0.size(); ++i){
      sol0[i]=solVec[adr0[i]];
    }
    geoCoord.resize(3*_nNodePerElem);
    _mesh->getCoord(_intSpace->getCncGeoTag(), iElm, geoCoord);
    ERROR = 0.0;
    x.resize(3);
    for(int k = 0; k < nQuad; ++k) {
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      solInt = _intSpace->interpolateFieldAtQuadNode(sol0, k);
      solRef = (_referenceSolution != nullptr) ? _referenceSolution->eval(t, x) : 0.0;
      L2Error +=
        (solInt - solRef) * (solInt - solRef) * _cnc->getJacobians()[nQuad * iElm + k] * w[k];
      ERROR+=(solInt - solRef) * _cnc->getJacobians()[nQuad * iElm + k] * w[k];
    }
    // printf("%24.18e ", ERROR);
    // std::cout<<",";
  }
  return sqrt(L2Error);
}

double feComputer::computeL2ErrorNormVec(feSolution *sol)
{
  double normL2 = 0.0, solIntU, solIntV, t = sol->getCurrentTime();
  int nElm = _VecfeSpace[0]->getNbElm();
  int nQuad = _geoSpace->getNbQuadPoints();
  std::vector<double> &w = _VecfeSpace[0]->getQuadratureWeights();
  std::vector<double> geoCoord;
  std::vector<double> x;
  std::vector<double> solRef;

  std::vector<feInt> adr0;
  std::vector<double> sol0;
  std::vector<feInt> adr1;
  std::vector<double> sol1;
  std::vector<double> &solVec = sol->getSolutionReference();

#if defined(HAVE_OMP)
#pragma omp parallel for private(geoCoord, adr0, adr1, sol0, sol1, solIntU, solIntV, solRef, x) reduction(+ : normL2) schedule(dynamic)
#endif
  for(int iElm = 0; iElm < nElm; ++iElm) {
    adr0.resize(_VecfeSpace[0]->getNbFunctions());
    _VecfeSpace[0]->initializeAddressingVector(
      _metaNumber->getNumbering(_VecfeSpace[0]->getFieldID()), iElm, adr0);
    sol0.resize(_VecfeSpace[0]->getNbFunctions());
    for(size_t i = 0; i<sol0.size(); ++i){
      sol0[i]=solVec[adr0[i]];
    }

    adr1.resize(_VecfeSpace[1]->getNbFunctions()); 
    _VecfeSpace[1]->initializeAddressingVector(
      _metaNumber->getNumbering(_VecfeSpace[1]->getFieldID()), iElm, adr1);
    sol1.resize(_VecfeSpace[1]->getNbFunctions());
    for(size_t i = 0; i<sol1.size(); ++i){
      sol1[i]=solVec[adr1[i]];
    }

    geoCoord.resize(3*_nNodePerElem);
    _mesh->getCoord(_VecfeSpace[0]->getCncGeoTag(), iElm, geoCoord);
    x.resize(3);
    for(int k = 0; k < nQuad; ++k) {
      solIntU = _VecfeSpace[0]->interpolateFieldAtQuadNode(sol0, k);
      solIntV = _VecfeSpace[1]->interpolateFieldAtQuadNode(sol1, k);  
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      solRef.resize(6);
      if(_referenceVectorSolution != nullptr) {
        _referenceVectorSolution->eval(t, x, solRef);
      }
      // std::cout<< "========"<<std::endl;
      // std::cout<<"    SolU   "<<solIntU<<std::endl;
      // std::cout<<"    SolRefU   "<<solRef[0]<<std::endl;
      normL2 += ((solIntU - solRef[0]) * (solIntU - solRef[0]) +
                 (solIntV - solRef[1]) * (solIntV - solRef[1])) *
                _cnc->getJacobians()[nQuad * iElm + k] * w[k];
    }
  }
  return sqrt(normL2);
}

double feComputer::computeH1ErrorNormVec(feSolution *sol)
{
  double normH1 = 0.0, solIntU, solIntV, J, t = sol->getCurrentTime(), E11, E12, E22;
  int nElm = _VecfeSpace[0]->getNbElm();
  int nQuad = _geoSpace->getNbQuadPoints();
  std::vector<double> &w = _VecfeSpace[0]->getQuadratureWeights();
  std::vector<double> geoCoord(3*_nNodePerElem);
  std::vector<feInt> adr0(_VecfeSpace[0]->getNbFunctions());
  std::vector<feInt> adr1(_VecfeSpace[1]->getNbFunctions());
  std::vector<double> sol0(adr0.size());
  std::vector<double> sol1(adr1.size());
  std::vector<double> &solVec = sol->getSolutionReference();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _VecfeSpace[0]->initializeAddressingVector(
      _metaNumber->getNumbering(_VecfeSpace[0]->getFieldID()), iElm, adr0);
    for(size_t i = 0; i<sol0.size(); ++i){
      sol0[i]=solVec[adr0[i]];
    }
    
    _VecfeSpace[1]->initializeAddressingVector(
      _metaNumber->getNumbering(_VecfeSpace[1]->getFieldID()), iElm, adr1);
    for(size_t i = 0; i<sol1.size(); ++i){
      sol1[i]=solVec[adr1[i]];
    }
    
    _mesh->getCoord(_VecfeSpace[0]->getCncGeoTag(), iElm, geoCoord);
    for(int k = 0; k < nQuad; ++k) {
      solIntU = _VecfeSpace[0]->interpolateFieldAtQuadNode(sol0, k);
      solIntV = _VecfeSpace[1]->interpolateFieldAtQuadNode(sol1, k);

      std::vector<double> x(3, 0.0);
      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      J = _cnc->getJacobians()[nQuad * iElm + k];
      double drdx = dxds[1] / J;
      double drdy = -dxds[0] / J;
      double dsdx = -dxdr[1] / J;
      double dsdy = dxdr[0] / J;
      double dudx, dudy, dvdx, dvdy;
      dudx = _VecfeSpace[0]->interpolateFieldAtQuadNode_rDerivative(sol0, k) * drdx +
             _VecfeSpace[0]->interpolateFieldAtQuadNode_sDerivative(sol0, k) * dsdx;
      dudy = _VecfeSpace[0]->interpolateFieldAtQuadNode_rDerivative(sol0, k) * drdy +
             _VecfeSpace[0]->interpolateFieldAtQuadNode_sDerivative(sol0, k) * dsdy;
      dvdx = _VecfeSpace[1]->interpolateFieldAtQuadNode_rDerivative(sol1, k) * drdx +
             _VecfeSpace[1]->interpolateFieldAtQuadNode_sDerivative(sol1, k) * dsdx;
      dvdy = _VecfeSpace[1]->interpolateFieldAtQuadNode_rDerivative(sol1, k) * drdy +
             _VecfeSpace[1]->interpolateFieldAtQuadNode_sDerivative(sol1, k) * dsdy;
      std::vector<double> solRef(6, 0);
      if(_referenceVectorSolution != nullptr) {
        _referenceVectorSolution->eval(t, x, solRef);
      }
      E11 = 2 * dudx - 2 * solRef[2]; // solRef[2] = dUref/dx
      E22 = 2 * dvdy - 2 * solRef[5]; // solRef[2] = dVref/dy
      E12 = (dudy + dvdx) - (solRef[3] + solRef[4]);
      normH1 += ((solIntU - solRef[0]) * (solIntU - solRef[0]) + (solIntV - solRef[1]) * (solIntV - solRef[1]) + (E11 * E11 + 2 * E12 * E12 + E22 * E22)) * J * w[k];
    }
  }
  return sqrt(normH1);
}

double feComputer::computeSemiH1ErrorNormVec(feSolution *sol)
{
  double semiNormH1 = 0.0, solIntU, solIntV, J, t = sol->getCurrentTime(), E11, E12, E22;
  int nElm = _VecfeSpace[0]->getNbElm();
  int nQuad = _geoSpace->getNbQuadPoints();
  std::vector<double> &w = _VecfeSpace[0]->getQuadratureWeights();
  std::vector<feInt> adr0(_VecfeSpace[0]->getNbFunctions());
  std::vector<feInt> adr1(_VecfeSpace[1]->getNbFunctions());
  std::vector<double> sol0(adr0.size());
  std::vector<double> sol1(adr1.size());
  std::vector<double> &solVec = sol->getSolutionReference();
  std::vector<double> geoCoord(3*_nNodePerElem);

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _VecfeSpace[0]->initializeAddressingVector(
      _metaNumber->getNumbering(_VecfeSpace[0]->getFieldID()), iElm, adr0);
    for(size_t i = 0; i<sol0.size(); ++i){
      sol0[i]=solVec[adr0[i]];
    }
    
    _VecfeSpace[1]->initializeAddressingVector(
      _metaNumber->getNumbering(_VecfeSpace[1]->getFieldID()), iElm, adr1);
    for(size_t i = 0; i<sol1.size(); ++i){
      sol1[i]=solVec[adr1[i]];
    }
    
    _mesh->getCoord(_VecfeSpace[0]->getCncGeoTag(), iElm, geoCoord);
    for(int k = 0; k < nQuad; ++k) {
      solIntU = _VecfeSpace[0]->interpolateFieldAtQuadNode(sol0, k);
      solIntV = _VecfeSpace[1]->interpolateFieldAtQuadNode(sol1, k);

      std::vector<double> x(3, 0.0);
      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      J = _cnc->getJacobians()[nQuad * iElm + k];
      double drdx = dxds[1] / J;
      double drdy = -dxds[0] / J;
      double dsdx = -dxdr[1] / J;
      double dsdy = dxdr[0] / J;
      double dudx, dudy, dvdx, dvdy;
      dudx = _VecfeSpace[0]->interpolateFieldAtQuadNode_rDerivative(sol0, k) * drdx +
             _VecfeSpace[0]->interpolateFieldAtQuadNode_sDerivative(sol0, k) * dsdx;
      dudy = _VecfeSpace[0]->interpolateFieldAtQuadNode_rDerivative(sol0, k) * drdy +
             _VecfeSpace[0]->interpolateFieldAtQuadNode_sDerivative(sol0, k) * dsdy;
      dvdx = _VecfeSpace[1]->interpolateFieldAtQuadNode_rDerivative(sol1, k) * drdx +
             _VecfeSpace[1]->interpolateFieldAtQuadNode_sDerivative(sol1, k) * dsdx;
      dvdy = _VecfeSpace[1]->interpolateFieldAtQuadNode_rDerivative(sol1, k) * drdy +
             _VecfeSpace[1]->interpolateFieldAtQuadNode_sDerivative(sol1, k) * dsdy;
      std::vector<double> solRef(6, 0);
      if(_referenceVectorSolution != nullptr) {
        _referenceVectorSolution->eval(t, x, solRef);
      }
      E11 = 2 * dudx - 2 * solRef[2]; // solRef[2] = dUref/dx
      E22 = 2 * dvdy - 2 * solRef[5]; // solRef[2] = dVref/dy
      E12 = (dudy + dvdx) - (solRef[3] + solRef[4]);
      semiNormH1 += (E11 * E11 + 2 * E12 * E12 + E22 * E22) * J * w[k];
    }
  }
  return sqrt(semiNormH1);
}

double feComputer::computeNormL2Lambda(feSolution *sol)
{
  double normL2Lambda = 0.0, solInt, J, t = sol->getCurrentTime();
  std::vector<double> geoCoord(3*_nNodePerElem);
  std::vector<double> &w = _intSpace->getQuadratureWeights();
  int nElm = _intSpace->getNbElm();
  int nQuad = _geoSpace->getNbQuadPoints();
  std::vector<feInt> adr0(_intSpace->getNbFunctions());
  std::vector<double> sol0(adr0.size());
  std::vector<double> &solVec = sol->getSolutionReference();

  if(_referenceVectorSolution == nullptr) feWarning("Reference solution is NULL.");

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()), iElm, adr0);
    for(size_t i = 0; i<sol0.size(); ++i){
      sol0[i]=solVec[adr0[i]];
    }
    _mesh->getCoord(_intSpace->getCncGeoTag(), iElm, geoCoord);
    for(int k = 0; k < nQuad; ++k) {
      std::vector<double> x(3, 0.0);
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      solInt = _intSpace->interpolateFieldAtQuadNode(sol0, k);
      std::vector<double> gradf(3, 0);
      double Kd;
      if(_referenceVectorSolution != nullptr) {
        _referenceVectorSolution->eval(t, x, gradf);
        Kd = _referenceVectorSolution->getParam()[0];
      }
      J = _cnc->getJacobians()[nQuad * iElm + k];
      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      double Nx = dxdr[1];
      double Ny = -dxdr[0];
      double N = sqrt(Nx * Nx + Ny * Ny);
      Nx /= N;
      Ny /= N;

      normL2Lambda += (solInt + Kd * (gradf[0] * Nx + gradf[1] * Ny)) *
                      (solInt + Kd * (gradf[0] * Nx + gradf[1] * Ny)) * J * w[k];
    }
  }
  return sqrt(normL2Lambda);
}

double feComputer::computeIntLagrangeMultiplier(feSolution *sol)
{
  double IntLambda = 0.0, solInt, J, t = sol->getCurrentTime();
  std::vector<double> geoCoord(3*_nNodePerElem);
  std::vector<double> &w = _intSpace->getQuadratureWeights();
  int nElm = _intSpace->getNbElm();
  int nQuad = _geoSpace->getNbQuadPoints();
  std::vector<feInt> adr0(_intSpace->getNbFunctions());
  std::vector<double> sol0(adr0.size());
  std::vector<double> &solVec = sol->getSolutionReference();

  if(_referenceVectorSolution == nullptr) feWarning("Reference solution is NULL.");

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()), iElm, adr0);
    for(size_t i = 0; i<sol0.size(); ++i){
      sol0[i]=solVec[adr0[i]];
    }
    _mesh->getCoord(_intSpace->getCncGeoTag(), iElm, geoCoord);
    for(int k = 0; k < nQuad; ++k) {
      std::vector<double> x(3, 0.0);
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      solInt = _intSpace->interpolateFieldAtQuadNode(sol0, k);
      std::vector<double> gradf(3, 0);
      double Kd;
      if(_referenceVectorSolution != nullptr) {
        _referenceVectorSolution->eval(t, x, gradf);
        Kd = _referenceVectorSolution->getParam()[0];
      }
      J = _cnc->getJacobians()[nQuad * iElm + k];
      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      double Nx = dxdr[1];
      double Ny = -dxdr[0];
      double N = sqrt(Nx * Nx + Ny * Ny);
      Nx /= N;
      Ny /= N;
      IntLambda += (solInt + Kd * (gradf[0] * Nx + gradf[1] * Ny)) * J * w[k];
    }
  }
  return IntLambda;
}

double feComputer::computeIntMultiplierNS(feSolution *sol)
{
  double IntLambda = 0.0, solInt, p, J, t = sol->getCurrentTime(), T11, T12, T22;
  std::vector<double> geoCoord(3*_nNodePerElem);
  std::vector<double> &w = _intSpace->getQuadratureWeights();
  int nElm = _intSpace->getNbElm();
  int nQuad = _geoSpace->getNbQuadPoints();
  std::vector<feInt> adr0(_intSpace->getNbFunctions());
  std::vector<double> sol0(adr0.size());
  std::vector<double> &solVec = sol->getSolutionReference();

  if(_referenceVectorSolution == nullptr) feWarning("Reference solution is NULL.");

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _intSpace->initializeAddressingVector(
      _metaNumber->getNumbering(_intSpace->getFieldID()), iElm, adr0);
    for(size_t i = 0; i<sol0.size(); ++i){
      sol0[i]=solVec[adr0[i]];
    }
    _mesh->getCoord(_intSpace->getCncGeoTag(), iElm, geoCoord);
    for(int k = 0; k < nQuad; ++k) {
      solInt = _intSpace->interpolateFieldAtQuadNode(sol0, k);

      std::vector<double> x(3, 0.0);
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      J = _cnc->getJacobians()[nQuad * iElm + k];
      p = (_referenceSolution != nullptr) ? _referenceSolution->eval(t, x) : 0.0; // here is the pressure
      double mu;
      std::vector<double> gradf(3, 0);
      if(_referenceVectorSolution != nullptr) {
        _referenceVectorSolution->eval(t, x,
                         gradf); // gradf[0] = dudx gradf[1] = dudy gradf[2] = dvdx gradf[3] = dvdy
        mu = _referenceVectorSolution->getParam()[1];
      }
      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      double Nx = dxdr[1];
      double Ny = -dxdr[0];
      double N = sqrt(Nx * Nx + Ny * Ny);
      Nx /= N;
      Ny /= N;

      T11 = 2 * gradf[0];
      T12 = gradf[1] + gradf[2];
      T22 = 2 * gradf[3];
      // std::cout<<"solInt   "<<solInt<<std::endl;
      // std::cout<<"J   "<<J<<std::endl;
      // std::cout<<"IntForce   "<<- p * Nx + mu * (T11 * Nx + T12 * Ny)<<std::endl;
      // std::cout<<"IntForce   "<< mu * (T11 * Nx + T12 * Ny)<<std::endl;

      // Integral +=  (solInt + mu * (T11*Nx + T12*Ny)) * J * w[k];
      // std::cout<<"ErrInt   "<<(-solInt - p * Nx + mu * (T11 * Nx + T12 * Ny))<<std::endl;
      IntLambda += (solInt - p * Nx + mu * (T11 * Nx + T12 * Ny)) * J * w[k];
      // std::cout<<"IntLambda   "<<IntLambda<<std::endl;
    }
  }
  return IntLambda;
}

// double feComputer::computeIntMultiplierNS(feSolution *sol)
// {
//   double IntLambda = 0.0, solIntU, solIntV, p, J, t = sol->getCurrentTime(), T11, T12, T22;
//   std::vector<double> geoCoord;
//   std::vector<double> &w = _intSpace->getQuadratureWeights();
//   int nElm = _VecfeSpace[0]->getNbElm();
//   int nQuad = _geoSpace->getNbQuadPoints();

//   if(_referenceVectorSolution == nullptr) feWarning("Reference solution is NULL.");

//   for(int iElm = 0; iElm < nElm; ++iElm) {
//     _VecfeSpace[0]->initializeAddressingVector(
//       _metaNumber->getNumbering(_VecfeSpace[0]->getFieldID()), iElm);
//     _VecfeSpace[0]->initializeSolution(sol);
//     _VecfeSpace[1]->initializeAddressingVector(
//       _metaNumber->getNumbering(_VecfeSpace[1]->getFieldID()), iElm);
//     _VecfeSpace[1]->initializeSolution(sol);
//     geoCoord = _mesh->getCoord(_intSpace->getCncGeoTag(), iElm);
//     for(int k = 0; k < nQuad; ++k) {
//       solIntU = _VecfeSpace[0]->interpolateSolutionAtQuadNode(k);
//       solIntV = _VecfeSpace[1]->interpolateSolutionAtQuadNode(k);

//       std::vector<double> x(3, 0.0);
//       _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
//       J = _cnc->getJacobians()[nQuad * iElm + k];
//       p = (_referenceSolution != nullptr) ? _referenceSolution->eval(t, x) : 0.0; // here is the pressure
//       double mu;
//       std::vector<double> gradf(3, 0);
//       if(_referenceVectorSolution != nullptr) {
//         _referenceVectorSolution->eval(t, x,
//                          gradf); // gradf[0] = dudx gradf[1] = dudy gradf[2] = dvdx gradf[3] = dvdy
//         mu = _referenceVectorSolution->getParam()[0];
//       }
//       std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
//       std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
//       _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
//       _geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
//       double Nx = dxdr[1];
//       double Ny = -dxdr[0];
//       double N = sqrt(Nx * Nx + Ny * Ny);
//       Nx /= N;
//       Ny /= N;

//       T11 = 2 * gradf[0];
//       T12 = gradf[1] + gradf[2];
//       T22 = 2 * gradf[3];
//       // std::cout<<"solInt   "<<solInt<<std::endl;
//       // std::cout<<"J   "<<J<<std::endl;
//       // std::cout<<"IntForce   "<<- p * Nx + mu * (T11 * Nx + T12 * Ny)<<std::endl;

//       // Integral +=  (solInt + mu * (T11*Nx + T12*Ny)) * J * w[k];
//       // std::cout<<"ErrInt   "<<(-solInt - p * Nx + mu * (T11 * Nx + T12 * Ny))<<std::endl;
//       IntLambda += (-solInt - p * Nx + mu * (T11 * Nx + T12 * Ny)) * J * w[k];
//       std::cout<<"IntLambda   "<<IntLambda<<std::endl;
//     }
//   }
//   return IntLambda;
// }
