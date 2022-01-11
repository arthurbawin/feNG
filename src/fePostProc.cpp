#include "fePostProc.h"
#include <cmath>
#include "feQuadrature.h"

/* Compute the measure (length, area, volume) of the associated geometric connectivity */
double fePostProc::computeMeasure()
{
  double meas = 0.0;
  std::vector<double> geoCoord;
  std::vector<double> &w = _intSpace->getQuadratureWeights();
  int nElm = _intSpace->getNbElm();
  int nQuad = _geoSpace->getNbQuadPoints();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    geoCoord = _mesh->getCoord(_intSpace->getCncGeoTag(), iElm);
    for(int k = 0; k < nQuad; ++k) {
      meas += _cnc->getJacobians()[nQuad * iElm + k] * w[k];
    }
  }
  return meas;
}

/* Compute the integral of the discrete solution sol on the associated geometric connectivity */
double fePostProc::computeSolutionIntegral(feSolution *sol)
{
  double integral = 0.0, solInt;
  std::vector<double> geoCoord;
  std::vector<double> &w = _intSpace->getQuadratureWeights();
  int nElm = _intSpace->getNbElm();
  int nQuad = _geoSpace->getNbQuadPoints();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()), iElm);
    _intSpace->initializeSolution(sol);
    geoCoord = _mesh->getCoord(_intSpace->getCncGeoTag(), iElm);
    for(int k = 0; k < nQuad; ++k) {
      solInt = _intSpace->interpolateSolutionAtQuadNode(k);
      integral += solInt * _cnc->getJacobians()[nQuad * iElm + k] * w[k];
    }
  }
  return integral;
}

/* Compute the integral at time t of a given function fun on the associated geometric connectivity
 */
double fePostProc::computeFunctionIntegral(feFunction *fun, double t)
{
  double integral = 0.0;
  std::vector<double> geoCoord;
  std::vector<double> &w = _intSpace->getQuadratureWeights();
  int nElm = _intSpace->getNbElm();
  int nQuad = _geoSpace->getNbQuadPoints();

  for(int iElm = 0; iElm < nElm; ++iElm) {
    geoCoord = _mesh->getCoord(_intSpace->getCncGeoTag(), iElm);
    for(int k = 0; k < nQuad; ++k) {
      std::vector<double> x(3, 0.0);
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      integral += fun->eval(t, x) * _cnc->getJacobians()[nQuad * iElm + k] * w[k];
    }
  }
  return integral;
}

/* Compute the L2 norm of the error function e = uh - u at current solution time,
where uh is the discrete solution stored in sol and u is the solution given in referenceSolution. */
double fePostProc::computeL2ErrorNorm(feSolution *sol)
{
  double L2Error = 0.0, solInt, solRef, t = sol->getCurrentTime();
  std::vector<double> geoCoord;
  std::vector<double> &w = _intSpace->getQuadratureWeights();
  int nElm = _intSpace->getNbElm();
  int nQuad = _geoSpace->getNbQuadPoints();

  if(_referenceSolution == nullptr) feWarning("Reference solution is NULL.");

  for(int iElm = 0; iElm < nElm; ++iElm) {
    _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()), iElm);
    _intSpace->initializeSolution(sol);
    geoCoord = _mesh->getCoord(_intSpace->getCncGeoTag(), iElm);
    for(int k = 0; k < nQuad; ++k) {
      std::vector<double> x(3, 0.0);
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
      solInt = _intSpace->interpolateSolutionAtQuadNode(k);
      solRef = (_referenceSolution != nullptr) ? _referenceSolution->eval(t, x) : 0.0;
      L2Error +=
        (solInt - solRef) * (solInt - solRef) * _cnc->getJacobians()[nQuad * iElm + k] * w[k];
    }
  }
  return sqrt(L2Error);
}

// void feNorm::computeL2Norm0D(feSolution *sol) {
//   double normL2 = 0.0, solInt, solRef, t = sol->getCurrentTime();

//   // solRef = _intSpace->evalFun(t, x);
//   solRef = (_solRef != nullptr) ? _solRef->eval(t, x) : 0.0;
//   std::cout << "la solution de ref est " << solRef << " au temps " << t << std::endl;
//   solInt = sol->getSolAtDOF(0);
//   std::cout << "la solution de numérique est " << solInt << " au temps " << t << std::endl;
//   // printf("Solution at (%f,%f,%f) = %10.10f - ref = %10.10f\n", x[0], x[1], x[2], solInt,
//   // solRef);
//   normL2 += (solInt - solRef) * (solInt - solRef);

//   norm = sqrt(normL2);
//   std::cout << "norme L2   " << norm << std::endl;
// }

// /* Estimates the L2 norm of the error taking an external solution as the reference solution. */
// void feNorm::computeErrorNormFromExternalSolution(feMetaNumber *metaNumber, feSolution *sol,
// feMesh *mesh,
//   feMetaNumber *refMN, feSolution *refSol, feMesh *refMesh, const std::vector<feSpace *>
//   refSpaces)
// {
//   double normL2 = 0.0, solInt, solRef;
//   int nElm = _intSpace->getNbElm();

//   if(fabs(sol->getCurrentTime() - refSol->getCurrentTime()) > 1e-3){
//     printf("In feNorm::computeErrorNormFromExternalSolution : Warning - Solutions time differ by
//     more than 1e-3 : tRef = %f - tSol = %f.\n", refSol->getCurrentTime(), sol->getCurrentTime());
//     printf("In feNorm::computeErrorNormFromExternalSolution : Warning - Maybe the solution file
//     does not match the current solution.\n");
//   }

//   std::vector<double> x(3, 0.0);
//   std::vector<double> r(3, 0.0);
//   std::vector<double> J = mesh->getCncGeoByTag(cncGeoTag)->getJacobians();

//   // Identify the feSpace from among the set of spaces from the reference solution (we should use
//   a map)
//   // Aussi : la numérotation doit être recréée pour assigner correctement les ddl aux espaces.
//   // Donc soit elle est la même que pour le problème courant, ce qui suppose que les fespaces
//   sont les memes, ce qui a du sens, soit il faut redéfinir
//   // un nouveau set de feSpace à partir d'une autre simulation, ce qui semble hasardeux.
//   bool matchingSpace = false;
//   for(feSpace *fS : refSpaces){
//     if(fS->getFieldID() == _intSpace->getFieldID() && fS->getCncGeoID() ==
//     _intSpace->getCncGeoID()){
//       matchingSpace = true;

//       feNumber *n = refMN->getNumbering(fS->getFieldID());

//       for(int iElm = 0; iElm < nElm; ++iElm) {
//         _intSpace->initializeAddressingVector(metaNumber->getNumbering(_intSpace->getFieldID()),
//         iElm); _intSpace->initializeSolution(sol); geoCoord = mesh->getCoord(cncGeoTag, iElm);
//         for(int k = 0; k < _nQuad; ++k) {
//           geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
//           // Current solution
//           solInt = _intSpace->interpolateSolutionAtQuadNode(k);
//           // Reference solution evaluated at quad node
//           int elm = -1;
//           bool isFound = static_cast<feMesh2DP1*>(refMesh)->locateVertex(x,elm,r);
//           if(!isFound){
//             printf("In feNorm::computeErrorNormFromExternalSolution : Warning - Point (%f, %f,
//             %f) was not found in the mesh.\n", x[0], x[1], x[2]); solRef = solInt; // Points
//             outside the mesh do not contribute
//           } else{
//             fS->initializeAddressingVector(n, elm);
//             fS->initializeSolution(refSol);
//             solRef = fS->interpolateSolution(r.data());
//             normL2 += (solInt - solRef) * (solInt - solRef) * J[_nQuad * iElm + k] * w[k];
//           }

//           // normL2 += (solInt - solRef) * (solInt - solRef) * J[_nQuad * iElm + k] * w[k];
//         }
//       }

//       norm = sqrt(normL2);

//       break;
//     }
//   }

//   if(!matchingSpace)
//     printf("In feNorm::computeErrorNormFromExternalSolution : Error - No finite element space in
//     the reference set matches the target space on this connectivity.\n");
// }