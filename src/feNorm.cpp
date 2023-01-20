#include "feNorm.h"
#include <cmath>

feNorm::feNorm(const std::vector<feSpace*> &spaces, feMesh *mesh, int degQuad, feFunction *scalarSolution,
               feVectorFunction *vectorSolution)
  : _spaces(spaces)
  , _cnc(spaces[0]->getCncGeo())
  , _nNodePerElem(spaces[0]->getNbNodePerElem())
  , _scalarSolution(scalarSolution)
  , _vectorSolution(vectorSolution)
  , _J(_cnc->getJacobians())
{
  _nQuad = spaces[0]->getNbQuadPoints();
  _w = spaces[0]->getQuadratureWeights();
  _geoSpace = _cnc->getFeSpace();
  _nElm = _spaces[0]->getNbElm();
  _pos.resize(3);
  _geoCoord.resize(3 * _cnc->getNbNodePerElem());

  _adr.resize(spaces.size());
  for(size_t k = 0; k < spaces.size(); ++k)
    _adr[k].resize(spaces[k]->getNbFunctions());
  _sol.resize(spaces.size());
  for(size_t k = 0; k < spaces.size(); ++k)
    _sol[k].resize(spaces[k]->getNbFunctions());
}

void feNorm::computeL2Norm_uh(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh, double scaling)
{
  double normL2 = 0.0, solInt;

  std::vector<feInt> adr0;
  std::vector<double> sol0;
  std::vector<double> &solVec = sol->getSolutionReference();

  #if defined(HAVE_OMP)
  #pragma omp parallel for private(_geoCoord, adr0, sol0, solInt) reduction(+ : normL2) schedule(dynamic)
  #endif
  for(int iElm = 0; iElm < _nElm; ++iElm) {
    adr0.resize(_spaces[0]->getNbFunctions());
    _spaces[0]->initializeAddressingVector(iElm, adr0);
    sol0.resize(adr0.size());
    for(size_t i = 0; i < adr0.size(); ++i) {
      sol0[i] = solVec[adr0[i]];
    }
    
    mesh->getCoord(_cnc, iElm, _geoCoord);

    for(int k = 0; k < _nQuad; ++k) {
      solInt = _spaces[0]->interpolateFieldAtQuadNode(sol0, k);

      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);

      normL2 += (solInt * scaling) * (solInt * scaling) * _J[_nQuad * iElm + k] * _w[k];
    }
  }
  _norm = sqrt(normL2);
}

void feNorm::computeL2Norm0D(feSolution *sol)
{
  double normL2 = 0.0, solInt, solRef, t = sol->getCurrentTime();

  std::vector<double> x;
  solRef = (_scalarSolution != nullptr) ? _scalarSolution->eval(t, x) : 0.0;
  solInt = sol->getSolAtDOF(0);
  normL2 += (abs(solInt) - abs(solRef)) * (abs(solInt) - abs(solRef));
  _norm = sqrt(normL2);
}

// This is actually the L2 norm of the error u-uh, with u = solRef
void feNorm::computeL2Norm(feSolution *sol, feMesh *mesh)
{
  double res = 0.0, uh, u, t = sol->getCurrentTime();

  #if defined(HAVE_OMP)
  #pragma omp parallel for private(_geoCoord, uh, u) reduction(+ : res) schedule(dynamic)
  #endif

  for(int iElm = 0; iElm < _nElm; ++iElm) {

    _spaces[0]->initializeAddressingVector(iElm, _adr[0]);
    sol->getSolAtDOF(_adr[0], _sol[0]);
    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for(int k = 0; k < _nQuad; ++k) {
      uh = _spaces[0]->interpolateFieldAtQuadNode(_sol[0], k);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);
      u = (_scalarSolution != nullptr) ? _scalarSolution->eval(t, _pos) : 0.0;
      res += (u - uh) * (u - uh) * _J[_nQuad * iElm + k] * _w[k];
    }
  }
  _norm = sqrt(res);
}

void feNorm::computeArea()
{
  double res = 0.0;
  for(int iElm = 0; iElm < _nElm; ++iElm) {
    for(int k = 0; k < _nQuad; ++k) {
      res += _J[_nQuad * iElm + k] * _w[k];
    }
  }
  _norm = res;
}

void feNorm::computeIntegral(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh,
                             feFunction *fun)
{
  double area = 0.0, solInt, J;
  
  std::vector<feInt> adr0(_spaces[0]->getNbFunctions());
  std::vector<double> sol0(adr0.size());
  std::vector<double> &solVec = sol->getSolutionReference();

  for(int iElm = 0; iElm < _nElm; ++iElm) {
    _spaces[0]->initializeAddressingVector(iElm, adr0);
    for(size_t i = 0; i < adr0.size(); ++i) {
      sol0[i] = solVec[adr0[i]];
    }

    mesh->getCoord(_cnc, iElm, _geoCoord);
    for(int k = 0; k < _nQuad; ++k) {
      solInt = _spaces[0]->interpolateFieldAtQuadNode(sol0, k);
      std::vector<double> x(3, 0.0);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, x);
      // J = _geoSpace->interpolateFieldAtQuadNode_rDerivative(xGeoCoord, k);
      // std::vector<double> j(3, 0.0);
      // _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(_geoCoord, k, j); // TODO : complete
      // for _dim > 1 J = j[0];
      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(_geoCoord, k, dxdr);
      _geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(_geoCoord, k, dxds);
      J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];

      // area += fun->eval(t, x) * J * _w[k];
      area += solInt * J * _w[k];
    }
  }
  _norm = area;
}

// Compute the integral of uh * f, where f is a user-defined function 
void feNorm::computeScalarProduct(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh,
                             feFunction *fun)
{
  double I = 0.0, solInt, J, t = sol->getCurrentTime();
  
  std::vector<feInt> adr0(_spaces[0]->getNbFunctions());
  std::vector<double> sol0(adr0.size());
  std::vector<double> &solVec = sol->getSolutionReference();

  
  for(int iElm = 0; iElm < _nElm; ++iElm) {
    _spaces[0]->initializeAddressingVector(iElm, adr0);
    for(size_t i = 0; i < adr0.size(); ++i) {
      sol0[i] = solVec[adr0[i]];
    }

    mesh->getCoord(_cnc, iElm, _geoCoord);
    for(int k = 0; k < _nQuad; ++k) {
      solInt = _spaces[0]->interpolateFieldAtQuadNode(sol0, k);
      std::vector<double> x(3, 0.0);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, x);
      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(_geoCoord, k, dxdr);
      _geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(_geoCoord, k, dxds);
      J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];
      I += fun->eval(t, x) * solInt * J * _w[k];
    }
  }
  _norm = I;
}

// Compute ||f - uh||_L2, where f is a user-defined function 
void feNorm::computeL2Norm_f_minus_uh(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh,
                             feFunction *fun, double scaling)
{
  double I = 0.0, solInt, J, t = sol->getCurrentTime();
  
  std::vector<feInt> adr0(_spaces[0]->getNbFunctions());
  std::vector<double> sol0(adr0.size());
  std::vector<double> &solVec = sol->getSolutionReference();

  
  for(int iElm = 0; iElm < _nElm; ++iElm) {
    _spaces[0]->initializeAddressingVector(iElm, adr0);
    for(size_t i = 0; i < adr0.size(); ++i) {
      sol0[i] = solVec[adr0[i]];
    }

    mesh->getCoord(_cnc, iElm, _geoCoord);
    for(int k = 0; k < _nQuad; ++k) {
      solInt = _spaces[0]->interpolateFieldAtQuadNode(sol0, k);
      std::vector<double> x(3, 0.0);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, x);
      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(_geoCoord, k, dxdr);
      _geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(_geoCoord, k, dxds);
      J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];
      I += (fun->eval(t, x) - solInt * scaling) * (fun->eval(t, x) - solInt * scaling) * J * _w[k];
    }
  }
  _norm = sqrt(I);
}

/* Estimates the L2 norm of the error taking an external solution as the reference solution. */
feStatus feNorm::computeErrorNormFromExternalSolution(feMetaNumber *metaNumber, feSolution *sol,
                                                      feMesh *mesh, feMetaNumber *refMN,
                                                      feSolution *refSol, feMesh *refMesh,
                                                      const std::vector<feSpace *> refSpaces)
{
  double normL2 = 0.0, solInt, solRef;
  std::vector<feInt> adr(_spaces[0]->getNbFunctions());
  std::vector<double> sol0(adr.size());
  std::vector<double> &solVec = sol->getSolutionReference();
  std::vector<double> &solVecS = refSol->getSolutionReference();

  if(fabs(sol->getCurrentTime() - refSol->getCurrentTime()) > 1e-3) {
    feWarning("Solutions time differ by  more than 1e-3 : tRef "
              "= %f - tSol = %f.\n",
              refSol->getCurrentTime(), sol->getCurrentTime());
    feWarning("Maybe the solution file does not match the current solution.\n");
  }

  std::vector<double> x(3, 0.0);
  double r[3];

  // Identify the feSpace from among the set of spaces from the reference solution (we should use a
  // map) Aussi : la numérotation doit être recréée pour assigner correctement les ddl aux espaces.
  // Donc soit elle est la même que pour le problème courant, ce qui suppose que les fespaces sont
  // les memes, ce qui a du sens, soit il faut redéfinir un nouveau set de feSpace à partir d'une
  // autre simulation, ce qui semble hasardeux.
  bool matchingSpace = false;
  for(feSpace *fS : refSpaces) {
    std::vector<feInt> adrS(fS->getNbFunctions());
    std::vector<double> solS(adrS.size());
    if(fS->getFieldID() == _spaces[0]->getFieldID() &&
       fS->getCncGeoID() == _spaces[0]->getCncGeoID()) {
      matchingSpace = true;

      // feNumber *n = refMN->getNumbering(fS->getFieldID());

      for(int iElm = 0; iElm < _nElm; ++iElm) {
        _spaces[0]->initializeAddressingVector(iElm, adr);
        for(size_t i = 0; i < adr.size(); ++i) {
          sol0[i] = solVec[adr[i]];
        }

        mesh->getCoord(_cnc, iElm, _geoCoord);
        for(int k = 0; k < _nQuad; ++k) {
          _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, x);
          // Current solution
          solInt = _spaces[0]->interpolateFieldAtQuadNode(sol0, k);
          // Reference solution evaluated at quad node
          int elm = -1;
          bool isFound = static_cast<feMesh2DP1 *>(refMesh)->locateVertex(x.data(), elm, r);
          if(!isFound) {
            feWarning("Point (%f, %f, %f) was not found in the mesh.\n", x[0], x[1], x[2]);
            solRef = solInt; // Points outside the mesh do not contribute
          } else {
            fS->initializeAddressingVector(elm, adrS);
            for(size_t i = 0; i < adr.size(); ++i) {
              solS[i] = solVecS[adrS[i]];
            }
            solRef = fS->interpolateField(solS, r);
            normL2 += (solInt - solRef) * (solInt - solRef) * _J[_nQuad * iElm + k] * _w[k];
          }

          // normL2 += (solInt - solRef) * (solInt - solRef) * J[_nQuad * iElm + k] * _w[k];
        }
      }

      _norm = sqrt(normL2);

      break;
    }
  }

  if(!matchingSpace)
    return feErrorMsg(FE_STATUS_ERROR, "No finite element space in the reference set "
                                       "matches the target space on this connectivity.\n");

  return FE_STATUS_OK;
}
